#include "../afivo/src/cpp_macros.h"
!> Fluid model module
module m_fluid_lfa
  use m_af_all

  implicit none
  private

  ! Public methods
  public :: forward_euler
  public :: fluxes_elec
  public :: update_solution

contains

  !> Advance fluid model using forward Euler step. If the equation is written as
  !> y' = f(y), the result is: y(s_out) = y(s_prev) + f(y(s_dt)), where the
  !> s_... refer to temporal states.
  subroutine forward_euler(tree, dt, dt_lim, time, s_deriv, s_prev, s_out, &
       i_step, n_steps)
    use m_chemistry
    use m_streamer
    use m_field
    use m_dt
    use m_transport_data
    use m_dielectric
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt     !< Time step
    real(dp), intent(inout)   :: dt_lim !< Time step limitation
    real(dp), intent(in)      :: time
    !> Time state to compute derivatives from
    integer, intent(in)       :: s_deriv
    !> Time state to add derivatives to
    integer, intent(in)       :: s_prev
    !< Time state to store result in
    integer, intent(in)       :: s_out
    integer, intent(in)       :: i_step, n_steps
    integer                   :: lvl, i, id, p_id, nc, ix, id_out
    logical                   :: set_dt

    nc = tree%n_cell

    set_dt = (i_step == n_steps)

    ! Use a shared array to determine maximum time step
    dt_matrix(1:dt_num_cond, :) = dt_max

    ! So that ghost cells can be computed properly near refinement boundaries
    call af_restrict_ref_boundary(tree, flux_species+s_deriv)

    ! Since field_compute is called after performing time integration, we don't
    ! have to call it again for the first sub-step of the next iteration
    if (i_step > 1) call field_compute(tree, mg, s_deriv, time, .true.)

    ! First calculate fluxes
    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call fluxes_elec(tree, id, nc, dt, s_deriv, set_dt)

          if (transport_data_ions%n_mobile_ions > 0) then
             call fluxes_ions(tree, id, nc, dt, s_deriv, set_dt)
          end if
       end do
       !$omp end do
    end do
    !$omp end parallel

    call af_consistent_fluxes(tree, flux_variables)

    if (transport_data_ions%n_mobile_ions > 0 .and. &
         ion_se_yield > 0.0_dp) then
       ! Handle secondary electron emission from ions
       call af_loop_box(tree, handle_ion_se_flux, .true.)
    end if

    ! Update the solution
    !$omp parallel private(lvl, i, id, p_id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call update_solution(tree%boxes(id), nc, dt, s_deriv, &
               s_prev, s_out, set_dt)

       end do
       !$omp end do
    end do
    !$omp end parallel

    if (ST_use_dielectric) then
       ! Update surface charge and handle photon emission
       ! @todo For parallelization, think about corner cells with two surfaces
       do ix = 1, diel%max_ix
          if (diel%surfaces(ix)%in_use) then
             id_out = diel%surfaces(ix)%id_out

             ! Convert fluxes onto dielectric to surface charge, and handle
             ! secondary emission
             call dielectric_update_surface_charge(tree%boxes(id_out), &
                  diel%surfaces(ix), dt, s_deriv, s_out)

             ! Add secondary emission from photons hitting the surface
             call dielectric_photon_emission(tree%boxes(id_out), &
                  diel%surfaces(ix), dt, s_deriv, s_out)
          end if
       end do
    end if

    if (set_dt) then
       dt_lim = min(2 * global_dt, dt_safety_factor * &
            minval(dt_matrix(1:dt_num_cond, :)))
    end if
  end subroutine forward_euler

  !> Compute the electron fluxes due to drift and diffusion
  subroutine fluxes_elec(tree, id, nc, dt, s_in, set_dt)
    use m_af_flux_schemes
    use m_units_constants
    use omp_lib
    use m_streamer
    use m_gas
    use m_dt
    use m_lookup_table
    use m_transport_data
    type(af_t), intent(inout) :: tree
    integer, intent(in)        :: id
    integer, intent(in)        :: nc   !< Number of cells per dimension
    real(dp), intent(in)       :: dt
    integer, intent(in)        :: s_in !< Input time state
    logical, intent(in)        :: set_dt
    real(dp)                   :: inv_dr(NDIM), fld, Td
    ! Velocities at cell faces
    real(dp)                   :: v(DTIMES(1:nc+1), NDIM)
    ! Diffusion coefficients at cell faces
    real(dp)                   :: dc(DTIMES(1:nc+1), NDIM)
    ! Maximal fluxes at cell faces
    real(dp)                   :: fmax(DTIMES(1:nc+1), NDIM)
    ! Cell-centered densities
    real(dp)                   :: cc(DTIMES(-1:nc+2), 1)
    real(dp)                   :: mu, max_mu_ion, fld_face, drt_fac, tmp
    real(dp)                   :: nsmall, N_inv
    real(dp)                   :: dt_cfl, dt_drt, dt_dif
    real(dp)                   :: vmean(NDIM)
    real(dp), parameter        :: eps = 1e-100_dp
    integer                    :: n, IJK, tid
#if NDIM > 1
    integer                    :: m
#endif
#if NDIM == 3
    integer                    :: l
#endif

    ! Inside the dielectric, set the flux to zero. We later determine the
    ! boundary flux onto dielectrics
    if (ST_use_dielectric) then
       if (tree%boxes(id)%cc(DTIMES(1), i_eps) > 1.0_dp) then
          tree%boxes(id)%fc(DTIMES(:), :, flux_elec) = 0.0_dp
          return
       end if
    end if

    inv_dr  = 1/tree%boxes(id)%dr
    drt_fac = UC_eps0 / max(1e-100_dp, UC_elem_charge * dt)
    nsmall  = 1.0_dp ! A small density
    N_inv   = 1/gas_number_density ! For constant gas densities

    v  = 0.0_dp
    dc = 0.0_dp

    if (ST_drt_limit_flux) then
       fmax = 0.0_dp
    end if

    ! Fill cc with interior data plus two layers of ghost cells
    call af_gc2_box(tree, id, [i_electron+s_in], cc)

    associate(box => tree%boxes(id))
      ! We use the average field to compute the mobility and diffusion coefficient
      ! at the interface
      do n = 1, nc+1
#if NDIM == 1
         if (.not. gas_constant_density) then
            N_inv = 2 / (box%cc(n-1, i_gas_dens) + &
                 box%cc(n, i_gas_dens))
         end if

         fld = 0.5_dp * (box%cc(n-1, i_electric_fld) + &
              box%cc(n, i_electric_fld))
         Td = fld * SI_to_Townsend * N_inv
         mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
         fld_face = box%fc(n, 1, electric_fld)
         v(n, 1)  = -mu * fld_face
         dc(n, 1) = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

         if (ST_drt_limit_flux) then
            tmp = abs(cc(n-1, 1) - cc(n, 1)) / &
                 max(cc(n-1, 1), cc(n, 1), nsmall)
            tmp = max(fld, tmp * inv_dr(1) * dc(n, 1) / mu)
            fmax(n, 1) = drt_fac * tmp
         end if

         if (abs(fld_face) > ST_diffusion_field_limit .and. &
              (cc(n-1, 1) - cc(n, 1)) * fld_face > 0.0_dp) then
            dc(n, 1) = 0.0_dp
         end if
#else
         do m = 1, nc
#if NDIM == 2
            if (.not. gas_constant_density) then
               N_inv = 2 / (box%cc(n-1, m, i_gas_dens) + &
                    box%cc(n, m, i_gas_dens))
            end if

            fld = 0.5_dp * (box%cc(n-1, m, i_electric_fld) + &
                 box%cc(n, m, i_electric_fld))
            Td = fld * SI_to_Townsend * N_inv
            mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
            fld_face = box%fc(n, m, 1, electric_fld)
            v(n, m, 1)  = -mu * fld_face
            dc(n, m, 1) = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

            if (ST_drt_limit_flux) then
               tmp = abs(cc(n-1, m, 1) - cc(n, m, 1)) / &
                    max(cc(n-1, m, 1), cc(n, m, 1), nsmall)
               tmp = max(fld, tmp * inv_dr(1) * dc(n, m, 1) / mu)
               fmax(n, m, 1) = drt_fac * tmp
            end if

            if (abs(fld_face) > ST_diffusion_field_limit .and. &
                 (cc(n-1, m, 1) - cc(n, m, 1)) * fld_face > 0.0_dp) then
               dc(n, m, 1) = 0.0_dp
            end if

            if (.not. gas_constant_density) then
               N_inv = 2 / (box%cc(m, n-1, i_gas_dens) + &
                    box%cc(m, n, i_gas_dens))
            end if

            fld = 0.5_dp * (box%cc(m, n-1, i_electric_fld) + &
                 box%cc(m, n, i_electric_fld))
            Td = fld * SI_to_Townsend * N_inv
            mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
            fld_face = box%fc(m, n, 2, electric_fld)
            v(m, n, 2)  = -mu * fld_face
            dc(m, n, 2) = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

            if (ST_drt_limit_flux) then
               tmp = abs(cc(m, n-1, 1) - cc(m, n, 1)) / &
                    max(cc(m, n-1, 1), cc(m, n, 1), nsmall)
               tmp = max(fld, tmp * inv_dr(2) * dc(m, n, 2) / mu)
               fmax(m, n, 2) = drt_fac * tmp
            end if

            if (abs(fld_face) > ST_diffusion_field_limit .and. &
                 (cc(m, n-1, 1) - cc(m, n, 1)) * fld_face > 0.0_dp) then
               dc(m, n, 2) = 0.0_dp
            end if
#elif NDIM == 3
            do l = 1, nc
               if (.not. gas_constant_density) then
                  N_inv = 2 / (box%cc(n-1, m, l, i_gas_dens) + &
                       box%cc(n, m, l, i_gas_dens))
               end if

               fld = 0.5_dp * (&
                    box%cc(n-1, m, l, i_electric_fld) + &
                    box%cc(n, m, l, i_electric_fld))
               Td = fld * SI_to_Townsend * N_inv
               mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
               fld_face = box%fc(n, m, l, 1, electric_fld)
               v(n, m, l, 1)  = -mu * fld_face
               dc(n, m, l, 1) = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

               if (ST_drt_limit_flux) then
                  tmp = abs(cc(n-1, m, l, 1) - cc(n, m, l, 1)) / &
                       max(cc(n-1, m, l, 1), cc(n, m, l, 1), nsmall)
                  tmp = max(fld, tmp * inv_dr(1) * dc(n, m, l, 1) / mu)
                  fmax(n, m, l, 1) = drt_fac * tmp
               end if

               if (abs(fld_face) > ST_diffusion_field_limit .and. &
                    (cc(n-1, m, l, 1) - cc(n, m, l, 1)) * fld_face > 0.0_dp) then
                  dc(n, m, l, 1) = 0.0_dp
               end if

               if (.not. gas_constant_density) then
                  N_inv = 2 / (box%cc(m, n-1, l, i_gas_dens) + &
                       box%cc(m, n, l, i_gas_dens))
               end if

               fld = 0.5_dp * (&
                    box%cc(m, n-1, l, i_electric_fld) + &
                    box%cc(m, n, l, i_electric_fld))
               Td = fld * SI_to_Townsend * N_inv
               mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
               fld_face = box%fc(m, n, l, 2, electric_fld)
               v(m, n, l, 2)  = -mu * fld_face
               dc(m, n, l, 2) = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

               if (ST_drt_limit_flux) then
                  tmp = abs(cc(m, n-1, l, 1) - cc(m, n, l, 1)) / &
                       max(cc(m, n-1, l, 1), cc(m, n, l, 1), nsmall)
                  tmp = max(fld, tmp * inv_dr(2) * dc(m, n, l, 2) / mu)
                  fmax(m, n, l, 2) = drt_fac * tmp
               end if

               if (abs(fld_face) > ST_diffusion_field_limit .and. &
                    (cc(m, n-1, l, 1) - cc(m, n, l, 1)) * fld_face > 0.0_dp) then
                  dc(m, n, l, 2) = 0.0_dp
               end if

               if (.not. gas_constant_density) then
                  N_inv = 2 / (box%cc(m, l, n-1, i_gas_dens) + &
                       box%cc(m, l, n, i_gas_dens))
               end if

               fld = 0.5_dp * (&
                    box%cc(m, l, n-1, i_electric_fld) + &
                    box%cc(m, l, n, i_electric_fld))
               Td = fld * SI_to_Townsend * N_inv
               mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
               fld_face = box%fc(m, l, n, 3, electric_fld)
               v(m, l, n, 3)  = -mu * fld_face
               dc(m, l, n, 3) = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

               if (ST_drt_limit_flux) then
                  tmp = abs(cc(m, l, n-1, 1) - cc(m, l, n, 1)) / &
                       max(cc(m, l, n-1, 1), cc(m, l, n, 1), nsmall)
                  tmp = max(fld, tmp * inv_dr(3) * dc(m, l, n, 3) / mu)
                  fmax(m, l, n, 3) = drt_fac * tmp
               end if

               if (abs(fld_face) > ST_diffusion_field_limit .and. &
                    (cc(m, l, n-1, 1) - cc(m, l, n, 1)) * fld_face > 0.0_dp) then
                  dc(m, l, n, 3) = 0.0_dp
               end if
            end do
#endif
         end do
#endif
      end do
    end associate

    if (ST_max_velocity > 0.0_dp) then
       where (abs(v) > ST_max_velocity)
          v = sign(ST_max_velocity, v)
       end where
    end if

    if (set_dt) then
       tid    = omp_get_thread_num() + 1
       dt_cfl = dt_max
       dt_drt = dt_matrix(dt_ix_cfl, tid)

       do KJI_DO(1,nc)
#if NDIM == 1
          vmean = 0.5_dp * (v(IJK, :) + [v(i+1, 1)])
#elif NDIM == 2
          vmean = 0.5_dp * (v(IJK, :) + &
               [v(i+1, j, 1), v(i, j+1, 2)])
#elif NDIM == 3
          vmean = 0.5_dp * (v(IJK, :) + &
               [v(i+1, j, k, 1), v(i, j+1, k, 2), v(i, j, k+1, 3)])
#endif
          ! CFL condition
          dt_cfl = 1.0_dp/sum(max(abs(vmean), eps) * inv_dr)

          ! Diffusion condition
          dt_dif = minval(tree%boxes(id)%dr)**2 / &
               max(2 * NDIM * maxval(dc(IJK, :)), eps)

          ! Take the combined CFL-diffusion condition with Courant number 0.5
          dt_cfl = 0.5_dp/(1/dt_cfl + 1/dt_dif)

          if (gas_constant_density) then
             N_inv = 1 / gas_number_density
          else
             N_inv = 1 / tree%boxes(id)%cc(IJK, i_gas_dens)
          end if

          ! Ion mobility
          if (size(transport_data_ions%mobilities) > 0) then
             max_mu_ion = maxval(abs(transport_data_ions%mobilities)) * N_inv
          else
             max_mu_ion = 0.0_dp
          end if

          ! Electron mobility
          if (ST_drt_limit_flux) then
             mu = 0.0_dp
          else
             Td = tree%boxes(id)%cc(IJK, i_electric_fld) * SI_to_Townsend * N_inv
             mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
          end if

          ! Take maximum of electron and ion mobility, in the future we should
          ! probably apply a weighted sum (weighted with species densities)
          mu = max(mu, max_mu_ion)

          ! Dielectric relaxation time
          dt_drt = UC_eps0 / max(UC_elem_charge * mu * cc(IJK, 1), eps)

          dt_matrix(dt_ix_drt, tid) = min(dt_matrix(dt_ix_drt, tid), dt_drt)
          dt_matrix(dt_ix_cfl, tid) = min(dt_matrix(dt_ix_cfl, tid), dt_cfl)
          dt_matrix(dt_ix_diff, tid) = min(dt_matrix(dt_ix_diff, tid), dt_dif)
       end do; CLOSE_DO
    end if

#if NDIM == 1
    call flux_koren_1d(cc, v, nc, 2)
    call flux_diff_1d(cc, dc, inv_dr(1), nc, 2)
#elif NDIM == 2
    call flux_koren_2d(cc, v, nc, 2)
    call flux_diff_2d(cc, dc, inv_dr, nc, 2)
#elif NDIM == 3
    call flux_koren_3d(cc, v, nc, 2)
    call flux_diff_3d(cc, dc, inv_dr, nc, 2)
#endif

    tree%boxes(id)%fc(DTIMES(:), :, flux_elec) = v + dc

    if (ST_drt_limit_flux) then
       where (abs(tree%boxes(id)%fc(DTIMES(:), :, flux_elec)) > fmax)
          tree%boxes(id)%fc(DTIMES(:), :, flux_elec) = &
               sign(fmax, tree%boxes(id)%fc(DTIMES(:), :, flux_elec))
       end where
    end if

  end subroutine fluxes_elec

  subroutine fluxes_ions(tree, id, nc, dt, s_in, set_dt)
    use m_af_flux_schemes
    use m_streamer
    use m_gas
    use m_transport_data
    type(af_t), intent(inout) :: tree
    integer, intent(in)        :: id
    integer, intent(in)        :: nc   !< Number of cells per dimension
    real(dp), intent(in)       :: dt
    integer, intent(in)        :: s_in !< Input time state
    logical, intent(in)        :: set_dt
    real(dp)                   :: inv_dr(NDIM)
    ! Velocities at cell faces
    real(dp)                   :: v(DTIMES(1:nc+1), NDIM)
    ! Cell-centered densities
    real(dp)                   :: cc(DTIMES(-1:nc+2), &
         transport_data_ions%n_mobile_ions)
    real(dp)                   :: mu
    integer                    :: n, i_ion, i_flux, ix
#if NDIM > 1
    integer                    :: m
#endif
#if NDIM == 3
    integer                    :: l
#endif

    inv_dr  = 1/tree%boxes(id)%dr

    call af_gc2_box(tree, id, [flux_species(2:)+s_in], cc)

    do ix = 1, transport_data_ions%n_mobile_ions
       i_ion  = flux_species(ix+1)
       i_flux = flux_variables(ix+1)
       mu     = transport_data_ions%mobilities(ix)

#if NDIM == 1
       do n = 1, nc+1
          v(n, 1) = mu * get_N_inv_face(tree%boxes(id), n, 1) * &
               tree%boxes(id)%fc(n, 1, electric_fld)
       end do
#elif NDIM == 2
       do n = 1, nc+1
          do m = 1, nc
             v(n, m, 1) = mu * get_N_inv_face(tree%boxes(id), n, m, 1) * &
                  tree%boxes(id)%fc(n, m, 1, electric_fld)
             v(m, n, 2) = mu * get_N_inv_face(tree%boxes(id), m, n, 2) * &
                  tree%boxes(id)%fc(m, n, 2, electric_fld)
          end do
       end do
#elif NDIM == 3
       do n = 1, nc+1
          do m = 1, nc
             do l = 1, nc
                v(n, m, l, 1) = mu * get_N_inv_face(tree%boxes(id), n, m, l, 1) * &
                     tree%boxes(id)%fc(n, m, l, 1, electric_fld)
                v(m, n, l, 2) = mu * get_N_inv_face(tree%boxes(id), m, n, l, 2) * &
                     tree%boxes(id)%fc(m, n, l, 2, electric_fld)
                v(m, l, n, 3) = mu * get_N_inv_face(tree%boxes(id), m, l, n, 3) * &
                     tree%boxes(id)%fc(m, l, n, 3, electric_fld)
             end do
          end do
       end do
#endif

#if NDIM == 1
       call flux_koren_1d(cc(DTIMES(:), ix), v, nc, 2)
#elif NDIM == 2
       call flux_koren_2d(cc(DTIMES(:), ix), v, nc, 2)
#elif NDIM == 3
       call flux_koren_3d(cc(DTIMES(:), ix), v, nc, 2)
#endif

       tree%boxes(id)%fc(DTIMES(:), :, i_flux) = v
    end do

  end subroutine fluxes_ions

  !> Get inverse gas density at a cell face, between cell-centered index i-1 and
  !> i along dimension idim
  pure real(dp) function get_N_inv_face(box, IJK, idim)
    use m_gas
    type(box_t), intent(in) :: box
    integer, intent(in)     :: IJK
    integer, intent(in)     :: idim !< Direction of flux through cell face

    if (gas_constant_density) then
       get_N_inv_face = gas_inverse_number_density
    else
#if NDIM == 1
       get_N_inv_face = 2 / (box%cc(i-1, i_gas_dens) + &
            box%cc(i, i_gas_dens))
#elif NDIM == 2
       select case (idim)
       case (1)
          get_N_inv_face = 2 / (box%cc(i-1, j, i_gas_dens) + &
            box%cc(i, j, i_gas_dens))
       case default
          get_N_inv_face = 2 / (box%cc(i, j-1, i_gas_dens) + &
            box%cc(i, j, i_gas_dens))
       end select
#elif NDIM == 3
       select case (idim)
       case (1)
          get_N_inv_face = 2 / (box%cc(i-1, j, k, i_gas_dens) + &
            box%cc(i, j, k, i_gas_dens))
       case (2)
          get_N_inv_face = 2 / (box%cc(i, j-1, k, i_gas_dens) + &
               box%cc(i, j, k, i_gas_dens))
       case default
          get_N_inv_face = 2 / (box%cc(i, j, k-1, i_gas_dens) + &
               box%cc(i, j, k, i_gas_dens))
       end select
#endif
    end if
  end function get_N_inv_face

  !> Advance solution in a box over dt based on the fluxes and reactions, using
  !> a forward Euler update
  subroutine update_solution(box, nc, dt, s_dt, s_prev, s_out, set_dt)
    use omp_lib
    use m_units_constants
    use m_gas
    use m_chemistry
    use m_streamer
    use m_photoi
    use m_dt
    use m_lookup_table
    use m_transport_data
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: nc     !< Box size
    real(dp), intent(in)       :: dt     !< Time step
    integer, intent(in)        :: s_dt   !< Time state to compute derivatives from
    integer, intent(in)        :: s_prev !< Time state to add derivatives to
    integer, intent(in)        :: s_out  !< Output time state
    logical, intent(in)        :: set_dt !< Whether to set new time step
    real(dp)                   :: inv_dr(NDIM)
    real(dp)                   :: tmp
    real(dp)                   :: rates(nc**NDIM, n_reactions)
    real(dp)                   :: derivs(nc**NDIM, n_species)
    real(dp)                   :: dens(nc**NDIM, n_species)
    real(dp)                   :: fields(nc**NDIM)
    real(dp)                   :: source_factor(nc**NDIM)
    real(dp)                   :: coords(nc, NDIM), r(NDIM)
#if NDIM == 2
    real(dp)                   :: rfac(2, box%n_cell)
#endif
    integer                    :: IJK, ix, n_cells, n, iv, i_flux
    integer                    :: tid
    real(dp), parameter        :: eps = 1e-100_dp

    n_cells = box%n_cell**NDIM
    inv_dr  = 1/box%dr
    do n = 1, NDIM
       coords(:, n) = box%r_min(n) + box%dr(n) * [(i-0.5_dp, i=1,nc)]
    end do

    ! Inside the dielectric, do not update the species densities, which should
    ! always be zero
    if (ST_use_dielectric) then
       if (box%cc(DTIMES(1), i_eps) > 1.0_dp) then
          return
       end if
    end if

    if (gas_constant_density) then
       ! Compute field in Townsends
       tmp = 1 / gas_number_density
       fields = SI_to_Townsend * tmp * &
            pack(box%cc(DTIMES(1:nc), i_electric_fld), .true.)
    else
       do n = 1, n_gas_species
          dens(:, n) = gas_fractions(n) * &
               pack(box%cc(DTIMES(1:nc), i_gas_dens), .true.)
       end do

       fields(:) = SI_to_Townsend * pack( &
            box%cc(DTIMES(1:nc), i_electric_fld) / &
            box%cc(DTIMES(1:nc), i_gas_dens), .true.)
    end if

    dens(:, n_gas_species+1:n_species) = reshape(box%cc(DTIMES(1:nc), &
         species_itree(n_gas_species+1:n_species)+s_dt), [n_cells, n_plasma_species])

    call get_rates(fields, rates, n_cells)

    if (ST_source_factor .or. ST_source_min_density > 0) then
       call compute_source_factor(box, nc, dens(:, ix_electron), &
            fields, source_factor)

       do n = 1, n_reactions
          if (reactions(n)%reaction_type == ionization_reaction) then
             rates(:, n) = rates(:, n) * source_factor
          end if
       end do
    end if

    call get_derivatives(dens, rates, derivs, n_cells)

    ! TODO clean this up in the future (when we set electrode Neumann boundary
    ! conditions in actual flux computation)
    if (box%tag == mg_lsf_box) then
       ix = 0
       do KJI_DO(1,nc)
          ix = ix + 1
          if (box%cc(IJK, i_lsf) <= 0.0_dp) then
             derivs(ix, :) = 0
          end if
       end do; CLOSE_DO
    end if

    if (set_dt) then
       tid = omp_get_thread_num() + 1
       dt_matrix(dt_ix_rates, tid) = min(dt_matrix(dt_ix_rates, tid), &
            minval((abs(dens) + dt_chemistry_nmin) / max(abs(derivs), eps)))
    end if

#if NDIM == 2
    if (ST_cylindrical) then
       call af_cyl_flux_factors(box, rfac)
    else
       rfac = 1.0_dp
    end if
#endif

    ix = 0
    do KJI_DO(1,nc)
       ix = ix + 1

       if (ST_plasma_region_enabled) then
          r(1) = coords(i, 1)
#if NDIM > 1
          r(2) = coords(j, 2)
#endif
#if NDIM > 2
          r(3) = coords(k, 3)
#endif

          if (any(r < ST_plasma_region_rmin) .or. &
               any(r > ST_plasma_region_rmax)) then
             ! Disable all reactions
             derivs(ix, :) = 0.0_dp
          end if
       end if

       ! Contribution of flux
#if NDIM == 1
       derivs(ix, ix_electron) = derivs(ix, ix_electron) + &
            inv_dr(1) * (box%fc(i, 1, flux_elec) - &
            box%fc(i+1, 1, flux_elec))
#elif NDIM == 2
       derivs(ix, ix_electron) = derivs(ix, ix_electron) + &
            inv_dr(1) * (rfac(1, i) * box%fc(i, j, 1, flux_elec) - &
            rfac(2, i) * box%fc(i+1, j, 1, flux_elec)) + &
            inv_dr(2) * (box%fc(i, j, 2, flux_elec) - &
            box%fc(i, j+1, 2, flux_elec))
#elif NDIM == 3
       derivs(ix, ix_electron) = derivs(ix, ix_electron) + &
            inv_dr(1) * (box%fc(i, j, k, 1, flux_elec) - &
            box%fc(i+1, j, k, 1, flux_elec)) + &
            inv_dr(2) * (box%fc(i, j, k, 2, flux_elec) - &
            box%fc(i, j+1, k, 2, flux_elec)) + &
            inv_dr(3) * (box%fc(i, j, k, 3, flux_elec) - &
            box%fc(i, j, k+1, 3, flux_elec))
#endif

       if (photoi_enabled) then
          derivs(ix, ix_electron) = derivs(ix, ix_electron) + &
               box%cc(IJK, i_photo)
          derivs(ix, ix_1pos_ion) = derivs(ix, ix_1pos_ion) + &
               box%cc(IJK, i_photo)
       end if

       if (ST_use_electrode) then
          if (box%cc(IJK, i_lsf) <= 0.0_dp) then
             derivs(ix, :) = 0
          end if
       end if

       do n = n_gas_species+1, n_species
          iv = species_itree(n)
          box%cc(IJK, iv+s_out) = box%cc(IJK, iv+s_prev) + dt * derivs(ix, n)
       end do
    end do; CLOSE_DO

    ! Ion fluxes
    do n = 2, size(flux_species)
       iv     = flux_species(n)
       i_flux = flux_variables(n)

       do KJI_DO(1,nc)
          ! Contribution of ion flux
          if (ST_use_electrode) then
             if (box%cc(IJK, i_lsf) <= 0.0_dp) cycle
          end if
#if NDIM == 1
          tmp = inv_dr(1) * (box%fc(i, 1, i_flux) - &
               box%fc(i+1, 1, i_flux))
#elif NDIM == 2
          tmp = inv_dr(1) * (rfac(1, i) * box%fc(i, j, 1, i_flux) - &
               rfac(2, i) * box%fc(i+1, j, 1, i_flux)) + &
               inv_dr(2) * (box%fc(i, j, 2, i_flux) - &
               box%fc(i, j+1, 2, i_flux))
#elif NDIM == 3
          tmp = inv_dr(1) * (box%fc(i, j, k, 1, i_flux) - &
               box%fc(i+1, j, k, 1, i_flux)) + &
               inv_dr(2) * (box%fc(i, j, k, 2, i_flux) - &
               box%fc(i, j+1, k, 2, i_flux)) + &
               inv_dr(3) * (box%fc(i, j, k, 3, i_flux) - &
               box%fc(i, j, k+1, 3, i_flux))
#endif
          box%cc(IJK, iv+s_out) = box%cc(IJK, iv+s_out) + tmp * dt
       end do; CLOSE_DO
    end do

  end subroutine update_solution

  !> Compute adjustment factor for electron source terms. Used to reduce them in
  !> certain regimes.
  subroutine compute_source_factor(box, nc, elec_dens, fields, source_factor)
    use m_streamer
    use m_gas
    use m_transport_data
    use m_lookup_table
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nc
    real(dp), intent(in)    :: elec_dens(nc**NDIM)
    real(dp), intent(in)    :: fields(nc**NDIM)
    real(dp), intent(out)   :: source_factor(nc**NDIM)
    real(dp)                :: tmp, mobilities(nc**NDIM)
    integer                 :: ix, IJK

    source_factor(:) = 1.0_dp

    if (ST_source_factor) then
       if (.not. gas_constant_density) &
            error stop "source_factor: gas_constant_density is false"

       tmp = 1 / gas_number_density
       mobilities = LT_get_col(td_tbl, td_mobility, fields) * tmp

       ix = 0
       do KJI_DO(1,nc)
          ix = ix + 1

          ! Compute norm of flux at cell center
#if NDIM == 1
          source_factor(ix) = 0.5_dp * norm2([ &
               box%fc(i, 1, flux_elec) + box%fc(i+1, 1, flux_elec)])
#elif NDIM == 2
          source_factor(ix) = 0.5_dp * norm2([ &
               box%fc(i, j, 1, flux_elec) + box%fc(i+1, j, 1, flux_elec), &
               box%fc(i, j, 2, flux_elec) + box%fc(i, j+1, 2, flux_elec)])
#elif NDIM == 3
          source_factor(ix) = 0.5_dp * norm2([ &
               box%fc(i, j, k, 1, flux_elec) + box%fc(i+1, j, k, 1, flux_elec), &
               box%fc(i, j, k, 2, flux_elec) + box%fc(i, j+1, k, 2, flux_elec), &
               box%fc(i, j, k, 3, flux_elec) + box%fc(i, j, k+1, 3, flux_elec)])
#endif
       end do; CLOSE_DO

       ! Compute source factor as |flux|/(n_e * mu * E)
       source_factor = source_factor / max(1e-10_dp, &
            elec_dens * mobilities * &
            pack(box%cc(DTIMES(1:nc), i_electric_fld), .true.))
       source_factor = min(1.0_dp, source_factor)
       source_factor = max(0.0_dp, source_factor)
    end if

    if (ST_source_min_density > 0) then
       ! Factor goes linearly to zero between n_min and 2 * n_min
       source_factor = min(source_factor, &
            max(0.0_dp, elec_dens-ST_source_min_density) / &
            max(elec_dens, ST_source_min_density))
    end if

  end subroutine compute_source_factor

  !> Handle secondary emission from positive ions
  subroutine handle_ion_se_flux(box)
    use m_transport_data
    use m_streamer
    type(box_t), intent(inout) :: box
    integer                    :: nc, nb, n, ion_flux

    nc = box%n_cell

    ! Return if there is no physical boundary
    if (all(box%neighbors >= af_no_box)) return

    do nb = 1, af_num_neighbors
       ! Check for physical boundary
       if (box%neighbors(nb) < af_no_box) then
          ! Loop over positive ion species
          do n = 1, transport_data_ions%n_mobile_ions
             if (flux_species_charge(n+1) > 0.0_dp) then
                ion_flux = flux_variables(n+1)
                select case (nb)
#if NDIM == 1
                case (af_neighb_lowx)
                   box%fc(1, 1, flux_elec) = box%fc(1, 1, flux_elec) - &
                        ion_se_yield * min(0.0_dp, box%fc(1, 1, ion_flux))
                case (af_neighb_highx)
                   box%fc(nc+1, 1, flux_elec) = box%fc(nc+1, 1, flux_elec) - &
                        ion_se_yield * max(0.0_dp, box%fc(1, 1, ion_flux))
#elif NDIM == 2
                case (af_neighb_lowx)
                   box%fc(1, 1:nc, 1, flux_elec) = &
                        box%fc(1, 1:nc, 1, flux_elec) - ion_se_yield * &
                        min(0.0_dp, box%fc(1, 1:nc, 1, ion_flux))
                case (af_neighb_highx)
                   box%fc(nc+1, 1:nc, 1, flux_elec) = &
                        box%fc(nc+1, 1:nc, 1, flux_elec) - ion_se_yield * &
                        max(0.0_dp, box%fc(nc+1, 1:nc, 1, ion_flux))
                case (af_neighb_lowy)
                   box%fc(1:nc, 1, 2, flux_elec) = &
                        box%fc(1:nc, 1, 2, flux_elec) - ion_se_yield * &
                        min(0.0_dp, box%fc(1:nc, 1, 2, ion_flux))
                case (af_neighb_highy)
                   box%fc(1:nc, nc+1, 2, flux_elec) = &
                        box%fc(1:nc, nc+1, 2, flux_elec) - ion_se_yield * &
                        max(0.0_dp, box%fc(1:nc, nc+1, 2, ion_flux))
#elif NDIM == 3
                case (af_neighb_lowx)
                   box%fc(1, 1:nc, 1:nc, 1, flux_elec) = &
                        box%fc(1, 1:nc, 1:nc, 1, flux_elec) - ion_se_yield * &
                        min(0.0_dp, box%fc(1, 1:nc, 1:nc, 1, ion_flux))
                case (af_neighb_highx)
                   box%fc(nc+1, 1:nc, 1:nc, 1, flux_elec) = &
                        box%fc(nc+1, 1:nc, 1:nc, 1, flux_elec) - ion_se_yield * &
                        max(0.0_dp, box%fc(nc+1, 1:nc, 1:nc, 1, ion_flux))
                case (af_neighb_lowy)
                   box%fc(1:nc, 1:nc, 1, 2, flux_elec) = &
                        box%fc(1:nc, 1:nc, 1, 2, flux_elec) - ion_se_yield * &
                        min(0.0_dp, box%fc(1:nc, 1:nc, 1, 2, ion_flux))
                case (af_neighb_highy)
                   box%fc(1:nc, nc+1, 1:nc, 2, flux_elec) = &
                        box%fc(1:nc, nc+1, 1:nc, 2, flux_elec) - ion_se_yield * &
                        max(0.0_dp, box%fc(1:nc, nc+1, 1:nc, 2, ion_flux))
                case (af_neighb_lowz)
                   box%fc(1:nc, 1:nc, 1, 3, flux_elec) = &
                        box%fc(1:nc, 1:nc, 1, 3, flux_elec) - ion_se_yield * &
                        min(0.0_dp, box%fc(1:nc, 1:nc, 1, 3, ion_flux))
                case (af_neighb_highz)
                   box%fc(1:nc, 1:nc, nc+1, 3, flux_elec) = &
                        box%fc(1:nc, 1:nc, nc+1, 3, flux_elec) - ion_se_yield * &
                        max(0.0_dp, box%fc(1:nc, 1:nc, nc+1, 3, ion_flux))
#endif

                end select
             end if
          end do
       end if
    end do

  end subroutine handle_ion_se_flux

end module m_fluid_lfa
