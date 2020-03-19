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
  subroutine forward_euler(tree, dt, s_dt, s_prev, s_out, set_dt)
    use m_chemistry
    use m_streamer
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt     !< Time step
    !> Time state to compute derivatives from
    integer, intent(in)       :: s_dt
    !> Time state to add derivatives to
    integer, intent(in)       :: s_prev
    !< Time state to store result in
    integer, intent(in)       :: s_out
    logical, intent(in)       :: set_dt !< Whether to set new time step
    integer                   :: lvl, i, id, p_id, nc

    nc = tree%n_cell

    ! First calculate fluxes
    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call fluxes_elec(tree, id, nc, dt, s_dt, set_dt)
       end do
       !$omp end do
    end do
    !$omp end parallel

    call af_consistent_fluxes(tree, [flux_elec])

    ! Update the solution
    !$omp parallel private(lvl, i, id, p_id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call update_solution(tree%boxes(id), nc, dt, s_dt, &
               s_prev, s_out, set_dt)
       end do
       !$omp end do
    end do
    !$omp end parallel
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
    real(dp)                   :: mu, fld_face, drt_fac, tmp
    real(dp)                   :: nsmall, N_inv
    real(dp)                   :: dt_cfl, dt_drt, dt_dif
    real(dp)                   :: vmean(NDIM)
    integer                    :: n, m, IJK, tid
#if NDIM == 3
    integer                    :: l
#endif

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
#if NDIM == 2
          vmean = 0.5_dp * (v(IJK, :) + &
               [v(i+1, j, 1), v(i, j+1, 2)])
#elif NDIM == 3
          vmean = 0.5_dp * (v(IJK, :) + &
               [v(i+1, j, k, 1), v(i, j+1, k, 2), v(i, j, k+1, 3)])
#endif
          ! CFL condition
          dt_cfl = 1.0_dp/sum(abs(vmean) * inv_dr)

          ! Diffusion condition
          dt_dif = minval(tree%boxes(id)%dr)**2 / &
               max(2 * NDIM * maxval(dc(IJK, :)), epsilon(1.0_dp))

          ! Take the combined CFL-diffusion condition with Courant number 0.5
          dt_cfl = 0.5_dp/(1/dt_cfl + 1/dt_dif)

          if (ST_drt_limit_flux) then
             mu = ST_ion_mobility
          else
             if (gas_constant_density) then
                N_inv = 1 / gas_number_density
             else
                N_inv = 1 / tree%boxes(id)%cc(IJK, i_gas_dens)
             end if

             Td = tree%boxes(id)%cc(IJK, i_electric_fld) * SI_to_Townsend * N_inv
             mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
             mu = max(mu, ST_ion_mobility)
          end if

          ! Dielectric relaxation time
          dt_drt = UC_eps0 / max(UC_elem_charge * mu * cc(IJK, 1), epsilon(1.0_dp))

          dt_matrix(dt_ix_drt, tid) = min(dt_matrix(dt_ix_drt, tid), dt_drt)
          dt_matrix(dt_ix_cfl, tid) = min(dt_matrix(dt_ix_cfl, tid), dt_cfl)
          dt_matrix(dt_ix_diff, tid) = min(dt_matrix(dt_ix_diff, tid), dt_dif)
       end do; CLOSE_DO
    end if

#if NDIM == 2
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
#if NDIM == 2
    real(dp)                   :: rfac(2, box%n_cell)
#endif
    integer                    :: IJK, ix, n_cells, n, iv
    integer                    :: tid
    real(dp), parameter        :: eps = 1e-100_dp

    n_cells = box%n_cell**NDIM
    inv_dr  = 1/box%dr

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

       ! Contribution of flux
#if NDIM == 2
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

       if (ST_use_dielectric) then
          if (box%cc(IJK, i_eps) > 1.0_dp) then
             ! Convert electrons to 'minus' positive ions, so they remain stuck
             derivs(ix, ix_1pos_ion) = derivs(ix, ix_1pos_ion) - &
                  derivs(ix, ix_electron)
             derivs(ix, ix_electron) = 0.0_dp
          end if
       end if

       do n = n_gas_species+1, n_species
          iv = species_itree(n)
          box%cc(IJK, iv+s_out) = box%cc(IJK, iv+s_prev) + dt * derivs(ix, n)
       end do
    end do; CLOSE_DO

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
#if NDIM == 2
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

end module m_fluid_lfa
