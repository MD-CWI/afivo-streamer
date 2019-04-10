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

  !> Advance fluid model using forward Euler step
  subroutine forward_euler(tree, dt, s_in, s_out, set_dt)
    use m_chemistry
    use m_streamer
    use m_dielectric
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt     !< Time step
    integer, intent(in)       :: s_in   !< Input time state
    integer, intent(in)       :: s_out  !< Output time state
    logical, intent(in)       :: set_dt !< Whether to set new time step
    integer                   :: lvl, i, id, p_id

    ! First calculate fluxes
    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call fluxes_elec(tree%boxes, id, dt, s_in, set_dt)
          call fluxes_ions(tree, id, s_in)
       end do
       !$omp end do
    end do
    !$omp end parallel

    call af_consistent_fluxes(tree, flux_variables)

    ! Update the solution
    !$omp parallel private(lvl, i, id, p_id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call update_solution(tree%boxes(id), id, dt, s_in, s_out, set_dt)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine forward_euler

  !> Compute the electron fluxes due to drift and diffusion
  subroutine fluxes_elec(boxes, id, dt, s_in, set_dt)
    use m_af_flux_schemes
    use m_units_constants
    use omp_lib
    use m_streamer
    use m_gas
    use m_dt
    use m_lookup_table
    use m_transport_data
    type(box_t), intent(inout) :: boxes(:)
    integer, intent(in)        :: id
    real(dp), intent(in)       :: dt
    integer, intent(in)        :: s_in   !< Input time state
    logical, intent(in)        :: set_dt
    real(dp)                   :: inv_dr(NDIM), fld, Td
    ! Velocities at cell faces
    real(dp), allocatable      :: v(DTIMES(:), :)
    ! Diffusion coefficients at cell faces
    real(dp), allocatable      :: dc(DTIMES(:), :)
    ! Maximal fluxes at cell faces
    real(dp), allocatable      :: fmax(DTIMES(:), :)
    ! Cell-centered densities
    real(dp), allocatable      :: cc(DTIMES(:))
    real(dp)                   :: mu, fld_face, drt_fac, tmp
    real(dp)                   :: nsmall, N_inv
    real(dp)                   :: dt_cfl, dt_drt, dt_dif
    real(dp)                   :: vmean(NDIM)
    integer                    :: nc, n, m, IJK, tid
#if NDIM == 3
    integer                    :: l
#endif
    real(dp), parameter        :: small_value = 1e-100_dp

    ! Inside the dielectric, set the flux to zero. We later determine the
    ! boundary flux onto dielectrics
    if (ST_use_dielectric) then
       if (boxes(id)%cc(DTIMES(1), i_eps) > 1.0_dp) then
          boxes(id)%fc(DTIMES(:), :, flux_elec) = 0.0_dp
          return
       end if
    end if

    nc      = boxes(id)%n_cell
    inv_dr  = 1/boxes(id)%dr
    drt_fac = UC_eps0 / max(1e-100_dp, UC_elem_charge * dt)
    nsmall  = 1.0_dp ! A small density
    N_inv   = 1/gas_number_density ! For constant gas densities

    allocate(v(DTIMES(1:nc+1), NDIM))
    allocate(dc(DTIMES(1:nc+1), NDIM))
    allocate(cc(DTIMES(-1:nc+2)))
    allocate(fmax(DTIMES(1:nc+1), NDIM))
    if (ST_drt_limit_flux) then
       fmax = 0.0_dp
    end if

    ! Fill ghost cells on the sides of boxes (no corners)
    call af_gc_box(boxes, id, i_electron+s_in, af_gc_interp_lim, &
         af_bc_neumann_zero, .false.)

    ! Fill cc with interior data plus a second layer of ghost cells
    call af_gc2_box(boxes, id, i_electron+s_in, af_gc2_prolong_linear, &
         af_bc2_neumann_zero, cc, nc)

    ! We use the average field to compute the mobility and diffusion coefficient
    ! at the interface
    do n = 1, nc+1
       do m = 1, nc
#if NDIM == 2
          if (.not. gas_constant_density) then
             N_inv = 2 / (boxes(id)%cc(n-1, m, i_gas_dens) + &
                  boxes(id)%cc(n, m, i_gas_dens))
          end if

          fld = 0.5_dp * (boxes(id)%cc(n-1, m, i_electric_fld) + &
               boxes(id)%cc(n, m, i_electric_fld))
          Td = fld * SI_to_Townsend * N_inv
          mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
          fld_face = boxes(id)%fc(n, m, 1, electric_fld)
          v(n, m, 1)  = -mu * fld_face
          dc(n, m, 1) = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

          if (ST_drt_limit_flux) then
             tmp = abs(cc(n-1, m) - cc(n, m))/max(cc(n-1, m), cc(n, m), nsmall)
             tmp = max(fld, tmp * inv_dr(1) * dc(n, m, 1) / mu)
             fmax(n, m, 1) = drt_fac * tmp
          end if

          if (abs(fld_face) > ST_diffusion_field_limit .and. &
               (cc(n-1, m) - cc(n, m)) * fld_face > 0.0_dp) then
             dc(n, m, 1) = 0.0_dp
          end if

          if (.not. gas_constant_density) then
             N_inv = 2 / (boxes(id)%cc(m, n-1, i_gas_dens) + &
                  boxes(id)%cc(m, n, i_gas_dens))
          end if

          fld = 0.5_dp * (boxes(id)%cc(m, n-1, i_electric_fld) + &
               boxes(id)%cc(m, n, i_electric_fld))
          Td = fld * SI_to_Townsend * N_inv
          mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
          fld_face = boxes(id)%fc(m, n, 2, electric_fld)
          v(m, n, 2)  = -mu * fld_face
          dc(m, n, 2) = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

          if (ST_drt_limit_flux) then
             tmp = abs(cc(m, n-1) - cc(m, n))/max(cc(m, n-1), cc(m, n), nsmall)
             tmp = max(fld, tmp * inv_dr(2) * dc(m, n, 2) / mu)
             fmax(m, n, 2) = drt_fac * tmp
          end if

          if (abs(fld_face) > ST_diffusion_field_limit .and. &
               (cc(m, n-1) - cc(m, n)) * fld_face > 0.0_dp) then
             dc(m, n, 2) = 0.0_dp
          end if
#elif NDIM == 3
          do l = 1, nc
             if (.not. gas_constant_density) then
                N_inv = 2 / (boxes(id)%cc(n-1, m, l, i_gas_dens) + &
                     boxes(id)%cc(n, m, l, i_gas_dens))
             end if

             fld = 0.5_dp * (&
                  boxes(id)%cc(n-1, m, l, i_electric_fld) + &
                  boxes(id)%cc(n, m, l, i_electric_fld))
             Td = fld * SI_to_Townsend * N_inv
             mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
             fld_face = boxes(id)%fc(n, m, l, 1, electric_fld)
             v(n, m, l, 1)  = -mu * fld_face
             dc(n, m, l, 1) = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

             if (ST_drt_limit_flux) then
                tmp = abs(cc(n-1, m, l) - cc(n, m, l)) / &
                     max(cc(n-1, m, l), cc(n, m, l), nsmall)
                tmp = max(fld, tmp * inv_dr(1) * dc(n, m, l, 1) / mu)
                fmax(n, m, l, 1) = drt_fac * tmp
             end if

             if (abs(fld_face) > ST_diffusion_field_limit .and. &
                  (cc(n-1, m, l) - cc(n, m, l)) * fld_face > 0.0_dp) then
                dc(n, m, l, 1) = 0.0_dp
             end if

             if (.not. gas_constant_density) then
                N_inv = 2 / (boxes(id)%cc(m, n-1, l, i_gas_dens) + &
                     boxes(id)%cc(m, n, l, i_gas_dens))
             end if

             fld = 0.5_dp * (&
                  boxes(id)%cc(m, n-1, l, i_electric_fld) + &
                  boxes(id)%cc(m, n, l, i_electric_fld))
             Td = fld * SI_to_Townsend * N_inv
             mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
             fld_face = boxes(id)%fc(m, n, l, 2, electric_fld)
             v(m, n, l, 2)  = -mu * fld_face
             dc(m, n, l, 2) = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

             if (ST_drt_limit_flux) then
                tmp = abs(cc(m, n-1, l) - cc(m, n, l)) / &
                     max(cc(m, n-1, l), cc(m, n, l), nsmall)
                tmp = max(fld, tmp * inv_dr(2) * dc(m, n, l, 2) / mu)
                fmax(m, n, l, 2) = drt_fac * tmp
             end if

             if (abs(fld_face) > ST_diffusion_field_limit .and. &
                  (cc(m, n-1, l) - cc(m, n, l)) * fld_face > 0.0_dp) then
                dc(m, n, l, 2) = 0.0_dp
             end if

             if (.not. gas_constant_density) then
                N_inv = 2 / (boxes(id)%cc(m, l, n-1, i_gas_dens) + &
                     boxes(id)%cc(m, l, n, i_gas_dens))
             end if

             fld = 0.5_dp * (&
                  boxes(id)%cc(m, l, n-1, i_electric_fld) + &
                  boxes(id)%cc(m, l, n, i_electric_fld))
             Td = fld * SI_to_Townsend * N_inv
             mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
             fld_face = boxes(id)%fc(m, l, n, 3, electric_fld)
             v(m, l, n, 3)  = -mu * fld_face
             dc(m, l, n, 3) = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

             if (ST_drt_limit_flux) then
                tmp = abs(cc(m, l, n-1) - cc(m, l, n)) / &
                     max(cc(m, l, n-1), cc(m, l, n), nsmall)
                tmp = max(fld, tmp * inv_dr(3) * dc(m, l, n, 3) / mu)
                fmax(m, l, n, 3) = drt_fac * tmp
             end if

             if (abs(fld_face) > ST_diffusion_field_limit .and. &
                  (cc(m, l, n-1) - cc(m, l, n)) * fld_face > 0.0_dp) then
                dc(m, l, n, 3) = 0.0_dp
             end if
          end do
#endif
       end do
    end do

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
          dt_cfl = 1.0_dp/max(sum(abs(vmean) * inv_dr), small_value)

          ! Diffusion condition
          dt_dif = minval(boxes(id)%dr)**2 / &
               max(2 * NDIM * maxval(dc(IJK, :)), small_value)

          ! Take the combined CFL-diffusion condition with Courant number 0.5
          dt_cfl = 0.5_dp/(1/dt_cfl + 1/dt_dif)

          if (ST_drt_limit_flux) then
             mu = ST_ion_mobility
          else
             if (gas_constant_density) then
                N_inv = 1 / gas_number_density
             else
                N_inv = 1 / boxes(id)%cc(IJK, i_gas_dens)
             end if

             Td = boxes(id)%cc(IJK, i_electric_fld) * SI_to_Townsend * N_inv
             mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
             mu = max(mu, ST_ion_mobility)
          end if

          ! Dielectric relaxation time
          dt_drt = UC_eps0 / max(UC_elem_charge * mu * cc(IJK), small_value)

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

    boxes(id)%fc(DTIMES(:), :, flux_elec) = v + dc

    if (ST_drt_limit_flux) then
       where (abs(boxes(id)%fc(DTIMES(:), :, flux_elec)) > fmax)
          boxes(id)%fc(DTIMES(:), :, flux_elec) = &
               sign(fmax, boxes(id)%fc(DTIMES(:), :, flux_elec))
       end where
    end if

  end subroutine fluxes_elec

  !> Compute the ions fluxes
  subroutine fluxes_ions(tree, id, s_in)
    use m_af_flux_schemes
    use m_units_constants
    use m_streamer
    use m_gas
    use m_transport_data
    type(af_t), intent(inout)  :: tree
    integer, intent(in)        :: id
    integer, intent(in)        :: s_in   !< Input time state

    real(dp)                   :: inv_dr(NDIM)
    ! Velocities at cell faces
    real(dp), allocatable      :: v(DTIMES(:), :)
    ! Cell-centered densities
    real(dp), allocatable      :: cc(DTIMES(:))
    real(dp)                   :: mu
    integer                    :: nc, n, m, ix, i_flux, i_ion
#if NDIM == 3
    integer                    :: l
#endif

    nc      = tree%boxes(id)%n_cell
    inv_dr  = 1/tree%boxes(id)%dr

    ! Inside the dielectric, set the flux to zero
    if (ST_use_dielectric) then
       if (tree%boxes(id)%cc(DTIMES(1), i_eps) > 1.0_dp) then
          do ix = 2, size(flux_species)
             i_flux = flux_variables(ix)
             tree%boxes(id)%fc(DTIMES(:), :, i_flux) = 0.0_dp
          end do
          return
       end if
    end if

    ! Fill ghost cells on the sides of boxes (no corners)
    call af_gc_box(tree, id, flux_species(2:)+s_in, corners=.false.)

    allocate(v(DTIMES(1:nc+1), NDIM))
    allocate(cc(DTIMES(-1:nc+2)))

    do ix = 2, size(flux_species)
       i_ion  = flux_species(ix)
       i_flux = flux_variables(ix)
       mu     = transport_data_ions%mobilities(ix-1)

       ! Fill cc with interior data plus a second layer of ghost cells
       call af_gc2_box(tree%boxes, id, i_ion+s_in, af_gc2_prolong_linear, &
            af_bc2_neumann_zero, cc, nc)

       do n = 1, nc+1
          do m = 1, nc
#if NDIM == 2
             v(n, m, 1) = mu * tree%boxes(id)%fc(n, m, 1, electric_fld)
             v(m, n, 2) = mu * tree%boxes(id)%fc(m, n, 2, electric_fld)
#elif NDIM == 3
             do l = 1, nc
                v(n, m, l, 1) = mu * tree%boxes(id)%fc(n, m, l, 1, electric_fld)
                v(m, n, l, 2) = mu * tree%boxes(id)%fc(m, n, l, 2, electric_fld)
                v(m, l, n, 3) = mu * tree%boxes(id)%fc(m, l, n, 3, electric_fld)
             end do
#endif
          end do
       end do

#if NDIM == 2
       call flux_koren_2d(cc, v, nc, 2)
#elif NDIM == 3
       call flux_koren_3d(cc, v, nc, 2)
#endif

       tree%boxes(id)%fc(DTIMES(:), :, i_flux) = v
    end do

  end subroutine fluxes_ions

  !> Advance solution in a box over dt based on the fluxes and reactions, using
  !> a forward Euler update
  subroutine update_solution(box, id, dt, s_in, s_out, set_dt)
    use omp_lib
    use m_units_constants
    use m_gas
    use m_chemistry
    use m_streamer
    use m_photoi
    use m_dt
    use m_dielectric
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: id     !< id of the box
    real(dp), intent(in)       :: dt     !< Time step
    integer, intent(in)        :: s_in   !< Input time state
    integer, intent(in)        :: s_out  !< Output time state
    logical, intent(in)        :: set_dt !< Whether to set new time step
    real(dp)                   :: inv_dr(NDIM)
    real(dp)                   :: tmp
    real(dp), allocatable      :: rates(:, :)
    real(dp), allocatable      :: derivs(:, :)
    real(dp), allocatable      :: dens(:, :)
    real(dp), allocatable      :: fields(:)
#if NDIM == 2
    real(dp)                   :: rfac(2, box%n_cell)
#endif
    integer                    :: IJK, ix, nc, n_cells, n, iv
    integer                    :: tid, i_flux
    real(dp), parameter        :: eps = 1e-100_dp

    nc      = box%n_cell
    n_cells = box%n_cell**NDIM
    inv_dr  = 1/box%dr

    ! Inside the dielectric, do not update the species densities, which should
    ! always be zero
    if (ST_use_dielectric) then
       if (box%cc(DTIMES(1), i_eps) > 1.0_dp) then
          return
       end if
    end if

    allocate(rates(n_cells, n_reactions))
    allocate(derivs(n_cells, n_species))
    allocate(fields(n_cells))
    allocate(dens(n_cells, n_species))

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
         species_itree(n_gas_species+1:n_species)+s_in), [n_cells, n_plasma_species])

    call get_rates(fields, rates, n_cells)
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
       ! Contribution of electron flux
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

       do n = n_gas_species+1, n_species
          iv = species_itree(n)
          box%cc(IJK, iv+s_out) = box%cc(IJK, iv+s_in) + dt * derivs(ix, n)
       end do
    end do; CLOSE_DO

    ! Ion fluxes
    do n = 2, size(flux_species)
       iv     = flux_species(n)
       i_flux = flux_variables(n)

       do KJI_DO(1,nc)
          ! Contribution of ion flux
#if NDIM == 2
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

    if (ST_use_dielectric) then
       ! Check if the box is part of a surface
       ix = box_id_to_surface_id(id)
       if (ix > 0) then
          if (id == surface_list(ix)%id_gas) then
             ! Convert fluxes onto dielectric to surface charge, and handle
             ! secondary emission
             call dielectric_update_surface_charge(box, &
                  surface_list(ix), dt, s_in, s_out)

             ! Add secondary emission from photons hitting the surface
             call dielectric_photon_emission(box, &
                  surface_list(ix), dt, s_in, s_out)
          end if
       end if
    end if

  end subroutine update_solution

end module m_fluid_lfa
