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

  subroutine forward_euler(tree, dt, s_in, s_out, set_dt)
    use m_chemistry
    use m_streamer
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt
    integer, intent(in)       :: s_in
    integer, intent(in)       :: s_out
    logical, intent(in)       :: set_dt
    integer                   :: lvl, i, id, p_id

    ! First calculate fluxes
    !$omp parallel private(lvl, i, id)
    do lvl = tree%lowest_lvl, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call fluxes_elec(tree%boxes, id, dt, s_in, set_dt)
       end do
       !$omp end do
    end do
    !$omp end parallel

    call af_consistent_fluxes(tree, [flux_elec])

    ! Update the solution
    !$omp parallel private(lvl, i, id, p_id)
    do lvl = tree%lowest_lvl, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call update_solution(tree%boxes(id), dt, s_in, s_out, set_dt)

          ! This can be important for setting ghost cells
          p_id = tree%boxes(id)%parent
          if (p_id > af_no_box) then
             call af_restrict_box_vars(tree%boxes(id), tree%boxes(p_id), &
                  flux_species)
          end if
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
    integer, intent(in)        :: s_in
    logical, intent(in)        :: set_dt
    real(dp)                   :: inv_dr, fld, Td
    ! Velocities at cell faces
    real(dp), allocatable      :: v(DTIMES(:), :)
    ! Diffusion coefficients at cell faces
    real(dp), allocatable      :: dc(DTIMES(:), :)
    ! Maximal fluxes at cell faces
    real(dp), allocatable      :: fmax(DTIMES(:), :)
    ! Cell-centered densities
    real(dp), allocatable      :: cc(DTIMES(:))
    real(dp)                   :: mu, fld_face, drt_fac, tmp
    real(dp)                   :: nsmall, N_inv, dr
    real(dp)                   :: dt_cfl, dt_drt, dt_dif
    real(dp)                   :: vmean(NDIM)
    integer                    :: nc, n, m, IJK, tid
#if NDIM == 3
    integer                    :: l
#endif

    nc      = boxes(id)%n_cell
    dr      = boxes(id)%dr
    inv_dr  = 1/boxes(id)%dr
    drt_fac = UC_eps0 / max(1e-100_dp, UC_elem_charge * dt)
    nsmall  = 1.0_dp ! A small density
    N_inv   = 1/gas_number_density

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
          fld = 0.5_dp * (boxes(id)%cc(n-1, m, i_electric_fld) + &
               boxes(id)%cc(n, m, i_electric_fld))
          Td = fld * SI_to_Townsend
          mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
          fld_face = boxes(id)%fc(n, m, 1, electric_fld)
          v(n, m, 1)  = -mu * fld_face
          dc(n, m, 1) = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

          if (ST_drt_limit_flux) then
             tmp = abs(cc(n-1, m) - cc(n, m))/max(cc(n-1, m), cc(n, m), nsmall)
             tmp = max(fld, tmp * inv_dr * dc(n, m, 1) / mu)
             fmax(n, m, 1) = drt_fac * tmp
          end if

          if (abs(fld_face) > ST_diffusion_field_limit .and. &
               (cc(n-1, m) - cc(n, m)) * fld_face > 0.0_dp) then
             dc(n, m, 1) = 0.0_dp
          end if

          fld = 0.5_dp * (boxes(id)%cc(m, n-1, i_electric_fld) + &
               boxes(id)%cc(m, n, i_electric_fld))
          Td = fld * SI_to_Townsend
          mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
          fld_face = boxes(id)%fc(m, n, 2, electric_fld)
          v(m, n, 2)  = -mu * fld_face
          dc(m, n, 2) = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

          if (ST_drt_limit_flux) then
             tmp = abs(cc(m, n-1) - cc(m, n))/max(cc(m, n-1), cc(m, n), nsmall)
             tmp = max(fld, tmp * inv_dr * dc(m, n, 2) / mu)
             fmax(m, n, 2) = drt_fac * tmp
          end if

          if (abs(fld_face) > ST_diffusion_field_limit .and. &
               (cc(m, n-1) - cc(m, n)) * fld_face > 0.0_dp) then
             dc(m, n, 2) = 0.0_dp
          end if
#elif NDIM == 3
          do l = 1, nc
             fld = 0.5_dp * (&
                  boxes(id)%cc(n-1, m, l, i_electric_fld) + &
                  boxes(id)%cc(n, m, l, i_electric_fld))
             Td = fld * SI_to_Townsend
             mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
             fld_face = boxes(id)%fc(n, m, l, 1, electric_fld)
             v(n, m, l, 1)  = -mu * fld_face
             dc(n, m, l, 1) = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

             if (ST_drt_limit_flux) then
                tmp = abs(cc(n-1, m, l) - cc(n, m, l)) / &
                     max(cc(n-1, m, l), cc(n, m, l), nsmall)
                tmp = max(fld, tmp * inv_dr * dc(n, m, l, 1) / mu)
                fmax(n, m, l, 1) = drt_fac * tmp
             end if

             if (abs(fld_face) > ST_diffusion_field_limit .and. &
                  (cc(n-1, m, l) - cc(n, m, l)) * fld_face > 0.0_dp) then
                dc(n, m, l, 1) = 0.0_dp
             end if

             fld = 0.5_dp * (&
                  boxes(id)%cc(m, n-1, l, i_electric_fld) + &
                  boxes(id)%cc(m, n, l, i_electric_fld))
             Td = fld * SI_to_Townsend
             mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
             fld_face = boxes(id)%fc(m, n, l, 2, electric_fld)
             v(m, n, l, 2)  = -mu * fld_face
             dc(m, n, l, 2) = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

             if (ST_drt_limit_flux) then
                tmp = abs(cc(m, n-1, l) - cc(m, n, l)) / &
                     max(cc(m, n-1, l), cc(m, n, l), nsmall)
                tmp = max(fld, tmp * inv_dr * dc(m, n, l, 2) / mu)
                fmax(m, n, l, 2) = drt_fac * tmp
             end if

             if (abs(fld_face) > ST_diffusion_field_limit .and. &
                  (cc(m, n-1, l) - cc(m, n, l)) * fld_face > 0.0_dp) then
                dc(m, n, l, 2) = 0.0_dp
             end if

             fld = 0.5_dp * (&
                  boxes(id)%cc(m, l, n-1, i_electric_fld) + &
                  boxes(id)%cc(m, l, n, i_electric_fld))
             Td = fld * SI_to_Townsend
             mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
             fld_face = boxes(id)%fc(m, l, n, 3, electric_fld)
             v(m, l, n, 3)  = -mu * fld_face
             dc(m, l, n, 3) = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

             if (ST_drt_limit_flux) then
                tmp = abs(cc(m, l, n-1) - cc(m, l, n)) / &
                     max(cc(m, l, n-1), cc(m, l, n), nsmall)
                tmp = max(fld, tmp * inv_dr * dc(m, l, n, 3) / mu)
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
          dt_cfl = 1.0_dp/sum(abs(vmean) * inv_dr)

          ! Diffusion condition
          dt_dif = dr**2 / max(2 * NDIM * maxval(dc(IJK, :)), epsilon(1.0_dp))

          ! Take the combined CFL-diffusion condition with Courant number 0.5
          dt_cfl = 0.5_dp/(1/dt_cfl + 1/dt_dif)

          if (ST_drt_limit_flux) then
             mu = ST_ion_mobility
          else
             Td = boxes(id)%cc(IJK, i_electric_fld) * SI_to_Townsend
             mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
             mu = max(mu, ST_ion_mobility)
          end if

          ! Dielectric relaxation time
          dt_drt = UC_eps0 / max(UC_elem_charge * mu * cc(IJK), epsilon(1.0_dp))

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

    !     if (ST_update_ions) then
    !        ! Use a constant diffusion coefficient for ions
    !        dc = ST_ion_diffusion

    !        ! Use a constant mobility for ions
    !        v(DTIMES(:), 1:NDIM) = ST_ion_mobility * &
    !             boxes(id)%fc(DTIMES(:), 1:NDIM, electric_fld)

    !        ! Fill ghost cells on the sides of boxes (no corners)
    !        call af_gc_box(boxes, id, i_pos_ion, af_gc_interp_lim, &
    !             af_bc_neumann_zero, .false.)

    !        ! Fill cc with interior data plus a second layer of ghost cells
    !        call af_gc2_box(boxes, id, i_pos_ion, af_gc2_prolong_linear, &
    !             af_bc2_neumann_zero, cc, nc)

    ! #if NDIM == 2
    !        call flux_koren_2d(cc, v, nc, 2)
    !        call flux_diff_2d(cc, dc, inv_dr, nc, 2)
    ! #elif NDIM == 3
    !        call flux_koren_3d(cc, v, nc, 2)
    !        call flux_diff_3d(cc, dc, inv_dr, nc, 2)
    ! #endif

    !        boxes(id)%fc(DTIMES(:), :, flux_ion) = v + dc
    !     end if
  end subroutine fluxes_elec

  !> Advance solution over dt based on the fluxes / source term, using forward Euler
  subroutine update_solution(box, dt, s_in, s_out, set_dt)
    use omp_lib
    use m_units_constants
    use m_gas
    use m_chemistry
    use m_streamer
    use m_photoi
    use m_dt
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: dt
    integer, intent(in)        :: s_in
    integer, intent(in)        :: s_out
    logical, intent(in)        :: set_dt
    real(dp)                   :: inv_dr
    real(dp), allocatable      :: rates(:, :)
    real(dp), allocatable      :: derivs(:, :)
    real(dp), allocatable      :: dens(:, :)
    real(dp), allocatable      :: fields(:)
#if NDIM == 2
    real(dp)                   :: rfac(2)
    integer                    :: ioff
#endif
    integer                    :: IJK, ix, nc, n_cells, n, iv
    integer                    :: tid

    nc     = box%n_cell
    n_cells = box%n_cell**NDIM
    inv_dr = 1/box%dr
#if NDIM == 2
    ioff   = (box%ix(1)-1) * nc
#endif

    allocate(rates(n_cells, n_reactions))
    allocate(derivs(n_cells, n_species))

    fields = SI_to_Townsend * &
         reshape(box%cc(DTIMES(1:nc), i_electric_fld), [n_cells])

    dens = reshape(box%cc(DTIMES(1:nc), species_ix(1:n_species)+s_in), &
         [n_cells, n_species])

    call get_rates(fields, rates, n_cells)
    call get_derivatives(dens, rates, derivs, n_cells)

    if (set_dt) then
       tid = omp_get_thread_num() + 1
       dt_matrix(dt_ix_rates, tid) = min(dt_matrix(dt_ix_rates, tid), &
            1.0_dp / max(maxval(rates), epsilon(1.0_dp)))
    end if

    ix = 0
    do KJI_DO(1,nc)
       ix = ix + 1
       ! Contribution of flux
#if NDIM == 2
       if (ST_cylindrical) then
          ! Weighting of flux contribution for cylindrical coordinates
          rfac(:) = [i+ioff-1, i+ioff] / (i+ioff-0.5_dp)
       else
          rfac(:) = 1.0_dp
       end if

       derivs(ix, i_electron) = derivs(ix, i_electron) + inv_dr * ( &
            box%fc(i, j, 2, flux_elec) - &
            box%fc(i, j+1, 2, flux_elec) + &
            rfac(1) * box%fc(i, j, 1, flux_elec) - &
            rfac(2) * box%fc(i+1, j, 1, flux_elec))
#elif NDIM == 3
       derivs(ix, i_electron) = derivs(ix, i_electron) + inv_dr * ( &
            sum(box%fc(i, j, k, 1:3, flux_elec)) - &
            box%fc(i+1, j, k, 1, flux_elec) - &
            box%fc(i, j+1, k, 2, flux_elec) - &
            box%fc(i, j, k+1, 3, flux_elec))
#endif

       if (photoi_enabled) then
          derivs(ix, i_electron) = derivs(ix, i_electron) + &
               box%cc(IJK, i_photo)
          derivs(ix, i_1pos_ion) = derivs(ix, i_1pos_ion) + &
               box%cc(IJK, i_photo)
       end if

       do n = 1, n_species
          iv = species_ix(n)
          box%cc(IJK, iv+s_out) = box%cc(IJK, iv+s_in) + dt * derivs(ix, n)
       end do
    end do; CLOSE_DO

  end subroutine update_solution

end module m_fluid_lfa
