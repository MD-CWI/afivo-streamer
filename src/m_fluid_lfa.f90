#include "../afivo/src/cpp_macros.h"
!> Fluid model module
module m_fluid_lfa
  use m_af_all
  use m_streamer

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
  subroutine forward_euler(tree, dt, dt_lim, time, s_deriv, n_prev, s_prev, &
       w_prev, s_out, i_step, n_steps)
    use m_chemistry
    use m_field
    use m_dt
    use m_transport_data
    use m_dielectric
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt             !< Time step
    real(dp), intent(inout)   :: dt_lim         !< Computed time step limit
    real(dp), intent(in)      :: time           !< Current time
    integer, intent(in)       :: s_deriv        !< State to compute derivatives from
    integer, intent(in)       :: n_prev         !< Number of previous states
    integer, intent(in)       :: s_prev(n_prev) !< Previous states
    real(dp), intent(in)      :: w_prev(n_prev) !< Weights of previous states
    integer, intent(in)       :: s_out          !< Output state
    integer, intent(in)       :: i_step         !< Step of the integrator
    integer, intent(in)       :: n_steps        !< Total number of steps
    integer                   :: lvl, i, id, p_id, nc, ix, id_out
    logical                   :: last_step

    nc = tree%n_cell

    ! Set current rates to zero; they are summed below
    ST_current_rates = 0
    ST_current_JdotE = 0

    last_step = (i_step == n_steps)

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
          call fluxes_elec(tree, id, nc, dt, s_deriv, last_step)

          if (transport_data_ions%n_mobile_ions > 0) then
             call fluxes_ions(tree, id, nc, dt, s_deriv, last_step)
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
               n_prev, s_prev, w_prev, s_out, last_step)

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

    if (last_step) then
       dt_lim = minval(dt_matrix(1:dt_num_cond, :))
    end if
  end subroutine forward_euler

  !> Get velocity and diffusion coefficient for electron flux
  subroutine compute_flux_coeff_1d(nc, E_cc, E_x, ne, N_gas, dt, dx, v, dc, fmax)
    use m_gas
    use m_units_constants
    use m_lookup_table
    use m_transport_data
    integer, intent(in)     :: nc
    real(dp), intent(in)    :: E_cc(0:nc+1)  !< Cell-centered field strengths
    real(dp), intent(in)    :: E_x(nc+1)     !< Face-centered field components
    real(dp), intent(in)    :: ne(0:nc+1)    !< Electron densities
    real(dp), intent(in)    :: N_gas(0:nc+1) !< Gas number density at cell centers
    real(dp), intent(in)    :: dt            !< Current time step
    real(dp), intent(in)    :: dx            !< Grid spacing
    real(dp), intent(out)   :: v(nc+1)       !< Velocity at cell faces
    real(dp), intent(out)   :: dc(nc+1)      !< Diffusion coefficient at cell faces
    real(dp), intent(inout) :: fmax(nc+1)    !< Maximum allowed flux

    real(dp), parameter :: nsmall = 1.0_dp ! A small density
    integer             :: n
    real(dp)            :: E_face(nc+1), Td(nc+1), N_inv(nc+1)
    real(dp)            :: mu(nc+1), tmp, drt_fac, inv_dx

    if (gas_constant_density) then
       N_inv = 1/gas_number_density
    else
       ! Compute gas number density at cell faces
       do n = 1, nc+1
          N_inv(n) = 2 / (N_gas(n-1) + N_gas(n))
       end do
    end if

    ! Compute field strength at cell faces, which is used to compute the
    ! mobility and diffusion coefficient at the interface
    do n = 1, nc+1
       E_face(n) = 0.5_dp * (E_cc(n-1) + E_cc(n))
       Td(n) = E_face(n) * SI_to_Townsend * N_inv(n)
    end do

    mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
    dc = LT_get_col(td_tbl, td_diffusion, Td) * N_inv

    ! Compute velocity, -mu accounts for negative electron charge
    v = -mu * E_x

    if (ST_drt_limit_flux) then
       !> Compute maximal flux if the dielectric relaxation time is not taken
       !> into account, see https://doi.org/10.1088/1361-6595/ab6757
       drt_fac = UC_eps0 / max(1e-100_dp, UC_elem_charge * dt)
       inv_dx = 1/dx

       do n = 1, nc+1
          tmp = abs(ne(n-1) - ne(n)) / max(ne(n-1), ne(n), nsmall)
          tmp = max(E_face(n), tmp * inv_dx * dc(n) / mu(n))
          fmax(n) = drt_fac * tmp
       end do
    end if
  end subroutine compute_flux_coeff_1d

  !> Compute the electron fluxes due to drift and diffusion
  subroutine fluxes_elec(tree, id, nc, dt, s_in, last_step)
    use m_af_flux_schemes
    use m_units_constants
    use omp_lib
    use m_gas
    use m_dt
    use m_lookup_table
    use m_transport_data
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: id
    integer, intent(in)       :: nc   !< Number of cells per dimension
    real(dp), intent(in)      :: dt
    integer, intent(in)       :: s_in !< Input time state
    logical, intent(in)       :: last_step
    real(dp)                  :: dr(NDIM), inv_dr(NDIM), Td
    ! Velocities at cell faces
    real(dp)                  :: v(DTIMES(1:nc+1), NDIM)
    ! Diffusion coefficients at cell faces
    real(dp)                  :: dc(DTIMES(1:nc+1), NDIM)
    ! Maximal fluxes at cell faces
    real(dp)                  :: fmax(DTIMES(1:nc+1), NDIM)
    ! Cell-centered densities
    real(dp)                  :: cc(DTIMES(-1:nc+2), 1)

    real(dp)                  :: mu, max_mu_ion, N_inv
    real(dp)                  :: dt_cfl, dt_drt, dt_dif
    real(dp)                  :: vmean(NDIM), N_gas(0:nc+1)
    real(dp)                  :: v_x(nc+1), dc_x(nc+1), fmax_x(nc+1)
    real(dp)                  :: E_cc(0:nc+1), E_fc(nc+1), ne(0:nc+1)
    real(dp), parameter       :: eps = 1e-100_dp
    integer                   :: IJK, tid, dir
#if NDIM == 2
    integer                   :: m
#elif NDIM == 3
    integer                   :: m, n
#endif

    ! Inside the dielectric, set the flux to zero. We later determine the
    ! boundary flux onto dielectrics
    if (ST_use_dielectric) then
       if (tree%boxes(id)%cc(DTIMES(1), i_eps) > 1.0_dp) then
          tree%boxes(id)%fc(DTIMES(:), :, flux_elec) = 0.0_dp
          return
       end if
    end if

    dr     = tree%boxes(id)%dr
    inv_dr = 1/tree%boxes(id)%dr
    v      = 0.0_dp
    dc     = 0.0_dp
    fmax   = 0.0_dp

    ! Fill cc with interior data plus two layers of ghost cells
    call af_gc2_box(tree, id, [i_electron+s_in], cc)

    associate(box => tree%boxes(id))
#if NDIM == 1
      dir = 1                ! x-component
      if (.not. gas_constant_density) N_gas = box%cc(0:nc+1, i_gas_dens)

      E_cc = box%cc(0:nc+1, i_electric_fld)
      E_fc = box%fc(:, dir, electric_fld)
      ne = cc(0:nc+1, 1)
      call compute_flux_coeff_1d(nc, E_cc, E_fc, ne, N_gas, dt, dr(dir), &
           v_x, dc_x, fmax_x)
      v(:, dir) = v_x
      dc(:, dir) = dc_x
      fmax(:, dir) = fmax_x
#elif NDIM == 2
      do m = 1, nc
         dir = 1                ! x-component
         if (.not. gas_constant_density) N_gas = box%cc(0:nc+1, m, i_gas_dens)

         ! Avoid allocating array temporaries, but explicitly copy
         E_cc = box%cc(0:nc+1, m, i_electric_fld)
         E_fc = box%fc(:, m, dir, electric_fld)
         ne = cc(0:nc+1, m, 1)
         call compute_flux_coeff_1d(nc, E_cc, E_fc, ne, N_gas, dt, dr(dir), &
              v_x, dc_x, fmax_x)
         v(:, m, dir) = v_x
         dc(:, m, dir) = dc_x
         fmax(:, m, dir) = fmax_x

         dir = 2                ! y-component
         if (.not. gas_constant_density) N_gas = box%cc(m, 0:nc+1, i_gas_dens)

         E_cc = box%cc(m, 0:nc+1, i_electric_fld)
         E_fc = box%fc(m, :, dir, electric_fld)
         ne = cc(m, 0:nc+1, 1)
         call compute_flux_coeff_1d(nc, E_cc, E_fc, ne, N_gas, dt, dr(dir), &
              v_x, dc_x, fmax_x)
         v(m, :, dir) = v_x
         dc(m, :, dir) = dc_x
         fmax(m, :, dir) = fmax_x
      end do
#elif NDIM == 3
      do n = 1, nc
         do m = 1, nc
            dir = 1                ! x-component
            if (.not. gas_constant_density) N_gas = box%cc(0:nc+1, m, n, i_gas_dens)

            E_cc = box%cc(0:nc+1, m, n, i_electric_fld)
            E_fc = box%fc(:, m, n, dir, electric_fld)
            ne = cc(0:nc+1, m, n, 1)
            call compute_flux_coeff_1d(nc, E_cc, E_fc, ne, N_gas, dt, dr(dir), &
                 v_x, dc_x, fmax_x)
            v(:, m, n, dir) = v_x
            dc(:, m, n, dir) = dc_x
            fmax(:, m, n, dir) = fmax_x

            dir = 2                ! y-component
            if (.not. gas_constant_density) N_gas = box%cc(m, 0:nc+1, n, i_gas_dens)

            E_cc = box%cc(m, 0:nc+1, n, i_electric_fld)
            E_fc = box%fc(m, :, n, dir, electric_fld)
            ne = cc(m, 0:nc+1, n, 1)
            call compute_flux_coeff_1d(nc, E_cc, E_fc, ne, N_gas, dt, dr(dir), &
                 v_x, dc_x, fmax_x)
            v(m, :, n, dir) = v_x
            dc(m, :, n, dir) = dc_x
            fmax(m, :, n, dir) = fmax_x

            dir = 3                ! z-component
            if (.not. gas_constant_density) N_gas = box%cc(m, n, 0:nc+1, i_gas_dens)

            E_cc = box%cc(m, n, 0:nc+1, i_electric_fld)
            E_fc = box%fc(m, n, :, dir, electric_fld)
            ne = cc(m, n, 0:nc+1, 1)
            call compute_flux_coeff_1d(nc, E_cc, E_fc, ne, N_gas, dt, dr(dir), &
                 v_x, dc_x, fmax_x)
            v(m, n, :, dir) = v_x
            dc(m, n, :, dir) = dc_x
            fmax(m, n, :, dir) = fmax_x
         end do
      end do
#endif
    end associate

    if (last_step) then
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

          ! Take combined CFL-diffusion condition
          dt_cfl = dt_cfl_number/(1/dt_cfl + 1/dt_dif)

          if (gas_constant_density) then
             N_inv = 1 / gas_number_density
          else
             N_inv = 1 / tree%boxes(id)%cc(IJK, i_gas_dens)
          end if

          ! Ion mobility
          if (size(transport_data_ions%mobilities) > 0) then
             ! This should be an over-estimate
             max_mu_ion = maxval(abs(transport_data_ions%mobilities)) * N_inv * &
                  ST_sheath_max_ion_mobility_factor
          else
             max_mu_ion = 0.0_dp
          end if

          ! Electron mobility
          Td = tree%boxes(id)%cc(IJK, i_electric_fld) * SI_to_Townsend * N_inv
          mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv

          ! Take sum of electron and ion mobility, in the future we should
          ! probably apply a weighted sum (weighted with species densities)
          mu = mu + max_mu_ion

          ! Dielectric relaxation time
          dt_drt = UC_eps0 / max(UC_elem_charge * mu * cc(IJK, 1), eps)

          if (ST_drt_limit_flux) then
             ! By limiting the flux, we reduce the conductivity of the cell. If
             ! the current through the cell is fixed, this ensures that the
             ! field in the cell will not exceed ST_drt_max_field
             dt_drt = dt_drt * max(1.0_dp, ST_drt_max_field / &
                  max(1e-10_dp, tree%boxes(id)%cc(IJK, i_electric_fld)))
          end if

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

    if (ST_source_factor == source_factor_original_flux) then
       ! Store approximation of E . [-D grad(n_e)] in temporary variable
       associate (box => tree%boxes(id))
         do KJI_DO(1, nc)
#if NDIM == 1
            box%cc(IJK, i_srcfac) = 0.5_dp * (&
                 box%fc(i, 1, electric_fld) * dc(i, 1) + &
                 box%fc(i+1, 1, electric_fld) * dc(i+1, 1))
#elif NDIM == 2
            box%cc(IJK, i_srcfac) = 0.5_dp * (&
                 box%fc(i, j, 1, electric_fld) * dc(i, j, 1) + &
                 box%fc(i+1, j, 1, electric_fld) * dc(i+1, j, 1) + &
                 box%fc(i, j, 2, electric_fld) * dc(i, j, 2) + &
                 box%fc(i, j+1, 2, electric_fld) * dc(i, j+1, 2))
#elif NDIM == 3
            box%cc(IJK, i_srcfac) = 0.5_dp * (&
                 box%fc(i, j, k, 1, electric_fld) * dc(i, j, k, 1) + &
                 box%fc(i+1, j, k, 1, electric_fld) * dc(i+1, j, k, 1) + &
                 box%fc(i, j, k, 2, electric_fld) * dc(i, j, k, 2) + &
                 box%fc(i, j+1, k, 2, electric_fld) * dc(i, j+1, k, 2) + &
                 box%fc(i, j, k, 3, electric_fld) * dc(i, j, k, 3) + &
                 box%fc(i, j, k+1, 3, electric_fld) * dc(i, j, k+1, 3))
#endif
         end do; CLOSE_DO
       end associate
    end if

    if (ST_drt_limit_flux) then
       where (abs(tree%boxes(id)%fc(DTIMES(:), :, flux_elec)) > fmax)
          tree%boxes(id)%fc(DTIMES(:), :, flux_elec) = &
               sign(fmax, tree%boxes(id)%fc(DTIMES(:), :, flux_elec))
       end where
    end if

  end subroutine fluxes_elec

  subroutine fluxes_ions(tree, id, nc, dt, s_in, last_step)
    use m_af_flux_schemes
    use m_gas
    use m_transport_data
    type(af_t), intent(inout) :: tree
    integer, intent(in)        :: id
    integer, intent(in)        :: nc   !< Number of cells per dimension
    real(dp), intent(in)       :: dt
    integer, intent(in)        :: s_in !< Input time state
    logical, intent(in)        :: last_step
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
       ! Account for ion charge in mobility
       mu = transport_data_ions%mobilities(ix) * flux_species_charge(ix+1)

       associate (box => tree%boxes(id), fc => tree%boxes(id)%fc, &
            cc => tree%boxes(id)%cc)
#if NDIM == 1
         do n = 1, nc+1
            v(n, 1) = get_ion_velocity(box, n, 1, mu)
         end do
#elif NDIM == 2
         do n = 1, nc+1
            do m = 1, nc
               v(n, m, 1) = get_ion_velocity(box, n, m, 1, mu)
               v(m, n, 2) = get_ion_velocity(box, m, n, 2, mu)
            end do
         end do
#elif NDIM == 3
         do n = 1, nc+1
            do m = 1, nc
               do l = 1, nc
                  v(n, m, l, 1) = get_ion_velocity(box, n, m, l, 1, mu)
                  v(m, n, l, 2) = get_ion_velocity(box, m, n, l, 2, mu)
                  v(m, l, n, 3) = get_ion_velocity(box, m, l, n, 3, mu)
               end do
            end do
         end do
#endif
       end associate

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

  real(dp) function get_ion_velocity(box, IJK, dim, mu) result(v)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: IJK !< Flux index
    integer, intent(in)     :: dim !< Flux dimension
    real(dp), intent(in)    :: mu  !< Original mobility
    real(dp)                :: field_norm, field_face, factor, lsf_face

    field_face = box%fc(IJK, dim, electric_fld)
    factor = 1.0_dp

    if (ST_sheath_max_ion_mobility_factor > 1 .and. ST_use_electrode) then
       ! Get lsf value at the cell face, which should be (approximately) the
       ! distance to the electrode
       lsf_face = cc_average_at_cell_face(box, IJK, dim, i_lsf)

       if (lsf_face < ST_sheath_max_lsf) then
          field_norm = cc_average_at_cell_face(box, IJK, dim, i_electric_fld)
          factor = max(1.0_dp, exp(field_norm/ST_sheath_field_threshold - 1))
          factor = min(factor, ST_sheath_max_ion_mobility_factor)
       end if
    end if

    v = mu * factor * field_face * get_N_inv_face(box, IJK, dim)
  end function get_ion_velocity

  !> Get average of cell-centered quantity at a cell face
  pure function cc_average_at_cell_face(box, IJK, idim, iv) result(avg)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: IJK  !< Face index
    integer, intent(in)     :: idim !< Direction of the cell face
    integer, intent(in)     :: iv   !< Index of cell-centered variable
    real(dp)                :: avg

#if NDIM == 1
    avg = 0.5_dp * (box%cc(i-1, iv) + box%cc(i, iv))
#elif NDIM == 2
    select case (idim)
    case (1)
       avg = 0.5_dp * (box%cc(i-1, j, iv) + box%cc(i, j, iv))
    case default
       avg = 0.5_dp * (box%cc(i, j-1, iv) + box%cc(i, j, iv))
    end select
#elif NDIM == 3
    select case (idim)
    case (1)
       avg = 0.5_dp * (box%cc(i-1, j, k, iv) + box%cc(i, j, k, iv))
    case (2)
       avg = 0.5_dp * (box%cc(i, j-1, k, iv) + box%cc(i, j, k, iv))
    case default
       avg = 0.5_dp * (box%cc(i, j, k-1, iv) + box%cc(i, j, k, iv))
    end select
#endif
  end function cc_average_at_cell_face

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
       get_N_inv_face = 1 / cc_average_at_cell_face(box, IJK, idim, i_gas_dens)
    end if
  end function get_N_inv_face

  !> Advance solution in a box over dt based on the fluxes and reactions, using
  !> a forward Euler update
  subroutine update_solution(box, nc, dt, s_dt, n_prev, s_prev, w_prev, &
       s_out, last_step)
    use omp_lib
    use m_units_constants
    use m_gas
    use m_chemistry
    use m_photoi
    use m_dt
    use m_lookup_table
    use m_transport_data
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: nc             !< Box size
    real(dp), intent(in)       :: dt             !< Time step
    integer, intent(in)        :: s_dt           !< Time state to compute derivatives from
    integer, intent(in)        :: n_prev         !< Number of previous states
    integer, intent(in)        :: s_prev(n_prev) !< Time state to add derivatives to
    real(dp), intent(in)       :: w_prev(n_prev) !< Weights of previous states
    integer, intent(in)        :: s_out          !< Output time state
    logical, intent(in)        :: last_step         !< Whether to set new time step
    real(dp)                   :: inv_dr(NDIM)
    real(dp)                   :: tmp
    real(dp)                   :: rates(nc**NDIM, n_reactions)
    real(dp)                   :: derivs(nc**NDIM, n_species)
    real(dp)                   :: dens(nc**NDIM, n_species)
    real(dp)                   :: fields(nc**NDIM), box_rates(n_reactions)
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

    if (ST_source_factor /= source_factor_none .or. &
         ST_source_min_density > 0) then
       if (ST_source_factor /= source_factor_none) then
          call compute_source_factor(box, nc, dens(:, ix_electron), &
               fields, s_dt, source_factor)
       else
          source_factor(:) = 1.0_dp
       end if

       if (ST_source_min_density > 0) then
          where (dens(:, ix_electron) < ST_source_min_density)
             source_factor = 0.0_dp
          end where
       end if

       if (i_srcfac > 0) then
          ! Write source factor to variable
          box%cc(DTIMES(1:nc), i_srcfac) = &
               reshape(source_factor, [DTIMES(nc)])
       end if

       do n = 1, n_reactions
          if (reactions(n)%reaction_type == ionization_reaction) then
             rates(:, n) = rates(:, n) * source_factor
          end if
       end do
    end if

    call get_derivatives(dens, rates, derivs, n_cells)

    ! TODO clean this up in the future (when we set electrode Neumann boundary
    ! conditions in actual flux computation)
    if (iand(box%tag, mg_lsf_box) > 0) then
       ix = 0
       do KJI_DO(1,nc)
          ix = ix + 1
          if (box%cc(IJK, i_lsf) <= 0.0_dp) then
             derivs(ix, :) = 0
          end if
       end do; CLOSE_DO
    end if

    if (last_step) then
       tid = omp_get_thread_num() + 1
       dt_matrix(dt_ix_rates, tid) = min(dt_matrix(dt_ix_rates, tid), &
            minval((abs(dens) + dt_chemistry_nmin) / max(abs(derivs), eps)))

       ! Keep track of chemical production at last time integration step
       call chemical_rates_box(box, nc, rates, box_rates)

       !> Integrate rates over space and time into global storage
       ST_current_rates(1:n_reactions, tid) = &
         ST_current_rates(1:n_reactions, tid) + box_rates

       ! Keep track of J.E
       call sum_global_JdotE(box, tid)
    end if

#if NDIM == 2
    if (ST_cylindrical) then
       call af_cyl_flux_factors(box, rfac)
    else
       rfac = 1.0_dp
    end if
#endif

    if (ST_plasma_region_enabled) then
       ! Compute box coordinates
       do n = 1, NDIM
          coords(:, n) = box%r_min(n) + box%dr(n) * [(i-0.5_dp, i=1,nc)]
       end do
    end if

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
             ! Disable all reactions outside the plasma region
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
          box%cc(IJK, iv+s_out) = sum(w_prev * box%cc(IJK, iv+s_prev)) + &
               dt * derivs(ix, n)
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
  subroutine compute_source_factor(box, nc, elec_dens, fields, s_dt, source_factor)
    use m_gas
    use m_transport_data
    use m_lookup_table
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: nc
    real(dp), intent(in)       :: elec_dens(nc**NDIM)
    real(dp), intent(in)       :: fields(nc**NDIM)
    integer, intent(in)        :: s_dt
    real(dp), intent(out)      :: source_factor(nc**NDIM)
    real(dp)                   :: mobilities(nc**NDIM), diffc(nc**NDIM)
    real(dp)                   :: tmp, inv_dr(NDIM), f(2*NDIM)
    real(dp)                   :: Evec(NDIM), gradn(NDIM)
    real(dp), parameter        :: small_flux = 1.0e-9_dp ! A small flux
    real(dp), parameter        :: harmonic_factor = 2.0_dp
    integer                    :: ix, IJK

    if (.not. gas_constant_density) &
         error stop "source_factor: gas_constant_density is false"

    inv_dr = 1/box%dr
    tmp = 1 / gas_number_density
    mobilities = LT_get_col(td_tbl, td_mobility, fields) * tmp

    select case (ST_source_factor)
    case (source_factor_flux)
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
       source_factor = (source_factor + small_flux) / (small_flux + &
            elec_dens * mobilities * &
            pack(box%cc(DTIMES(1:nc), i_electric_fld), .true.))
    case (source_factor_flux_hmean)
       ix = 0
       do KJI_DO(1,nc)
          ix = ix + 1

          ! Compute norm of flux at cell center, but taking harmonic mean of
          ! flux components
#if NDIM == 1
          f = abs([box%fc(i, 1, flux_elec), box%fc(i+1, 1, flux_elec)])
          source_factor(ix) = 2*f(1)*f(2)/(f(1)+f(2)+small_flux)
#elif NDIM == 2
          f = abs([box%fc(i, j, 1, flux_elec), box%fc(i+1, j, 1, flux_elec), &
               box%fc(i, j, 2, flux_elec), box%fc(i, j+1, 2, flux_elec)])
          source_factor(ix) = norm2([2*f(1)*f(2)/(f(1)+f(2)+small_flux), &
               2*f(3)*f(4)/(f(3)+f(4)+small_flux)])
#elif NDIM == 3
          f = abs([&
               box%fc(i, j, k, 1, flux_elec), box%fc(i+1, j, k, 1, flux_elec), &
               box%fc(i, j, k, 2, flux_elec), box%fc(i, j+1, k, 2, flux_elec), &
               box%fc(i, j, k, 3, flux_elec), box%fc(i, j, k+1, 3, flux_elec)])
          source_factor(ix) = norm2([&
               2*f(1)*f(2)/(f(1)+f(2)+small_flux), &
               2*f(3)*f(4)/(f(3)+f(4)+small_flux), &
               2*f(5)*f(6)/(f(5)+f(6)+small_flux)])
#endif
       end do; CLOSE_DO

       ! Compute source factor as harmonic_factor * |flux|/(n_e * mu * E). The
       ! harmonic_factor is used to ensure 'regular' solutions are not affected;
       ! only when source_factor < 1/harmonic_factor does it start to affect the
       ! solution.
       source_factor = harmonic_factor * (source_factor + small_flux) / &
            (small_flux + elec_dens * mobilities * &
            pack(box%cc(DTIMES(1:nc), i_electric_fld), .true.))
    case (source_factor_original_cc)
       ! This is the 'original' scheme, 1 - (E_hat . F_diff)/F_flux
       diffc = LT_get_col(td_tbl, td_diffusion, fields) * tmp
       ix = 0
       do KJI_DO(1,nc)
          ix = ix + 1
#if NDIM == 1
          Evec = 0.5_dp * [box%fc(i, 1, electric_fld) + &
               box%fc(i+1, 1, electric_fld)]
          gradn = 0.5_dp * inv_dr * [&
               box%cc(i+1, i_electron+s_dt) &
               - box%cc(i-1, i_electron+s_dt)]
#elif NDIM == 2
          Evec = 0.5_dp * [box%fc(i, j, 1, electric_fld) + &
               box%fc(i+1, j, 1, electric_fld), &
               box%fc(i, j, 2, electric_fld) + &
               box%fc(i, j+1, 2, electric_fld)]
          gradn = 0.5_dp * inv_dr * [&
               box%cc(i+1, j, i_electron+s_dt) &
               - box%cc(i-1, j, i_electron+s_dt), &
               box%cc(i, j+1, i_electron+s_dt) &
               - box%cc(i, j-1, i_electron+s_dt)]
#elif NDIM == 3
          Evec = 0.5_dp * [box%fc(i, j, k, 1, electric_fld) + &
               box%fc(i+1, j, k, 1, electric_fld), &
               box%fc(i, j, k, 2, electric_fld) + &
               box%fc(i, j+1, k, 2, electric_fld), &
               box%fc(i, j, k, 3, electric_fld) + &
               box%fc(i, j, k+1, 3, electric_fld)]
          gradn = 0.5_dp * inv_dr * [&
               box%cc(i+1, j, k, i_electron+s_dt) &
               - box%cc(i-1, j, k, i_electron+s_dt), &
               box%cc(i, j+1, k, i_electron+s_dt) &
               - box%cc(i, j-1, k, i_electron+s_dt), &
               box%cc(i, j, k+1, i_electron+s_dt) &
               - box%cc(i, j, k-1, i_electron+s_dt)]
#endif
          source_factor(ix) = 1 + diffc(ix) * sum(Evec * gradn) / &
               (small_flux + elec_dens(ix) * mobilities(ix) * &
               box%cc(IJK, i_electric_fld)**2)
       end do; CLOSE_DO
    case (source_factor_original_flux)
       ! Compute source factor as 1 - (E . F_diff)/F_drift
       source_factor = 1 - pack(box%cc(DTIMES(1:nc), i_srcfac), .true.) / &
            (small_flux + elec_dens * mobilities * &
            pack(box%cc(DTIMES(1:nc), i_electric_fld)**2, .true.))
    case default
       error stop
    end select

    source_factor = min(1.0_dp, source_factor)
    source_factor = max(0.0_dp, source_factor)
  end subroutine compute_source_factor

  !> Handle secondary emission from positive ions
  subroutine handle_ion_se_flux(box)
    use m_transport_data
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

  !> Volume integrate chemical reaction rates
  subroutine chemical_rates_box(box, nc, rates, box_rates)
    use m_chemistry
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nc
    real(dp), intent(in)    :: rates(nc**NDIM, n_reactions)
    real(dp), intent(out)   :: box_rates(n_reactions)
#if NDIM == 2
    integer                 :: i, n
    real(dp)                :: w(nc), tmp(nc, nc)
#endif

    if (box%coord_t == af_xyz) then
       box_rates = sum(rates, dim=1) * product(box%dr)
#if NDIM == 2
    else if (box%coord_t == af_cyl) then
       box_rates(:) = 0

       ! Get volume versus radius
       do i = 1, nc
          w(i) = af_cyl_volume_cc(box, i)
       end do

       do n = 1, n_reactions
          tmp = reshape(rates(:, n), [nc, nc])
          do i = 1, nc
             tmp(i, :) = w(i) * tmp(i, :)
          end do
          box_rates(n) = box_rates(n) + sum(tmp)
       end do
#endif
    else
       error stop "Unknown box coordinates"
    end if
  end subroutine chemical_rates_box

  !> Integrate J.E over space into global storage
  subroutine sum_global_JdotE(box, tid)
    use m_units_constants
    type(box_t), intent(in) :: box
    integer, intent(in)     :: tid !< Thread id
    integer                 :: IJK, nc
    real(dp)                :: JdotE, tmp
    real(dp)                :: volume(box%n_cell)

    JdotE = 0.0_dp

    volume = product(box%dr)

#if NDIM == 2
    if (box%coord_t == af_cyl) then
       ! Cylindrical case
       do i = 1, box%n_cell
          volume(i) = af_cyl_volume_cc(box, i)
       end do
    end if
#endif

    nc = box%n_cell
    do KJI_DO(1, nc)
       ! Compute inner product flux * field over the cell faces
       tmp = 0.5_dp * sum(box%fc(IJK, :, flux_elec) * box%fc(IJK, :, electric_fld))
#if NDIM == 1
       tmp = tmp + 0.5_dp * (&
            box%fc(i+1, 1, flux_elec) * box%fc(i+1, 1, electric_fld))
#elif NDIM == 2
       tmp = tmp + 0.5_dp * (&
            box%fc(i+1, j, 1, flux_elec) * box%fc(i+1, j, 1, electric_fld) + &
            box%fc(i, j+1, 2, flux_elec) * box%fc(i, j+1, 2, electric_fld))
#elif NDIM == 3
       tmp = tmp + 0.5_dp * (&
            box%fc(i+1, j, k, 1, flux_elec) * box%fc(i+1, j, k, 1, electric_fld) + &
            box%fc(i, j+1, k, 2, flux_elec) * box%fc(i, j+1, k, 2, electric_fld) + &
            box%fc(i, j, k+1, 3, flux_elec) * box%fc(i, j, k+1, 3, electric_fld))
#endif
       JdotE = JdotE + tmp * volume(i)
    end do; CLOSE_DO

    ST_current_JdotE(1, tid) = ST_current_JdotE(1, tid) + &
         JdotE * UC_elec_charge
  end subroutine sum_global_JdotE

end module m_fluid_lfa
