#include "../afivo/src/cpp_macros.h"
!> Fluid model module
module m_fluid_lfa
  use m_af_all
  use m_streamer

  implicit none
  private

  logical, private :: last_step

  ! Public methods
  public :: forward_euler

contains

  !> Advance fluid model using forward Euler step. If the equation is written as
  !> y' = f(y), the result is: y(s_out) = y(s_prev) + f(y(s_dt)), where the
  !> s_... refer to temporal states.
  subroutine forward_euler(tree, dt, dt_stiff, dt_lim, time, s_deriv, n_prev, &
       s_prev, w_prev, s_out, i_step, n_steps)
    use m_chemistry
    use m_field
    use m_dt
    use m_transport_data
    use m_dielectric
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt             !< Time step
    real(dp), intent(in)      :: dt_stiff       !< Time step for stiff terms (IMEX)
    real(dp), intent(inout)   :: dt_lim         !< Computed time step limit
    real(dp), intent(in)      :: time           !< Current time
    integer, intent(in)       :: s_deriv        !< State to compute derivatives from
    integer, intent(in)       :: n_prev         !< Number of previous states
    integer, intent(in)       :: s_prev(n_prev) !< Previous states
    real(dp), intent(in)      :: w_prev(n_prev) !< Weights of previous states
    integer, intent(in)       :: s_out          !< Output state
    integer, intent(in)       :: i_step         !< Step of the integrator
    integer, intent(in)       :: n_steps        !< Total number of steps
    integer                   :: ix, id_out
    real(dp)                  :: all_dt(2)

    ! Set current rates to zero; they are summed below
    ST_current_rates = 0
    ST_current_JdotE = 0

    last_step = (i_step == n_steps)

    ! Use a shared array to determine maximum time step
    dt_matrix(1:dt_num_cond, :) = dt_max

    ! Since field_compute is called after performing time integration, we don't
    ! have to call it again for the first sub-step of the next iteration
    if (i_step > 1) call field_compute(tree, mg, s_deriv, time, .true.)

    call flux_upwind_tree(tree, flux_num_species, flux_species, s_deriv, &
         flux_variables, 2, all_dt, flux_upwind, flux_direction, &
         flux_dummy_line_modify, af_limiter_koren_t)

    dt_matrix(dt_ix_cfl, 1) = all_dt(1) * dt_cfl_number
    dt_matrix(dt_ix_drt, 1) = all_dt(2)

    if (transport_data_ions%n_mobile_ions > 0 .and. &
         ion_se_yield > 0.0_dp) then
       ! Handle secondary electron emission from ions
       call af_loop_box(tree, handle_ion_se_flux, .true.)
    end if

    call flux_update_densities(tree, dt, n_species-n_gas_species, &
         species_itree(n_gas_species+1:n_species), flux_num_species, &
         flux_species, flux_variables, s_deriv, n_prev, s_prev, &
         w_prev, s_out, add_source_terms, set_box_mask)

    if (ST_use_dielectric) then
       ! Update surface charge and handle photon emission
       ! @todo For parallelization, think about corner cells with two surfaces
       do ix = 1, diel%max_ix
          if (diel%surfaces(ix)%in_use) then
             id_out = diel%surfaces(ix)%id_out

             ! Convert fluxes onto dielectric to surface charge, and handle
             ! secondary emission
             call dielectric_update_surface_charge(tree%boxes(id_out), &
                  diel%surfaces(ix), dt, n_prev, s_prev, w_prev, s_out)

             ! Add secondary emission from photons hitting the surface
             call dielectric_photon_emission(tree%boxes(id_out), &
                  diel%surfaces(ix), dt, s_out)
          end if
       end do
    end if

    if (last_step) then
       dt_lim = minval(dt_matrix(1:dt_num_cond, :))
    end if
  end subroutine forward_euler

  subroutine flux_upwind(nf, n_var, flux_dim, u, flux, cfl_sum, &
       n_other_dt, other_dt, box, line_ix, s_deriv)
    use m_af_flux_schemes
    use m_units_constants
    use m_gas
    use m_lookup_table
    use m_transport_data
    integer, intent(in)     :: nf              !< Number of cell faces
    integer, intent(in)     :: n_var           !< Number of variables
    integer, intent(in)     :: flux_dim        !< In which dimension fluxes are computed
    real(dp), intent(in)    :: u(nf, n_var)    !< Face values
    real(dp), intent(out)   :: flux(nf, n_var) !< Computed fluxes
    !> Terms per cell-center to be added to CFL sum, see flux_upwind_box
    real(dp), intent(out)   :: cfl_sum(nf-1)
    integer, intent(in)     :: n_other_dt !< Number of non-cfl time step restrictions
    real(dp), intent(inout) :: other_dt(n_other_dt) !< Non-cfl time step restrictions
    type(box_t), intent(in) :: box             !< Current box
    integer, intent(in)     :: line_ix(NDIM-1) !< Index of line for dim /= flux_dim
    integer, intent(in)     :: s_deriv        !< State to compute derivatives from

    real(dp) :: E_cc(0:nf)  !< Cell-centered field strengths
    real(dp) :: E_x(nf)     !< Face-centered field components
    real(dp) :: N_gas(0:nf) !< Gas number density at cell centers
    real(dp) :: ne_cc(0:nf) !< Electron density at cell centers
    real(dp) :: v(nf)       !< Velocity at cell faces
    real(dp) :: dc(nf)      !< Diffusion coefficient at cell faces
    real(dp) :: E_face(nf), Td(nf), N_inv(nf)
    real(dp) :: mu(nf), sigma(nf), inv_dx
    integer  :: n, nc

    nc = box%n_cell
    inv_dx = 1/box%dr(flux_dim)

    ! Inside dielectrics, set the flux to zero
    if (ST_use_dielectric) then
       if (box%cc(DTIMES(1), i_eps) > 1.0_dp) then
          flux = 0.0_dp
          return
       end if
    end if

    call flux_get_line_1fc(box, electric_fld, flux_dim, line_ix, E_x)
    call flux_get_line_1cc(box, i_electric_fld, flux_dim, line_ix, E_cc)
    call flux_get_line_1cc(box, i_electron+s_deriv, flux_dim, line_ix, ne_cc)

    if (gas_constant_density) then
       N_inv = 1/gas_number_density
    else
       ! Compute gas number density at cell faces
       call flux_get_line_1cc(box, i_gas_dens, flux_dim, line_ix, N_gas)
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

    ! Combine advective and diffusive flux
    flux(:, 1) = v * u(:, 1) - dc * inv_dx * (ne_cc(1:nc+1) - ne_cc(0:nc))

    ! Used to determine electron CFL time step
    cfl_sum = max(abs(v(2:)), abs(v(:nf-1))) * inv_dx + &
         2 * max(dc(2:), dc(:nf-1))  * inv_dx**2

    ! Electron conductivity
    sigma = mu * u(:, 1)

    ! Ion fluxes (note: ions are slow, so their CFL condition is ignored)
    do n = 2, flux_num_species
       mu = transport_data_ions%mobilities(n-1) * N_inv
       v = flux_species_charge_sign(n) * mu * E_x
       flux(:, n) = v * u(:, n)
       sigma = sigma + mu * u(:, n)
    end do

    ! Dielectric relaxation time
    other_dt(1) = UC_eps0 / (UC_elem_charge * max(maxval(sigma), 1e-100_dp))
  end subroutine flux_upwind

  subroutine flux_direction(box, line_ix, s_deriv, n_var, flux_dim, direction_positive)
    type(box_t), intent(in) :: box             !< Current box
    integer, intent(in)     :: line_ix(NDIM-1) !< Index of line for dim /= flux_dim
    integer, intent(in)     :: s_deriv         !< State to compute derivatives from
    integer, intent(in)     :: flux_dim        !< In which dimension fluxes are computed
    integer, intent(in)     :: n_var           !< Number of variables
    !> True means positive flow (to the "right"), false to the left
    logical, intent(out)    :: direction_positive(box%n_cell+1, n_var)
    real(dp)                :: E_x(box%n_cell+1)
    integer                 :: n

    call flux_get_line_1fc(box, electric_fld, flux_dim, line_ix, E_x)
    do n = 1, n_var
       direction_positive(:, n) = (flux_species_charge_sign(n) * E_x > 0)
    end do
  end subroutine flux_direction

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

  !> Add chemistry and photoionization source terms
  subroutine add_source_terms(box, dt, n_vars, i_cc, s_deriv, s_out, mask)
    use omp_lib
    use m_units_constants
    use m_gas
    use m_chemistry
    use m_photoi
    use m_dt
    use m_lookup_table
    use m_transport_data
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: dt
    integer, intent(in)        :: n_vars
    integer, intent(in)        :: i_cc(n_vars)
    integer, intent(in)        :: s_deriv
    integer, intent(in)        :: s_out
    logical, intent(in)        :: mask(DTIMES(box%n_cell))

    real(dp)                   :: tmp
    real(dp)                   :: rates(box%n_cell**NDIM, n_reactions)
    real(dp)                   :: derivs(box%n_cell**NDIM, n_species)
    real(dp)                   :: dens(box%n_cell**NDIM, n_species)
    real(dp)                   :: fields(box%n_cell**NDIM), box_rates(n_reactions)
    real(dp)                   :: source_factor(box%n_cell**NDIM)
    integer                    :: IJK, ix, nc, n_cells, n, iv
    integer                    :: tid
    real(dp), parameter        :: eps = 1e-100_dp

    nc      = box%n_cell
    n_cells = box%n_cell**NDIM

    ! Skip this routine if there are no cells to update
    if (.not. any(mask)) return

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
         species_itree(n_gas_species+1:n_species)+s_deriv), [n_cells, n_plasma_species])

    ! It is assumed that species densities should be non-negative. When
    ! computing the effect of chemical reactions, this can also help with
    ! stability, see e.g. http://dx.doi.org/10.1088/1749-4699/6/1/015001
    dens = max(dens, 0.0_dp)

    call get_rates(fields, rates, n_cells)

    if (ST_source_factor /= source_factor_none) then
       if (ST_source_factor /= source_factor_none) then
          call compute_source_factor(box, nc, dens(:, ix_electron), &
               fields, s_deriv, source_factor)
       else
          source_factor(:) = 1.0_dp
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

    ! Note that this routine updates its rates argument
    call get_derivatives(dens, rates, derivs, n_cells)

    if (last_step) then
       tid = omp_get_thread_num() + 1

       ! Update chemistry time step. Note that 'dens' is already non-negative.
       if (dt_chemistry_nmin > 0) then
          ! The time step is restricted by both the production and destruction
          ! rate of species
          tmp = minval((dens + dt_chemistry_nmin) / max(abs(derivs), eps))
       else
          ! Prevent negative values due to too much removal of a species
          tmp = minval(max(dens, eps) / max(-derivs, eps))
       end if

       dt_matrix(dt_ix_rates, tid) = min(dt_matrix(dt_ix_rates, tid), tmp)

       ! Keep track of chemical production at last time integration step
       call chemical_rates_box(box, nc, rates, box_rates)

       !> Integrate rates over space and time into global storage
       ST_current_rates(1:n_reactions, tid) = &
            ST_current_rates(1:n_reactions, tid) + box_rates

       ! Keep track of J.E
       call sum_global_JdotE(box, tid)
    end if

    ix = 0
    do KJI_DO(1,nc)
       ix = ix + 1
       if (.not. mask(IJK)) cycle

       if (photoi_enabled) then
          derivs(ix, ix_electron) = derivs(ix, ix_electron) + &
               box%cc(IJK, i_photo)
          derivs(ix, photoi_species_index) = &
               derivs(ix, photoi_species_index) + box%cc(IJK, i_photo)
       end if

       do n = n_gas_species+1, n_species
          iv = species_itree(n)
          box%cc(IJK, iv+s_out) = box%cc(IJK, iv+s_out) + dt * derivs(ix, n)
       end do
    end do; CLOSE_DO

  end subroutine add_source_terms

  !> Set a mask to true in the gas phase, where the solution should be updated
  subroutine set_box_mask(box, mask)
    type(box_t), intent(in) :: box
    logical, intent(out)    :: mask(DTIMES(box%n_cell))
    integer                 :: n, nc, IJK
    real(dp)                :: coords(box%n_cell, NDIM), r(NDIM)
    real(dp), parameter     :: eps = 1e-10_dp

    nc = box%n_cell
    mask = .true.

    ! Do no update chemistry inside electrode
    if (ST_use_electrode) then
       where (box%cc(DTIMES(1:nc), i_lsf) <= 0.0_dp)
          mask = .false.
       end where
    end if

    ! Inside a dielectric, do not update the species densities
    if (ST_use_dielectric) then
       where (abs(box%cc(DTIMES(1:nc), i_eps) - 1) > eps)
          mask = .false.
       end where
    end if

    ! Optionally limit chemistry to a particular region
    if (ST_plasma_region_enabled) then
       ! Compute box coordinates
       do n = 1, NDIM
          coords(:, n) = box%r_min(n) + box%dr(n) * [(i-0.5_dp, i=1,nc)]
       end do

       do KJI_DO(1,nc)
          r(1) = coords(i, 1)
#if NDIM > 1
          r(2) = coords(j, 2)
#endif
#if NDIM > 2
          r(3) = coords(k, 3)
#endif
          if (any(r < ST_plasma_region_rmin) .or. &
               any(r > ST_plasma_region_rmax)) then
             mask(IJK) = .false.
          end if
       end do; CLOSE_DO
    end if

  end subroutine set_box_mask

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
    real(dp)                   :: mobilities(nc**NDIM)
    real(dp)                   :: N_inv(nc**NDIM)
    real(dp)                   :: inv_dr(NDIM)
    real(dp), parameter        :: small_flux      = 1.0e-9_dp ! A small flux
    integer                    :: ix, IJK

    inv_dr = 1/box%dr

    if (gas_constant_density) then
       N_inv = 1 / gas_number_density
    else
       N_inv = pack(1 / box%cc(DTIMES(1:nc), i_gas_dens), .true.)
    end if

    mobilities = LT_get_col(td_tbl, td_mobility, fields) * N_inv

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
    case default
       error stop "This type of source factor not implemented"
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
