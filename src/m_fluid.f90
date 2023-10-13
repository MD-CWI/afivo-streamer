#include "../afivo/src/cpp_macros.h"
!> Fluid model module
module m_fluid_lfa
  use m_af_all
  use m_streamer
  use m_model

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

    ! Set current rates to zero; they are summed below
    ST_current_rates = 0
    ST_current_JdotE = 0

    last_step = (i_step == n_steps)

    ! Since field_compute is called after performing time integration, we don't
    ! have to call it again for the first sub-step of the next iteration
    if (i_step > 1) call field_compute(tree, mg, s_deriv, time, .true.)

    call flux_upwind_tree(tree, flux_num_species, flux_species, s_deriv, &
         flux_variables, 2, dt_limits(1:2), flux_upwind, flux_direction, &
         flux_dummy_line_modify, af_limiter_koren_t)

    if (transport_data_ions%n_mobile_ions > 0 .and. &
         ion_se_yield > 0.0_dp) then
       ! Handle secondary electron emission from ions
       call af_loop_box(tree, handle_ion_se_flux, .true.)
    end if

    call flux_update_densities(tree, dt, size(all_densities), &
         all_densities, flux_num_species, &
         flux_species, flux_variables, s_deriv, n_prev, s_prev, &
         w_prev, s_out, add_source_terms, 2, dt_limits(3:4), set_box_mask)

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

    ! Set time step limit
    dt_limits(1) = dt_limits(1) * dt_cfl_number
    dt_lim = min(dt_max, minval(dt_limits))
  end subroutine forward_euler

  !> Compute flux for the fluid model
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

    real(dp), parameter :: five_third = 5/3.0_dp

    real(dp) :: E_cc(0:nf)  !< Cell-centered field strengths
    real(dp) :: E_x(nf)     !< Face-centered field components
    real(dp) :: N_gas(0:nf) !< Gas number density at cell centers
    real(dp) :: ne_cc(0:nf) !< Electron density at cell centers
    real(dp) :: en_cc(0:nf) !< Electron energy density at cell centers
    real(dp) :: v(nf)       !< Velocity at cell faces
    real(dp) :: dc(nf)      !< Diffusion coefficient at cell faces
    real(dp) :: tmp_fc(nf), N_inv(nf)
    real(dp) :: mu(nf), sigma(nf), inv_dx, cfl_factor
    integer  :: n, nc, flux_ix

    nc = box%n_cell
    inv_dx = 1/box%dr(flux_dim)

    ! Inside dielectrics, set the flux to zero
    if (ST_use_dielectric) then
       if (box%cc(DTIMES(1), i_eps) > 1.0_dp) then
          flux = 0.0_dp
          return
       end if
    end if

    if (gas_constant_density) then
       N_inv = 1/gas_number_density
    else
       ! Compute gas number density at cell faces
       call flux_get_line_1cc(box, i_gas_dens, flux_dim, line_ix, N_gas)
       do n = 1, nc+1
          N_inv(n) = 2 / (N_gas(n-1) + N_gas(n))
       end do
    end if

    call flux_get_line_1fc(box, electric_fld, flux_dim, line_ix, E_x)
    call flux_get_line_1cc(box, i_electron+s_deriv, flux_dim, line_ix, ne_cc)

    if (model_has_energy_equation) then
       call flux_get_line_1cc(box, i_electron_energy+s_deriv, &
            flux_dim, line_ix, en_cc)

       ! Get mean electron energies at cell faces
       tmp_fc = mean_electron_energy(u(:, 2), u(:, 1))
       mu = LT_get_col(td_ee_tbl, td_ee_mobility, tmp_fc) * N_inv
       dc = LT_get_col(td_ee_tbl, td_ee_diffusion, tmp_fc) * N_inv
    else
       call flux_get_line_1cc(box, i_electric_fld, flux_dim, line_ix, E_cc)

       ! Compute field strength at cell faces, which is used to compute the
       ! mobility and diffusion coefficient at the interface
       tmp_fc = 0.5_dp * (E_cc(0:nc) + E_cc(1:nc+1)) * SI_to_Townsend * N_inv
       mu = LT_get_col(td_tbl, td_mobility, tmp_fc) * N_inv
       dc = LT_get_col(td_tbl, td_diffusion, tmp_fc) * N_inv
    end if

    ! Compute velocity, -mu accounts for negative electron charge
    v = -mu * E_x

    ! Combine advective and diffusive flux
    flux(:, 1) = v * u(:, 1) - dc * inv_dx * (ne_cc(1:nc+1) - ne_cc(0:nc))

    ! Electron conductivity
    sigma = mu * u(:, 1)

    if (model_has_energy_equation) then
       flux(:, 2) = five_third * (v * u(:, 2) - &
            dc * inv_dx * (en_cc(1:nc+1) - en_cc(0:nc)))
       cfl_factor = five_third
    else
       cfl_factor = 1.0_dp
    end if

    ! Used to determine electron CFL time step
    cfl_sum = cfl_factor * max(abs(v(2:)), abs(v(:nf-1))) * inv_dx + &
         2 * max(dc(2:), dc(:nf-1))  * inv_dx**2

    ! Ion fluxes (note: ions are slow, so their CFL condition is ignored)
    do n = 1, transport_data_ions%n_mobile_ions
       flux_ix = flux_num_electron_vars + n
       mu = transport_data_ions%mobilities(n) * N_inv
       v = flux_species_charge_sign(flux_ix) * mu * E_x
       flux(:, flux_ix) = v * u(:, flux_ix)
       sigma = sigma + mu * u(:, flux_ix)
    end do

    ! Dielectric relaxation time
    other_dt(1) = UC_eps0 / (UC_elem_charge * max(maxval(sigma), 1e-100_dp))
  end subroutine flux_upwind

  !> Determine the direction of fluxes
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

  !> Get inner product of face-centered variables
  pure function fc_inner_product(box, IJK, flux_a, flux_b) result(inprod)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: IJK  !< Face index
    integer, intent(in)     :: flux_a !< Index of first face-centered variable
    integer, intent(in)     :: flux_b !< Index of second face-centered variable
    real(dp)                :: inprod

    inprod = 0.5_dp * sum(box%fc(IJK, :, flux_a) * box%fc(IJK, :, flux_b))
#if NDIM == 1
    inprod = inprod + 0.5_dp * (&
         box%fc(i+1, 1, flux_a) * box%fc(i+1, 1, flux_b))
#elif NDIM == 2
    inprod = inprod + 0.5_dp * (&
         box%fc(i+1, j, 1, flux_a) * box%fc(i+1, j, 1, flux_b) + &
         box%fc(i, j+1, 2, flux_a) * box%fc(i, j+1, 2, flux_b))
#elif NDIM == 3
    inprod = inprod + 0.5_dp * (&
         box%fc(i+1, j, k, 1, flux_a) * box%fc(i+1, j, k, 1, flux_b) + &
         box%fc(i, j+1, k, 2, flux_a) * box%fc(i, j+1, k, 2, flux_b) + &
         box%fc(i, j, k+1, 3, flux_a) * box%fc(i, j, k+1, 3, flux_b))
#endif
  end function fc_inner_product

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
  subroutine add_source_terms(box, dt, n_vars, i_cc, s_deriv, s_out, n_dt, dt_lim, mask)
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
    integer, intent(in)        :: n_dt
    real(dp), intent(inout)    :: dt_lim(n_dt)
    logical, intent(in)        :: mask(DTIMES(box%n_cell))

    real(dp)                   :: tmp, gain, loss_rate
    real(dp)                   :: rates(box%n_cell**NDIM, n_reactions)
    real(dp)                   :: derivs(box%n_cell**NDIM, n_species)
    real(dp)                   :: dens(box%n_cell**NDIM, n_species)
    real(dp)                   :: fields(box%n_cell**NDIM), box_rates(n_reactions)
    real(dp)                   :: source_factor(box%n_cell**NDIM)
    real(dp)                   :: mean_energies(box%n_cell**NDIM)
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

    if (model_has_energy_equation) then
       mean_energies = pack(mean_electron_energy(&
            box%cc(DTIMES(1:nc), i_electron_energy+s_out), &
            box%cc(DTIMES(1:nc), i_electron+s_out)), .true.)
       call get_rates(fields, rates, n_cells, mean_energies)
    else
       mean_energies = 0.0_dp
       call get_rates(fields, rates, n_cells)
    end if

    if (ST_source_factor /= source_factor_none) then
       if (ST_source_factor /= source_factor_none) then
          call compute_source_factor(box, nc, dens(:, ix_electron), &
               fields, s_deriv, source_factor)
       else
          source_factor(:) = 1.0_dp
       end if

       if (ST_source_min_electrons_per_cell > 0) then
          ! Prevent ionization in cells with a low number of electrons. Note
          ! that that the radius is not taken into account for axisymmetric
          ! cases, as this would lead to artifacts.
          where (dens(:, ix_electron) * minval(box%dr)**3 < &
               ST_source_min_electrons_per_cell)
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

       dt_lim(1) = tmp

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

       if (model_has_energy_equation) then
          gain = -fc_inner_product(box, IJK, flux_elec, electric_fld)
          loss_rate = LT_get_col(td_ee_tbl, td_ee_loss, mean_energies(ix))
          box%cc(IJK, i_electron_energy+s_out) = box%cc(IJK, i_electron_energy+s_out) &
               + dt * (gain - loss_rate * box%cc(IJK, i_electron+s_out))
       end if

       do n = n_gas_species+1, n_species
          iv = species_itree(n)
          box%cc(IJK, iv+s_out) = box%cc(IJK, iv+s_out) + dt * derivs(ix, n)
       end do
    end do; CLOSE_DO

    tmp = maxval(mean_energies)
    if (tmp > 0) then
       ! Set time step restriction for energy loss
       dt_lim(2) = tmp/LT_get_col(td_ee_tbl, td_ee_loss, tmp)
    end if

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

  !> Get mean electron energy
  pure elemental real(dp) function mean_electron_energy(n_energy, n_e)
    real(dp), intent(in) :: n_energy, n_e
    mean_electron_energy = n_energy / max(n_e, 1.0_dp)
  end function mean_electron_energy

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

  !> Handle secondary emission from positive ions at the domain walls
  subroutine handle_ion_se_flux(box)
    use m_transport_data
    type(box_t), intent(inout) :: box
    integer                    :: nc, nb, n, ion_flux, flux_ix

    nc = box%n_cell

    ! Return if there is no physical boundary
    if (all(box%neighbors >= af_no_box)) return

    do nb = 1, af_num_neighbors
       ! Check for physical boundary
       if (box%neighbors(nb) < af_no_box) then

          ! Loop over positive ion species
          do n = 1, transport_data_ions%n_mobile_ions
             flux_ix = flux_num_electron_vars + n
             if (flux_species_charge(flux_ix) > 0.0_dp) then
                ion_flux = flux_variables(flux_ix)
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
       tmp = fc_inner_product(box, IJK, flux_elec, electric_fld)
       JdotE = JdotE + tmp * volume(i)
    end do; CLOSE_DO

    ST_current_JdotE(1, tid) = ST_current_JdotE(1, tid) + &
         JdotE * UC_elec_charge
  end subroutine sum_global_JdotE

end module m_fluid_lfa
