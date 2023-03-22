#include "../afivo/src/cpp_macros.h"
!> Module that stores parameters related to the gas
module m_gas
  use m_types
  use m_af_types
  use m_units_constants

  implicit none
  private

  !> Whether the gas has a constant density
  logical, public, protected :: gas_constant_density = .true.

  !> Whether the gas dynamics are simulated
  logical, public, protected :: gas_dynamics = .false.

  ! Pressure of the gas in bar
  real(dp), public, protected :: gas_pressure = 1.0_dp

  ! Gas temperature in Kelvin
  real(dp), public, protected :: gas_temperature = 300.0_dp

  ! List of gas fractions (the last one is 1.0 for the total density)
  real(dp), allocatable, public, protected :: gas_fractions(:)

  ! List of gas densities (the last one is the total density)
  real(dp), allocatable, public, protected :: gas_densities(:)

  ! List of gas components (the last one is M, the total density)
  character(len=comp_len), allocatable, public, protected :: gas_components(:)

  ! Gas number density (1/m3)
  real(dp), public, protected :: gas_number_density

  ! Inverse gas number density (1/m3)
  real(dp), public, protected :: gas_inverse_number_density

  ! Convert V/m to Townsend
  real(dp), parameter, public :: SI_to_Townsend = 1e21_dp

  ! Convert Townsend to V/m
  real(dp), parameter, public :: Townsend_to_SI = 1e-21_dp

  ! Index of the gas number density
  integer, public, protected :: i_gas_dens = -1

  ! Gas mean molecular weight (kg)
  real(dp), public, protected :: gas_molecular_weight = 28.8_dp * UC_atomic_mass

  ! Joule heating efficiency
  real(dp), public, protected :: gas_heating_efficiency  = 1.0_dp

  ! Ratio of heat capacities (polytropic index)
  real(dp), public, protected :: gas_euler_gamma = 1.4_dp

  ! Number of variables for the Euler equations
  integer, parameter, public :: n_vars_euler = 2+NDIM

  ! Indices of the Euler variables
  integer, public, protected :: gas_vars(n_vars_euler) = -1
  ! Indices of the primiteve Euler variables
  integer, public, protected :: gas_prim_vars(n_vars_euler+1) = -1

  ! Indices of the Euler fluxes
  integer, public, protected :: gas_fluxes(n_vars_euler) = -1

  ! Names of the Euler variables
  character(len=name_len), public, protected :: gas_var_names(n_vars_euler)

  ! Indices defining the order of the gas dynamics variables
  integer, parameter, public :: i_rho = 1
#if NDIM == 1
  integer, parameter, public :: i_mom(NDIM) = [2]
#elif NDIM == 2
  integer, parameter, public :: i_mom(NDIM) = [2, 3]
#elif NDIM == 3
  integer, parameter, public :: i_mom(NDIM) = [2, 3, 4]
#endif
  integer, parameter, public :: i_e = i_mom(NDIM) + 1

  public :: gas_initialize
  public :: gas_index
  public :: gas_forward_euler

contains

  !> Initialize this module
  subroutine gas_initialize(tree, cfg)
    use m_config
    use m_units_constants
    use m_user_methods
    use m_dt
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg
    integer                    :: n

    call CFG_add_get(cfg, "gas%dynamics", gas_dynamics, &
         "Whether the gas dynamics are simulated")

    if (gas_dynamics) then
       gas_var_names(i_rho) = "gas_rho"
       gas_var_names(i_mom(1)) = "gas_mom_x"
#if NDIM > 1
       gas_var_names(i_mom(2)) = "gas_mom_y"
#endif
#if NDIM == 3
       gas_var_names(i_mom(3)) = "gas_mom_z"
#endif
       gas_var_names(i_e)   = "gas_e"

       gas_constant_density = .false.
       call af_add_cc_variable(tree, "M", ix=i_gas_dens)

       do n = 1, n_vars_euler
          call af_add_cc_variable(tree, gas_var_names(n), ix=gas_vars(n), &
               n_copies=af_advance_num_steps(time_integrator))
          call af_add_fc_variable(tree, "flux", ix=gas_fluxes(n))
          ! @todo improve boundary conditions?
          call af_set_cc_methods(tree, gas_vars(n), af_bc_neumann_zero)
       end do
       call af_add_cc_variable(tree, "u", ix=gas_prim_vars(i_mom(1)))
#if NDIM > 1
       call af_add_cc_variable(tree, "v", ix=gas_prim_vars(i_mom(2)))
#endif
#if NDIM == 3
       call af_add_cc_variable(tree, "w", ix=gas_prim_vars(i_mom(3)))
#endif
       call af_add_cc_variable(tree, "pressure", ix=gas_prim_vars(i_e))
       call af_add_cc_variable(tree, "temperature", ix=gas_prim_vars(i_e+1))
    else if (associated(user_gas_density)) then
       gas_constant_density = .false.
       call af_add_cc_variable(tree, "M", ix=i_gas_dens)
    else
       gas_constant_density = .true.
    end if

    call CFG_add_get(cfg, "gas%pressure", gas_pressure, &
         "The gas pressure (bar)")
    call CFG_add_get(cfg, "gas%temperature", gas_temperature, &
         "The gas temperature (Kelvin)")
    call CFG_add_get(cfg, "gas%molecular_weight", gas_molecular_weight, &
         "Gas mean molecular weight (kg), for gas dynamics")
    call CFG_add_get(cfg, "gas%heating_efficiency", gas_heating_efficiency, &
         "Joule heating efficiency (between 0.0 and 1.0)")

    ! Ideal gas law
    gas_number_density = 1e5_dp * gas_pressure / &
         (UC_boltzmann_const * gas_temperature)
    gas_inverse_number_density = 1/gas_number_density

    call CFG_add(cfg, "gas%components", ["N2", "O2"], &
         "Gas component names", .true.)
    call CFG_add(cfg, "gas%fractions", [0.8_dp, 0.2_dp], &
         "Gas component fractions", .true.)

    call CFG_get_size(cfg, "gas%components", n)
    allocate(gas_components(n+1))
    allocate(gas_fractions(n+1))
    call CFG_get(cfg, "gas%components", gas_components(1:n))
    call CFG_get(cfg, "gas%fractions", gas_fractions(1:n))

    gas_components(n+1) = "M"
    gas_fractions(n+1)  = 1.0_dp

    if (any(gas_fractions < 0.0_dp)) &
         error stop "gas%fractions has negative value"
    if (abs(sum(gas_fractions(1:n)) - 1.0_dp) > 1e-4_dp) &
         error stop "gas%fractions not normalized"

    gas_densities = gas_fractions * gas_number_density

  end subroutine gas_initialize

  !> A forward-Euler method for the Euler equations
  subroutine gas_forward_euler(tree, dt, dt_stiff, dt_lim, time, s_deriv, n_prev, &
       s_prev, w_prev, s_out, i_step, n_steps)
    use m_af_flux_schemes
    use m_af_limiters
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
    real(dp)                  :: wmax(NDIM)

    call flux_generic_tree(tree, n_vars_euler, gas_vars, s_deriv, &
         gas_fluxes, wmax, max_wavespeed, get_fluxes, &
         flux_dummy_other, to_primitive, to_conservative, af_limiter_vanleer_t)
    if (tree%coord_t == af_cyl) then
       call flux_update_densities(tree, dt, n_vars_euler, gas_vars, gas_fluxes, &
            s_deriv, n_prev, s_prev, w_prev, s_out, add_geometric_source)
    else
       call flux_update_densities(tree, dt, n_vars_euler, gas_vars, gas_fluxes, &
            s_deriv, n_prev, s_prev, w_prev, s_out, flux_dummy_source)
    end if

    ! Compute new time step
    dt_lim = 1.0_dp / sum(wmax/af_lvl_dr(tree, tree%highest_lvl))
  end subroutine gas_forward_euler


  !> Add geometric source term for axisymmetric simulations
  subroutine add_geometric_source(box, dt, n_vars, i_cc, s_deriv, s_out)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: dt
    integer, intent(in)        :: n_vars
    integer, intent(in)        :: i_cc(n_vars)
    integer, intent(in)        :: s_deriv
    integer, intent(in)        :: s_out

#if NDIM == 2
    real(dp)                   :: pressure(DTIMES(box%n_cell))
    real(dp)                   :: inv_radius
    integer                    :: nc, i

    nc = box%n_cell
    pressure = get_pressure(box, s_deriv)

    do i = 1, nc
       inv_radius = 1/af_cyl_radius_cc(box, i)
       box%cc(i, 1:nc, i_cc(i_mom(1))+s_out) = &
            box%cc(i, 1:nc, i_cc(i_mom(1))+s_out) + dt * &
            pressure(i, :) * inv_radius
    end do
#endif
  end subroutine add_geometric_source

#if NDIM == 2
  pure function get_pressure(box, s_in) result(pressure)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: s_in
    real(dp)                :: pressure(DTIMES(box%n_cell))
    integer                 :: nc

    nc = box%n_cell
    pressure = (gas_euler_gamma-1.0_dp) * (&
         box%cc(DTIMES(1:nc), gas_vars(i_e)+s_in) - 0.5_dp * &
         sum(box%cc(DTIMES(1:nc), gas_vars(i_mom)+s_in)**2, dim=NDIM+1) / &
         box%cc(DTIMES(1:nc), gas_vars(i_rho)+s_in))
  end function get_pressure
#endif

  !> Find index of a gas component, return -1 if not found
  elemental integer function gas_index(name)
    character(len=*), intent(in) :: name
    do gas_index = 1, size(gas_components)
       if (gas_components(gas_index) == name) exit
    end do
    if (gas_index == size(gas_components)+1) gas_index = -1
  end function gas_index

  subroutine to_primitive(n_values, n_vars, u)
    integer, intent(in)     :: n_values, n_vars
    real(dp), intent(inout) :: u(n_values, n_vars)
    integer :: idim

    do idim = 1, NDIM
       u(:, i_mom(idim)) = u(:, i_mom(idim))/u(:, i_rho)
    end do

    u(:, i_e) = (gas_euler_gamma-1.0_dp) * (u(:, i_e) - &
         0.5_dp*u(:, i_rho)* sum(u(:, i_mom(:))**2, dim=2))
  end subroutine to_primitive

  subroutine to_conservative(n_values, n_vars, u)
    integer, intent(in)     :: n_values, n_vars
    real(dp), intent(inout) :: u(n_values, n_vars)
    real(dp)                :: kin_en(n_values)
    real(dp)                :: inv_fac
    integer                 :: i

    ! Compute kinetic energy (0.5 * rho * velocity^2)
    kin_en = 0.5_dp * u(:, i_rho) * sum(u(:, i_mom(:))**2, dim=2)

    ! Compute energy from pressure and kinetic energy
    inv_fac = 1/(gas_euler_gamma - 1.0_dp)
    u(:, i_e) = u(:, i_e) * inv_fac + kin_en

    ! Compute momentum from density and velocity components
    do i = 1, NDIM
       u(:, i_mom(i)) = u(:, i_rho) * u(:, i_mom(i))
    end do
  end subroutine to_conservative

  subroutine max_wavespeed(n_values, n_var, flux_dim, u, w)
    integer, intent(in)   :: n_values !< Number of cell faces
    integer, intent(in)   :: n_var    !< Number of variables
    integer, intent(in)   :: flux_dim !< In which dimension fluxes are computed
    real(dp), intent(in)  :: u(n_values, n_var) !< Primitive variables
    real(dp), intent(out) :: w(n_values) !< Maximum speed
    real(dp)              :: sound_speeds(n_values)

    sound_speeds = sqrt(gas_euler_gamma * u(:, i_e) / u(:, i_rho))
    w = sound_speeds + abs(u(:, i_mom(flux_dim)))
  end subroutine max_wavespeed

  subroutine get_fluxes(n_values, n_var, flux_dim, u, flux, box, line_ix, s_deriv)
    integer, intent(in)     :: n_values !< Number of cell faces
    integer, intent(in)     :: n_var    !< Number of variables
    integer, intent(in)     :: flux_dim !< In which dimension fluxes are computed
    real(dp), intent(in)    :: u(n_values, n_var)
    real(dp), intent(out)   :: flux(n_values, n_var)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: line_ix(NDIM-1)
    integer, intent(in)     :: s_deriv        !< State to compute derivatives from
    real(dp)                :: E(n_values), inv_fac
    integer                 :: i

    ! Compute left and right flux for conservative variables from the primitive
    ! reconstructed values.

    ! Density flux
    flux(:, i_rho) = u(:, i_rho) * u(:, i_mom(flux_dim))

    ! Momentum flux
    do i = 1, NDIM
       flux(:,  i_mom(i)) = u(:, i_rho) * &
            u(:, i_mom(i)) * u(:, i_mom(flux_dim))
    end do

    ! Add pressure term
    flux(:, i_mom(flux_dim)) = flux(:, i_mom(flux_dim)) + u(:, i_e)

    ! Compute energy
    inv_fac = 1/(gas_euler_gamma-1.0_dp)
    E = u(:, i_e) * inv_fac + 0.5_dp * u(:, i_rho) * &
         sum(u(:, i_mom(:))**2, dim=2)

    ! Energy flux
    flux(:, i_e) = u(:, i_mom(flux_dim)) * (E + u(:, i_e))

  end subroutine get_fluxes

end module m_gas
