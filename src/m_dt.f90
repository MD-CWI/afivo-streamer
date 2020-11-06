!> Module to set the time step
module m_dt
  use m_af_all
  use m_types

  implicit none
  private

  ! Number of time step restrictions
  integer, parameter, public :: dt_num_cond = 4

  ! Array of time step restrictions per thread
  real(dp), allocatable, public :: dt_matrix(:, :)

  ! Index of CFL condition
  integer, parameter, public :: dt_ix_cfl = 1

  ! Index of diffusion time step condition
  integer, parameter, public :: dt_ix_diff = 2

  ! Index of dielectric relaxation time step condition
  integer, parameter, public :: dt_ix_drt = 3

  ! Index of reaction rate time step condition
  integer, parameter, public :: dt_ix_rates = 4

  ! Safety factor for the time step
  real(dp), public, protected :: dt_safety_factor = 0.9_dp

  ! Small density for the chemistry time step
  real(dp), public, protected :: dt_chemistry_nmin = 1e15

  ! Maximum allowed time step
  real(dp), public, protected :: dt_max = 1.0e-11_dp

  ! Minimum allowed time step
  real(dp), public, protected :: dt_min = 1.0e-14_dp

  !> Which time integrator is used
  integer, public, protected :: time_integrator

  public :: dt_initialize

contains

  !> Initialize the time step module
  subroutine dt_initialize(cfg)
    use m_config
    use omp_lib
    type(CFG_t), intent(inout) :: cfg
    integer                    :: n_threads
    character(len=name_len)    :: integrator

    !> [relevant_parameters]
    call CFG_add_get(cfg, "dt_max", dt_max, &
         "The maximum timestep (s)")
    call CFG_add_get(cfg, "dt_min", dt_min, &
         "The minimum timestep (s)")
    call CFG_add_get(cfg, "dt_safety_factor", dt_safety_factor, &
         "Safety factor for the time step")
    call CFG_add_get(cfg, "dt_chemistry_nmin", dt_chemistry_nmin, &
         "Small density for the chemistry time step")
    !> [relevant_parameters]

    integrator = "heuns_method"
    call CFG_add_get(cfg, "time_integrator", integrator, &
         "Time integrator (forward_euler, heuns_method)")
    !> [integrators]
    select case (integrator)
    case ("forward_euler")
       time_integrator = af_forward_euler
    case ("rk2")
       time_integrator = af_midpoint_method
    case ("heuns_method")
       time_integrator = af_heuns_method
    case default
       print *, "Time integrator: ", trim(integrator)
       error stop "Invalid time integrator"
    end select
    !> [integrators]

    n_threads = af_get_max_threads()
    ! Prevent cache invalidation issues by enlarging the array
    allocate(dt_matrix(dt_num_cond+32, n_threads))
    dt_matrix(:, :) = 0.0_dp

  end subroutine dt_initialize

end module m_dt
