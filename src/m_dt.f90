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

  !> CFL number to use
  real(dp), public, protected :: dt_cfl_number = undefined_real

  !> If > 0, a density to control the accuracy of the chemistry time step
  real(dp), public, protected :: dt_chemistry_nmin = -1.0_dp

  ! Maximum allowed time step
  real(dp), public, protected :: dt_max = 1.0e-11_dp

  ! Minimum allowed time step
  real(dp), public, protected :: dt_min = 1.0e-14_dp

  !> Maximal relative increase dt for the next iteration
  real(dp), public, protected :: dt_max_growth_factor = 2.0_dp

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
    real(dp)                   :: default_cfl_number = 0.5_dp
    character(len=name_len)    :: integrator

    !> [relevant_parameters]
    call CFG_add_get(cfg, "dt_max", dt_max, &
         "The maximum timestep (s)")
    call CFG_add_get(cfg, "dt_min", dt_min, &
         "The minimum timestep (s)")
    call CFG_add_get(cfg, "dt_safety_factor", dt_safety_factor, &
         "Safety factor for the time step")
    call CFG_add_get(cfg, "dt_cfl_number", dt_cfl_number, &
         "CFL number to use")
    call CFG_add_get(cfg, "dt_chemistry_nmin", dt_chemistry_nmin, &
         "If > 0, a density to control the accuracy of the chemistry time step")
    !> [relevant_parameters]

    call CFG_add_get(cfg, "dt_max_growth_factor", dt_max_growth_factor, &
         "Maximal relative increase dt for the next iteration")

    integrator = "heuns_method"
    call CFG_add_get(cfg, "time_integrator", integrator, &
         "Time integrator (use arbitrary value to see options)")

    do time_integrator = 1, af_num_integrators
       if (integrator == af_integrator_names(time_integrator)) exit
    end do

    if (time_integrator == af_num_integrators+1) then
       print *, "Use one of the following time integrators:"
       do time_integrator = 1, af_num_integrators
          print *, trim(af_integrator_names(time_integrator))
       end do
       error stop "Unknown time integrator"
    end if

    ! Set CFL number automatically if not set
    if (dt_cfl_number <= undefined_real) dt_cfl_number = default_cfl_number

    n_threads = af_get_max_threads()
    ! Prevent cache invalidation issues by enlarging the array
    allocate(dt_matrix(dt_num_cond+32, n_threads))
    dt_matrix(:, :) = 0.0_dp

  end subroutine dt_initialize

end module m_dt
