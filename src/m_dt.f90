!> Module to set the time step
module m_dt
  use m_af_all

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
  real(dp), public, protected :: dt_chemistry_nmin = 1e12

  ! Maximum allowed time step
  real(dp), public, protected :: dt_max = 1.0e-11_dp

  ! Minimum allowed time step
  real(dp), public, protected :: dt_min = 1.0e-14_dp

  public :: dt_initialize

contains

  !> Initialize the time step module
  subroutine dt_initialize(cfg)
    use m_config
    use omp_lib
    type(CFG_t), intent(inout) :: cfg
    integer                    :: n_threads

    call CFG_add_get(cfg, "dt_max", dt_max, &
         "The maximum timestep (s)")
    call CFG_add_get(cfg, "dt_min", dt_min, &
         "The minimum timestep (s)")
    call CFG_add_get(cfg, "dt_safety_factor", dt_safety_factor, &
         "Safety factor for the time step")
    call CFG_add_get(cfg, "dt_chemistry_nmin", dt_chemistry_nmin, &
         "Small density for the chemistry time step")

    n_threads = af_get_max_threads()
    ! Prevent cache invalidation issues by enlarging the array
    allocate(dt_matrix(dt_num_cond+32, n_threads))

  end subroutine dt_initialize

end module m_dt
