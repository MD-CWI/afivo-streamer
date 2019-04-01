#include "../afivo/src/cpp_macros.h"
!> Base module for time integration
module m_advance_base
  implicit none
  public

  integer, protected :: time_integrator
  integer, parameter :: forward_euler_t = 1
  integer, parameter :: heuns_method_t  = 2
  integer, parameter :: rk2_t           = 3

  !> How many time states are required
  integer, public, protected :: advance_num_states

  ! Public methods
  public :: advance_base_initialize

contains

  subroutine advance_base_initialize(cfg)
    use m_config
    use m_types
    type(CFG_t), intent(inout) :: cfg
    character(len=name_len)    :: integrator

    integrator = "heuns_method"
    call CFG_add_get(cfg, "time_integrator", integrator, &
         "Time integrator (forward_euler, heuns_method)")
    select case (integrator)
    case ("forward_euler")
       time_integrator = forward_euler_t
       advance_num_states = 1
    case ("rk2")
       time_integrator = rk2_t
       advance_num_states = 2
    case ("heuns_method")
       time_integrator = heuns_method_t
       advance_num_states = 2
    case default
       print *, "Time integrator: ", trim(integrator)
       error stop "Invalid time integrator"
    end select
  end subroutine advance_base_initialize

end module m_advance_base
