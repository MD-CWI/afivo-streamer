!> Template for user code, this one simply uses the default routines
module m_user
  use m_af_all
  use m_config
  use m_user_methods

  implicit none
  private

  ! Public methods
  public :: user_initialize

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    user_evolve_electrons => electrons_active
  end subroutine user_initialize

  logical function electrons_active(tree, time)
    use m_field
    type(af_t), intent(in) :: tree
    real(dp), intent(in) :: time

    call field_set_voltage(tree, time)
    electrons_active = (abs(field_voltage) > 1e3_dp)

  end function electrons_active

end module m_user
