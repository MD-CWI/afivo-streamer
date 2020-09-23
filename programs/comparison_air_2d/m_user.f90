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

    user_initial_conditions => my_init_cond

  end subroutine user_initialize

  subroutine my_init_cond(box)
    type(box_t), intent(inout) :: box

    ! print *, box%ix
  end subroutine my_init_cond

end module m_user
