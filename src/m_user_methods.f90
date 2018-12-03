!> This module contains all the methods that users can customize
module m_user_methods
  use m_af_all

  implicit none
  public

  !> User-defined refinement routine
  procedure(af_subr_ref), pointer :: user_refine => null()

  !> If defined, call this routine after setting initial conditions
  procedure(af_subr), pointer :: user_initial_conditions => null()

end module m_user_methods
