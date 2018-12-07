!> This module contains all the methods that users can customize
module m_user_methods
  use m_af_all

  implicit none
  public

  !> User-defined refinement routine
  procedure(af_subr_ref), pointer :: user_refine => null()

  !> If defined, call this routine after setting initial conditions
  procedure(af_subr), pointer :: user_initial_conditions => null()

  procedure(log_subr), pointer :: user_write_log => null()

  interface
     subroutine log_subr(tree, filename, out_cnt, dir)
       import
       type(af_t), intent(in)      :: tree
       character(len=*), intent(in) :: filename
       integer, intent(in)          :: out_cnt
       character(len=*), intent(in) :: dir
     end subroutine log_subr
  end interface

end module m_user_methods
