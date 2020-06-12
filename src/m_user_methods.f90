#include "../afivo/src/cpp_macros.h"
!> This module contains all the methods that users can customize
module m_user_methods
  use m_af_all

  implicit none
  public

  !> User-defined refinement routine
  procedure(af_subr_ref), pointer :: user_refine => null()

  !> If defined, call this routine after setting initial conditions
  procedure(af_subr), pointer :: user_initial_conditions => null()

  !> To set custom boundary conditions for the electric potential
  procedure(af_subr_bc), pointer :: user_potential_bc => null()

  !> To set a user-defined gas number density
  procedure(gas_dens_func), pointer :: user_gas_density => null()

  !> To set the field amplitude manually
  procedure(field_func), pointer :: user_field_amplitude => null()

  procedure(log_subr), pointer :: user_write_log => null()

  interface
     subroutine log_subr(tree, filename, out_cnt)
       import
       type(af_t), intent(in)      :: tree
       character(len=*), intent(in) :: filename
       integer, intent(in)          :: out_cnt
     end subroutine log_subr

     function gas_dens_func(box, IJK) result(dens)
       import
       type(box_t), intent(in) :: box
       integer, intent(in)     :: IJK
       real(dp)                :: dens
     end function gas_dens_func

     function field_func(tree, time) result(amplitude)
       import
       type(af_t), intent(in) :: tree
       real(dp), intent(in)   :: time
       real(dp)               :: amplitude
     end function field_func
  end interface

end module m_user_methods
