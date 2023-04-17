#include "../afivo/src/cpp_macros.h"
!> This module contains all the methods that users can customize
module m_user_methods
  use m_af_all
  use m_types

  implicit none
  public

  !< [method_list]
  !> User-defined refinement routine
  procedure(af_subr_ref), pointer :: user_refine => null()

  !> If defined, call this routine after setting initial conditions
  procedure(af_subr), pointer :: user_initial_conditions => null()

  !> If defined, call this routine after a new voltage pulse starts
  procedure(af_subr), pointer :: user_new_pulse_conditions => null()

  !> To set custom boundary conditions for the electric potential
  procedure(af_subr_bc), pointer :: user_potential_bc => null()

  !> To set a user-defined gas number density
  procedure(gas_dens_func), pointer :: user_gas_density => null()

  !> To set the field amplitude manually
  procedure(field_func), pointer :: user_field_amplitude => null()

  !> Generic procedure that is called every time step
  procedure(generic_subr), pointer :: user_generic_method => null()

  !> To write a custom log file
  procedure(log_subr), pointer :: user_write_log => null()

  !> To add entries to the log file
  procedure(log_vars), pointer :: user_log_variables => null()

  !> Custom level-set function
  procedure(mg_func_lsf), pointer :: user_lsf => null()

  !> Function to get boundary value for level set function
  procedure(mg_func_lsf), pointer :: user_lsf_bc => null()
  !< [method_list]

  integer, parameter :: user_max_log_vars = 20

  interface
     !< [interface_list]
     subroutine log_subr(tree, filename, out_cnt)
       import
       type(af_t), intent(in)      :: tree
       character(len=*), intent(in) :: filename
       integer, intent(in)          :: out_cnt
     end subroutine log_subr

     subroutine log_vars(tree, n_vars, var_names, var_values)
       import
       type(af_t), intent(in)                 :: tree
       integer, intent(out)                   :: n_vars
       character(len=name_len), intent(inout) :: var_names(user_max_log_vars)
       real(dp), intent(inout)                :: var_values(user_max_log_vars)
     end subroutine log_vars

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

     logical function bool_subr(tree, time)
       import
       type(af_t), intent(in) :: tree
       real(dp), intent(in)   :: time
     end function bool_subr

     subroutine generic_subr(tree, time)
       import
       type(af_t), intent(in) :: tree
       real(dp), intent(in)   :: time
     end subroutine generic_subr
     !< [interface_list]
  end interface

end module m_user_methods
