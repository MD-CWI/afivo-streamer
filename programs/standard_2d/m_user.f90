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

    user_evolve_electrons => gas_active
    
    user_log_variables => total_power_deposited
  end subroutine user_initialize

  logical function electrons_active(tree, time)
    use m_field
    type(af_t), intent(in) :: tree
    real(dp), intent(in) :: time

    call field_set_voltage(tree, time)
    electrons_active = (abs(field_voltage) > 1e3_dp)

  end function electrons_active


  logical function gas_active(tree, time)
    use m_streamer
    type(af_t), intent(in) :: tree
    real(dp), intent(in) :: time
    real(dp)             :: max_field, r(2)
    type(af_loc_t)       :: loc_field

    call af_tree_max_cc(tree, i_electric_fld, max_field, loc_field)
    r = af_r_loc(tree, loc_field)
    gas_active = (abs(1_dp - r(2)/ST_domain_len(2)) < 0.97_dp)

    if (.not.(gas_active)) then
        print *, 'Stopping discharge evolution'
    end if

  end function gas_active
  
  subroutine total_power_deposited(tree, n_vars, var_names, var_values)
       use m_streamer
       type(af_t), intent(in)                 :: tree
       integer, intent(out)                   :: n_vars
       character(len=name_len), intent(inout) :: var_names(user_max_log_vars)
       real(dp), intent(inout)                :: var_values(user_max_log_vars)
       
       
       n_vars = 1
       var_names(1) = 'total_power'
       
       call af_tree_sum_cc(tree, i_power_density, var_values(1))
  end subroutine total_power_deposited
  
  

end module m_user
