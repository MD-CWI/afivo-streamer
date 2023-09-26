!> Template for user code, this one simply uses the default routines
module m_user
  use m_af_all
  use m_config
  use m_user_methods

  implicit none
  private

  real(dp) :: desired_current = 0.1_dp
  real(dp) :: relaxation_time = 5e-9_dp

  ! Public methods
  public :: user_initialize

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    user_field_amplitude => my_field_amplitude
    call CFG_add_get(cfg, "user_current", desired_current, &
        "Supplying the desired current")
    call CFG_add_get(cfg, "user_relaxation_time", relaxation_time, &
        "Supplying the relaxation time for current control")
    user_log_variables => add_log_variables

  end subroutine user_initialize

  real(dp) function my_field_amplitude(tree, time)
    use m_streamer
    use m_field
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: time

    real(dp)       :: resistance, goal_voltage, voltage_change, dt
    integer, save  :: counter   = 0
    real(dp), save :: prev_time = 0

    dt = time - prev_time

    if (time < 14e-9_dp) then
       my_field_amplitude = -24e3 / (ST_domain_len(NDIM))
       !my_field_amplitude = current_voltage / (ST_domain_len(NDIM))
       !print *, "Ampli", current_voltage
    else ! Estimate resistance

       resistance = current_voltage/ST_global_displ_current

       goal_voltage = desired_current * resistance
       voltage_change = (goal_voltage - current_voltage) * dt / relaxation_time
       my_field_amplitude = (current_voltage + voltage_change) / (-ST_domain_len(NDIM))

       counter = counter + 1
       if (modulo(counter, 100) == 0) then
          print *, time, ST_global_displ_current, my_field_amplitude, resistance
       end if
    end if

    prev_time = time

  end function my_field_amplitude

  subroutine add_log_variables(tree, n_vars, var_names, var_values)
       use m_streamer
       type(af_t), intent(in)                 :: tree
       integer, intent(out)                   :: n_vars
       character(len=name_len), intent(inout) :: var_names(user_max_log_vars)
       real(dp), intent(inout)                :: var_values(user_max_log_vars)


       n_vars = 2
       var_names(1) = 'current'
       var_values(1) = ST_global_displ_current
       var_names(2) = 'power_deposited'
       call af_tree_sum_cc(tree, i_power_density, var_values(2))
       !var_names(3) = 'Je_x'
       !call af_tree_sum_cc(tree, af_find_cc_variable(tree,"Je_1"), var_values(3))
       !var_names(4) = 'Je_y'
       !call af_tree_sum_cc(tree, af_find_cc_variable(tree,"Je_2"), var_values(4))

  end subroutine add_log_variables

end module m_user
