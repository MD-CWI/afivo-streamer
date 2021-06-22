!> Template for user code, this one simply uses the default routines
module m_user
  use m_af_all
  use m_config
  use m_user_methods
  use m_streamer

  implicit none
  private

  ! Public methods
  public :: user_initialize

  real(dp) :: initial_field = -1e6_dp
  real(dp) :: v_desired = 5e5_dp
  real(dp) :: adjust_start_position = 0.85_dp
  real(dp) :: control_factor = 1.0e+08_dp
  logical  :: velocity_control    = .false.

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree
    call CFG_add_get(cfg, "my%initial_field", initial_field, &
         "Initial field before adjusting (V/m)")
    call CFG_add_get(cfg, "my%desired_velocity", v_desired, &
         "Desired velocity (m/s)")
    call CFG_add_get(cfg, "my%adjust_start_position", adjust_start_position, &
         "Adjust start position (m)")
    call CFG_add_get(cfg, "my%control_factor", control_factor, &
         "Adjust control factor")
    call CFG_add_get(cfg, "my%velocity_control", velocity_control, &
         "Whether to control the velocity")

    user_initial_conditions => my_init_cond
    if (velocity_control) then
      user_field_amplitude => my_field_amplitude
    end if

  end subroutine user_initialize

  ! To control the streamer velocity
  function my_field_amplitude(tree, time, dt) result(amplitude)
    use m_streamer
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: time
    real(dp), intent(in)   :: dt
    real(dp)               :: amplitude, max_field, r_now(2), v_now, dt
    type(af_loc_t)         :: loc_field
    logical, save          :: first_time = .true.
    logical, save          :: second_time = .true.
    real(dp), save         :: r_pre(2)
    real(dp), save         :: v_pre
    real(dp), save         :: E_pre
    real(dp)               :: change_factor

    ! Get location of Emax
    call af_tree_max_cc(tree, i_electric_fld, max_field, loc_field)
    r_now = af_r_loc(tree, loc_field)
    !dt = global_dt

    if (first_time) then
       amplitude = initial_field
       first_time = .false.
       r_pre = r_now
       E_pre = amplitude
    else if (second_time) then
       amplitude = initial_field
       second_time = .false.
       v_now = abs((r_now(2) - r_pre(2))/dt)
       v_pre = v_now
       r_pre = r_now
       E_pre = amplitude
    else if (r_now(2) > adjust_start_position * ST_domain_len(2)) then
       amplitude = initial_field
       v_now = abs((r_now(2) - r_pre(2))/dt)
       v_pre = v_now
       r_pre = r_now
       E_pre = amplitude
    else
       v_now = abs((r_now(2) - r_pre(2))/dt)
       change_factor = abs(v_now - v_desired)/v_desired * control_factor * dt
       if (v_now > v_desired) then
         !if (v_now > v_pre) then
           amplitude = E_pre * (1 - change_factor)
         !else
           !amplitude = E_pre
         !end if
       else if (v_now < v_desired) then
         !if (v_now < v_pre) then
           amplitude = E_pre * (1 + change_factor)
         !else
           amplitude = E_pre
         !end if
       else
         amplitude = E_pre
       end if
       v_pre = v_now
       r_pre = r_now
       E_pre = amplitude
    end if
    !print *, amplitude, v_now, v_pre, change_factor
  end function my_field_amplitude

  subroutine my_init_cond(box)
    type(box_t), intent(inout) :: box

    ! print *, box%ix
  end subroutine my_init_cond

end module m_user
