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

  ! Whether pulse field versus time is enabled
  logical, protected, public :: pulse_field_enabled = .false.

  real(dp) :: start_time = 0.0e-9_dp
  real(dp) :: rise_time = 5.0e-9_dp
  real(dp) :: pulse_time = 30.0e-9_dp
  real(dp) :: fall_time = 5.0e-9_dp
  real(dp) :: start_field = -0.0e6_dp
  real(dp) :: pulse_field = -1.5e6_dp
  real(dp) :: post_field = -0.0e6_dp

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    call CFG_add_get(cfg, "field_waveform%enabled", pulse_field_enabled, &
          "Whether pulse field versus time is enabled")
    call CFG_add_get(cfg, "field_waveform%start_time", start_time, &
          "Start time of field (s)")
    call CFG_add_get(cfg, "field_waveform%rise_time", rise_time, &
          "Linear rise time of field (s)")
    call CFG_add_get(cfg, "field_waveform%pulse_time", pulse_time, &
          "Pulse time of field (s)")
    call CFG_add_get(cfg, "field_waveform%fall_time", fall_time, &
          "Linear fall time of field (s)")
    call CFG_add_get(cfg, "field_waveform%start_field", start_field, &
          "Start field (V/m)")
    call CFG_add_get(cfg, "field_waveform%pulse_field", pulse_field, &
          "Pulse field (V/m)")
    call CFG_add_get(cfg, "field_waveform%post_field", post_field, &
          "Post field (V/m)")

    if (pulse_field_enabled) then
      user_field_amplitude => my_field_amplitude
    end if

  end subroutine user_initialize

  ! Set electric field amplitude
  function my_field_amplitude(tree, time) result(amplitude)
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: time
    real(dp)               :: amplitude
  
    if (time <= start_time) then
      amplitude = start_field
    else if (time <= (start_time + rise_time)) then
      amplitude = start_field + (pulse_field - start_field)/rise_time * (time - start_time) 
    else if (time <= (start_time + rise_time + pulse_time)) then
      amplitude = pulse_field
    else if (time <= (start_time + rise_time + pulse_time + fall_time)) then
      amplitude = pulse_field + (post_field - pulse_field)/fall_time * (time - start_time - rise_time - pulse_time) 
    else 
      amplitude = post_field
    end if

  end function my_field_amplitude

end module m_user
