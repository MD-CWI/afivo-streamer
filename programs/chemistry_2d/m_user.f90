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

  ! Whether field waveform versus time is enabled
  logical, protected, public :: field_waveform_enabled = .false.

  real(dp) :: rise_time = 1.0e-9_dp
  real(dp) :: fall_time = 1.0e-9_dp
  real(dp) :: pulse_field = -1.4e6_dp
  real(dp) :: turn_off_z = 18.0e-3_dp
  real(dp) :: detection_density = 1e18_dp
  
contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    call CFG_add_get(cfg, "field_waveform%enabled", field_waveform_enabled, &
          "Whether pulse field versus time is enabled")
    call CFG_add_get(cfg, "field_waveform%rise_time", rise_time, &
          "Linear rise time of field (s)")
    call CFG_add_get(cfg, "field_waveform%fall_time", fall_time, &
          "Linear fall time of field (s)")
    call CFG_add_get(cfg, "field_waveform%pulse_field", pulse_field, &
          "Pulse field (V/m)")
    call CFG_add_get(cfg, "field_waveform%turn_off_z", turn_off_z, &
          "Z coordinate at which turn off the voltage (m)")
    call CFG_add_get(cfg, "field_waveform%detection_density", detection_density, &
         "Detection density (1/m3)")

    if (field_waveform_enabled) then
      user_field_amplitude => my_field_amplitude
    end if

  end subroutine user_initialize

  ! Set electric field amplitude
  function my_field_amplitude(tree, time) result(amplitude)
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: time
    real(dp)               :: amplitude, zmin
    real(dp), save         :: voltage_off_time, prev_time

    if (time <= rise_time) then
      amplitude = pulse_field / rise_time * time
      voltage_off_time = time
      prev_time = time
    else 
      ! Get lowest z-coordinate of a streamer
      call af_reduction(tree, streamer_min_z, reduce_min, 1e100_dp, zmin)
      if (zmin > turn_off_z .and. voltage_off_time >= prev_time) then
        amplitude = pulse_field
        voltage_off_time = time
      else if(time <= (voltage_off_time + fall_time)) then
        amplitude = pulse_field - pulse_field / fall_time * (time - voltage_off_time)
      else 
        amplitude = 0.0_dp
      end if
      prev_time = time
    end if
    
  end function my_field_amplitude

  ! Find cell with smallest z coordinate that has a density exceeding
  ! detection_density
  real(dp) function streamer_min_z(box)
    type(box_t), intent(in) :: box
    integer                 :: i, n, nc
    real(dp)                :: r(NDIM)

    nc = box%n_cell
    i = -1

    do n = 1, nc
       if (maxval(box%cc(1:nc, n, i_electron)) > detection_density) then
          i = n
          exit
       end if
    end do

    if (i /= -1) then
       r = af_r_cc(box, [1, i])
       streamer_min_z = r(NDIM)
    else
       streamer_min_z = 1e100_dp
    end if
  end function streamer_min_z

  real(dp) function reduce_min(a, b)
    real(dp), intent(in) :: a, b
    reduce_min = min(a, b)
  end function reduce_min
  
end module m_user
