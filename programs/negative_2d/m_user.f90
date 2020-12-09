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

  ! Whether field_decay versus length is enabled
  logical, protected, public :: field_decay_enabled = .false.

  ! Decay parameters of the applied electric field
  character(len=32) :: decay_profile = "exponential" ! (constant, exponential, linear_across_zero, linear_over_zero, step)
  real(dp) :: initial_field = 2.5e6_dp
  real(dp) :: min_field =0.8e6_dp
  real(dp) :: decay_distance = 10e-3_dp
  real(dp) :: decay_slope = -0.4e8_dp
  real(dp) :: decay_start_z = 60e-3_dp
  real(dp) :: detection_density = 1e16_dp

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

   call CFG_add_get(cfg, "field_decay%enabled", field_decay_enabled, &
         "Whether field decay versus length is enabled")
   call CFG_add_get(cfg, "field_decay%decay_profile", decay_profile, &
         "Decay type for applied electric field (constant, exponential, linear_across_zero, linear_over_zero, step)")
    call CFG_add_get(cfg, "field_decay%initial_field", initial_field, &
         "Initial electric field before any decay (V/m)")
    call CFG_add_get(cfg, "field_decay%min_field", min_field, &
         "Minimal electric field (V/m)")
    call CFG_add_get(cfg, "field_decay%decay_distance", decay_distance, &
         "Decay distance (m)")
    call CFG_add_get(cfg, "field_decay%decay_slope", decay_slope, &
         "Decay slope (V/m2)")
    call CFG_add_get(cfg, "field_decay%decay_start_z", decay_start_z, &
         "Decay starts from this z-coordinate (m)")
    call CFG_add_get(cfg, "field_decay%detection_density", detection_density, &
         "Detection density (1/m3)")

    if (field_decay_enabled) then
      user_field_amplitude => my_field_amplitude
    end if

  end subroutine user_initialize


  ! Get the electric field amplitude
  function my_field_amplitude(tree, time) result(amplitude)
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: time
    real(dp)               :: amplitude, zmin, dist

    ! Get lowest z-coordinate of a streamer
    call af_reduction(tree, streamer_min_z, reduce_min, 1e100_dp, zmin)

    dist = max(decay_start_z - zmin, 0.0_dp)

  ! Electric field proflie
  select case (decay_profile)
  case ("constant")  ! constant electric field versus length
    amplitude = initial_field
  case ("exponential")  ! exponential decay electric field versus length
    amplitude = min_field + (initial_field - min_field) * exp(-dist / decay_distance)
  case ("linear_across_zero")  ! linear decay electric field versus length (across zero)
    amplitude = initial_field+decay_slope * dist
    if (amplitude<=-initial_field) then
       amplitude = -initial_field
    end if
  case ("linear_over_zero")  ! linear decay electric field versus length (over zero)
    amplitude = initial_field+decay_slope * dist
    if (amplitude<=0) then
       amplitude = 0
    end if
  case ("step")  ! step decay electric field versus length
    if (dist == 0) then
       amplitude = initial_field
    else
       amplitude = min_field
    end if
   end select

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
