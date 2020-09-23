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

  real(dp) :: initial_field = 2.5e6_dp
  real(dp) :: min_field = 1.0e6_dp
  real(dp) :: decay_distance = 10e-3_dp
  real(dp) :: detection_density = 1e17_dp
  real(dp) :: decay_start_z = 30e-3_dp

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    call CFG_add_get(cfg, "my%initial_field", initial_field, &
         "Initial field before any decay (V/m)")
    call CFG_add_get(cfg, "my%min_field", min_field, &
         "Minimal field (V/m)")
    call CFG_add_get(cfg, "my%decay_distance", decay_distance, &
         "Decay distance (m)")
    call CFG_add_get(cfg, "my%decay_start_z", decay_start_z, &
         "Decay starts from this z-coordinate")

    user_field_amplitude => my_field_amplitude

  end subroutine user_initialize

  ! Get the electric field amplitude
  function my_field_amplitude(tree, time) result(amplitude)
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: time
    real(dp)               :: amplitude, zmin, dist

    ! Get lowest z-coordinate of a streamer
    call af_reduction(tree, streamer_min_z, reduce_min, 1e100_dp, zmin)

    dist = max(decay_start_z - zmin, 0.0_dp)

    ! This function can be customized
    amplitude = min_field + (initial_field - min_field) * &
         exp(-dist / decay_distance)

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
