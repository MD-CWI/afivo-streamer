!> Template for user code, this one simply uses the default routines
module m_user
  use m_af_all
  use m_config
  use m_user_methods

  implicit none
  private

  ! Public methods
  public :: user_initialize

  integer, parameter :: buffer_size = 5
  integer :: buffer_index = 0
  real(dp) :: vring(buffer_size) = 0
  real(dp) :: goal_velocity = 3.0e5_dp
  real(dp), parameter :: dfieldt = -2e14_dp

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    user_generic_method => my_velocity
    user_field_amplitude => my_field_amplitude
  end subroutine user_initialize

  real(dp) function my_field_amplitude(tree, time)
    use m_field
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: time
    real(dp)               :: v, diff_amplitude
    real(dp), save         :: prev_field, prev_time

    v = sum(vring)/buffer_size

    if (time < 1e-9_dp) then
       my_field_amplitude = field_amplitude
       prev_field         = field_amplitude
       prev_time          = time
    else
       diff_amplitude     = (goal_velocity - v)/goal_velocity * dfieldt * &
            (time - prev_time)
       my_field_amplitude = prev_field + diff_amplitude
       prev_time          = time
       prev_field         = my_field_amplitude
       print *, time, my_field_amplitude, (goal_velocity - v)/goal_velocity
    end if

  end function my_field_amplitude

  subroutine my_velocity(tree, time)
    use m_streamer
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: time
    real(dp)               :: max_field
    type(af_loc_t)         :: loc_field
    real(dp)               :: xnew(NDIM), n_grid_cells
    logical, save          :: first_time = .true.
    real(dp)               :: v
    real(dp), save         :: prev_time, x_prev(NDIM)

    call af_tree_max_cc(tree, i_electric_fld, max_field, loc_field)

    if (first_time) then
       x_prev = af_r_loc(tree, loc_field)
       prev_time  = time
       first_time = .false.
       buffer_index = 1
    else
       xnew = af_r_loc(tree, loc_field)
       n_grid_cells = abs(xnew(NDIM) - x_prev(NDIM))/af_min_dr(tree)

       if (n_grid_cells > 7.5_dp) then
          v         = abs(xnew(NDIM) - x_prev(NDIM)) / (time - prev_time)
          x_prev    = xnew
          prev_time = time

          buffer_index = buffer_index+1
          if (buffer_index > buffer_size) buffer_index = 1
          vring(buffer_index) = v
       end if
    end if

  end subroutine my_velocity

end module m_user
