!> Template for user code, this one simply uses the default routines
module m_user
  use m_af_all
  use m_config
  use m_user_methods

  implicit none
  private

  ! Public methods
  public :: user_initialize

  real(dp) :: initial_field = -2e6_dp
  real(dp) :: min_field = -5e5_dp
  real(dp) :: decay_distance = 10e-3_dp
  real(dp) :: decay_start_time = 10.0e-9_dp

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
    call CFG_add_get(cfg, "my%decay_start_time", decay_start_time, &
         "Decay start time (s)")

    user_field_amplitude => my_field_amplitude
  end subroutine user_initialize

  function my_field_amplitude(tree, time) result(amplitude)
    use m_streamer
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: time
    real(dp)               :: amplitude, max_field, r(NDIM), dist
    type(af_loc_t)         :: loc_field
    logical, save          :: first_time = .true.
    real(dp), save         :: initial_coord

    ! Get location of Emax
    call af_tree_max_cc(tree, i_electric_fld, max_field, loc_field)
    r = af_r_loc(tree, loc_field)

    if (first_time .or. time < decay_start_time) then
       initial_coord = r(NDIM)
       amplitude = initial_field
       first_time = .false.
    else
       dist = abs(r(NDIM) - initial_coord)
       amplitude = min_field + (initial_field - min_field) * &
            exp(-dist / decay_distance)
    end if

    ! print *, time, r(NDIM), amplitude
  end function my_field_amplitude

end module m_user
