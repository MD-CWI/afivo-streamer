!> Module with constants used in AFiVO that are relevant for users
module m_afivo_constants

  implicit none
  public

  !> Value indicating a box should be refined
  integer, parameter :: a5_do_ref = 1

  !> Value indicating a box should be derefined
  integer, parameter :: a5_rm_ref = 2

  !> Value indicating a box should keep its refinement
  integer, parameter :: a5_kp_ref = 3

  !> Value indicating that the children of a box can be removed (for internal use)
  integer, parameter :: a5_rm_children = 4

  !> Special value indicating there is no box
  integer, parameter :: a5_no_box = 0

  ! Each box contains a tag, for which bits can be set.

  !> Bit that indicates that the box is in use
  integer, parameter :: a5_bit_in_use = 1

  !> Bit that indicates that a box has gotten new children
  integer, parameter :: a5_bit_new_children = 2

end module m_afivo_constants
