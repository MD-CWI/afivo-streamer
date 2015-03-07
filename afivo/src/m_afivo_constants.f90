!> Module with constants used in AFiVO that are relevant for users
module m_afivo_constants

  implicit none
  public

  ! Be aware: the code assumes that rm < kp < do
  !> Value indicating you want to derefine a box
  integer, parameter :: a5_rm_ref = -1

  !> Value indicating you want to keep a box's refinement
  integer, parameter :: a5_kp_ref = 0

  !> Value indicating you want to refine a box
  integer, parameter :: a5_do_ref = 1

  !> The children of a box are removed (for internal use)
  integer, parameter :: a5_derefine = -2

  !> A box will be refined (for internal use)
  integer, parameter :: a5_refine = 2

  !> Special value indicating there is no box
  integer, parameter :: a5_no_box = 0

  ! Each box contains a tag, for which bits can be set.

  !> Bit that indicates that the box is in use
  integer, parameter :: a5_bit_in_use = 1

end module m_afivo_constants
