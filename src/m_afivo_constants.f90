module m_afivo_constants

  implicit none
  public

  ! Tags that can be set for a block
  integer, parameter :: a5_do_ref = 1
  integer, parameter :: a5_rm_ref = 2
  integer, parameter :: a5_kp_ref = 3
  integer, parameter :: a5_rm_children = 4

  ! Special value indicating there is no box
  integer, parameter :: a5_no_box = 0

  ! Each box contains a tag, for which bits can be set.
  integer, parameter :: a5_bit_in_use = 1 ! For internal use
  integer, parameter :: a5_bit_new_children = 2 ! For the user

end module m_afivo_constants