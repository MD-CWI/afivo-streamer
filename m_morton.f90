! This module contains methods to convert indices to morton numbers.
module m_morton
  
  implicit none
  private
  
  integer, parameter :: dp = kind(0.0d0)

  type morton_t
     ! 64 bit integer
     integer(kind=) :: 
  end type morton_t
  
  ! Public methods
  public :: 
  
contains

  subroutine morton_from_ix2(ix, m_ix)
    integer, intent(in) :: ix(2)
    type(morton_t), intent(out) :: m_ix
    
  end subroutine morton_from_ix2

  subroutine morton_to_ix2(m_ix, ix)
    type(morton_t), intent(in) :: m_ix
    integer, intent(out) :: ix(2)
    
  end subroutine morton_to_ix2
  
end module m_morton