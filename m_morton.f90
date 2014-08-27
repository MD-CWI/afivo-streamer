! This module contains methods to convert indices to morton numbers.
module m_morton

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: k15 = selected_int_kind(15)

  type morton_t
     ! 64 bit integer
     integer(k15) :: mrtn
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

  function bit_space_2(a) result(x)
    integer, intent(in) :: a
    integer(kind=k15)   :: x

    ! We only look at the first 21 bits
    x = iand(a,                    z'1fffff')
    x = iand(ior(x, ishft(x, 32)), z'1f00000000ffff')
    x = iand(ior(x, ishft(x, 16)), z'1f0000ff0000ff')
    x = iand(ior(x, ishft(x, 8)),  z'100f00f00f00f00f')
    x = iand(ior(x, ishft(x, 4)),  z'10c30c30c30c30c3')
    x = iand(ior(x, ishft(x, 2)),  z'1249249249249249')
  end function bit_space_2

  function bit_despace_2(a) result(y)
    integer(k15), intent(in) :: a
    integer(k15)             :: x
    integer                  :: y

    x = iand(a,                     z'1249249249249249')
    x = iand(ior(x, ishft(x, -2)),  z'10c30c30c30c30c3')
    x = iand(ior(x, ishft(x, -4)),  z'100f00f00f00f00f')
    x = iand(ior(x, ishft(x, -8)),  z'1f0000ff0000ff')
    x = iand(ior(x, ishft(x, -16)), z'1f00000000ffff')
    x = iand(ior(x, ishft(x, -32)), z'1fffff')
    y = x
  end function bit_despace_2

  function bit_space_1(a) result(x)
    integer, intent(in) :: a
    integer(kind=k15)   :: x

    x = a
    x = iand(ior(x, ishft(x, 32)), z'00000000ffffffff')
    x = iand(ior(x, ishft(x, 16)), z'0000ffff0000ffff')
    x = iand(ior(x, ishft(x, 8)),  z'00ff00ff00ff00ff')
    x = iand(ior(x, ishft(x, 4)),  z'0f0f0f0f0f0f0f0f')
    x = iand(ior(x, ishft(x, 2)),  z'3333333333333333')
    x = iand(ior(x, ishft(x, 1)),  z'5555555555555555')
  end function bit_space_1

  function bit_despace_1(a) result(y)
    integer(k15), intent(in) :: a
    integer(k15)             :: x
    integer                  :: y

    x = a
    x = iand(x,                     z'5555555555555555')
    x = iand(ior(x, ishft(x, -1)),  z'3333333333333333')
    x = iand(ior(x, ishft(x, -2)),  z'0f0f0f0f0f0f0f0f')
    x = iand(ior(x, ishft(x, -4)),  z'00ff00ff00ff00ff')
    x = iand(ior(x, ishft(x, -8)),  z'0000ffff0000ffff')
    x = iand(ior(x, ishft(x, -16)), z'00000000ffffffff')
    y = x
  end function bit_despace_1

  subroutine print_bits(x)
    integer, intent(in) :: x
    integer :: n
    do n = 0, bit_size(x)-1
       if (btest(x, n)) then
          write(*, "(A)", advance="NO") "1"
       else
          write(*, "(A)", advance="NO") "0"
       end if
    end do
    write(*, *) ""
  end subroutine print_bits

  subroutine print_bits_k15(x)
    integer(k15), intent(in) :: x
    integer :: n
    do n = 0, bit_size(x)-1
       if (btest(x, n)) then
          write(*, "(A)", advance="NO") "1"
       else
          write(*, "(A)", advance="NO") "0"
       end if
    end do
    write(*, *) ""
  end subroutine print_bits_k15
  
end module m_morton