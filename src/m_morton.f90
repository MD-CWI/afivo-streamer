! This module contains methods to convert indices to morton numbers.
!
! Because fortran does not support unsigned integers, you can only use these
! routines for positive integers.
module m_morton

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: k15 = selected_int_kind(15)

  type morton_t
     ! 64 bit integer
     integer(k15) :: mrtn
  end type morton_t

  ! Public types
  public :: k15
  public :: morton_t

  ! Public methods
  public :: morton_from_ix2
  public :: morton_to_ix2
  public :: morton_bsearch
  public :: print_bits
  public :: print_bits_k15

contains

  function morton_from_ix2(ix) result(m_ix)
    integer, intent(in) :: ix(2)
    type(morton_t)      :: m_ix
    integer(k15) :: a, b
    a = bit_space_1(ix(1))
    b = bit_space_1(ix(2))
    m_ix%mrtn = ior(a, ishft(b, 1))
  end function morton_from_ix2

  function morton_to_ix2(m_ix) result(ix)
    type(morton_t), intent(in) :: m_ix
    integer                    :: ix(2)
    ix(1) = bit_despace_1(m_ix%mrtn)
    ix(2) = bit_despace_1(ishft(m_ix%mrtn, -1))
  end function morton_to_ix2

  function morton_bsearch(list, val) result(ix)
    integer(kind=k15), intent(in) :: list(:)
    integer(kind=k15), intent(in) :: val
    integer                       :: ix, i_min, i_max, i_middle

    i_min = 1
    i_max = size(list)

    do while (i_min < i_max)
       i_middle = i_min + (i_max - i_min) / 2

       if (val <= list(i_middle)) then
          i_max = i_middle
       else
          i_min = i_middle + 1
       end if
    end do

    ix = i_min
    if (val > list(ix)) ix = -1
  end function morton_bsearch

  ! Add two "zero" bits between each bit of the input. Because the result has 64
  ! bits available, only the first 21 bits from the input are spaced.
  function bit_space_2(a) result(x)
    integer, intent(in) :: a
    integer(kind=k15)   :: x

    ! We only look at the first 21 bits
    x = iand(a,                    int(z'1fffff', k15))
    x = iand(ior(x, ishft(x, 32)), int(z'1f00000000ffff', k15))
    x = iand(ior(x, ishft(x, 16)), int(z'1f0000ff0000ff', k15))
    x = iand(ior(x, ishft(x, 8)),  int(z'100f00f00f00f00f', k15))
    x = iand(ior(x, ishft(x, 4)),  int(z'10c30c30c30c30c3', k15))
    x = iand(ior(x, ishft(x, 2)),  int(z'1249249249249249', k15))
  end function bit_space_2

  ! Invert bit_space_2
  function bit_despace_2(a) result(y)
    integer(k15), intent(in) :: a
    integer(k15)             :: x
    integer                  :: y

    x = iand(a,                     int(z'1249249249249249', k15))
    x = iand(ior(x, ishft(x, -2)),  int(z'10c30c30c30c30c3', k15))
    x = iand(ior(x, ishft(x, -4)),  int(z'100f00f00f00f00f', k15))
    x = iand(ior(x, ishft(x, -8)),  int(z'1f0000ff0000ff', k15))
    x = iand(ior(x, ishft(x, -16)), int(z'1f00000000ffff', k15))
    x = iand(ior(x, ishft(x, -32)), int(z'1fffff', k15))
    y = int(x)
  end function bit_despace_2

  ! Add one "zero" bit between each bit of the input.
  function bit_space_1(a) result(x)
    integer, intent(in) :: a
    integer(kind=k15)   :: x

    x = a
    x = iand(ior(x, ishft(x, 32)), int(z'00000000ffffffff', k15))
    x = iand(ior(x, ishft(x, 16)), int(z'0000ffff0000ffff', k15))
    x = iand(ior(x, ishft(x, 8)),  int(z'00ff00ff00ff00ff', k15))
    x = iand(ior(x, ishft(x, 4)),  int(z'0f0f0f0f0f0f0f0f', k15))
    x = iand(ior(x, ishft(x, 2)),  int(z'3333333333333333', k15))
    x = iand(ior(x, ishft(x, 1)),  int(z'5555555555555555', k15))
  end function bit_space_1

  ! Invert bit_space_1
  function bit_despace_1(a) result(y)
    integer(k15), intent(in) :: a
    integer(k15)             :: x
    integer                  :: y

    x = a
    x = iand(x,                     int(z'5555555555555555', k15))
    x = iand(ior(x, ishft(x, -1)),  int(z'3333333333333333', k15))
    x = iand(ior(x, ishft(x, -2)),  int(z'0f0f0f0f0f0f0f0f', k15))
    x = iand(ior(x, ishft(x, -4)),  int(z'00ff00ff00ff00ff', k15))
    x = iand(ior(x, ishft(x, -8)),  int(z'0000ffff0000ffff', k15))
    x = iand(ior(x, ishft(x, -16)), int(z'00000000ffffffff', k15))
    y = int(x)
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