!> Module for cubic spline interpolation
!>
!> Author: Jannis Teunissen
!>
!> This code is based on:
!> 1. spline.f (https://www.netlib.org/fmm/spline.f)
!> 2. spline.f90 (translation of spline.f,
!> https://ww2.odu.edu/~agodunov/computing/programs/book2/Ch01/spline.f90)
!>
!> License: the above two files do not specify a license, so I assume they are
!> in the public domain, as is this code.
module m_spline_interp
  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  type spline_t
     !> Number of tabulated points
     integer               :: n
     !> Tabulated x values
     real(dp), allocatable :: x(:)
     !> Tabulated y values
     real(dp), allocatable :: y(:)
     !> Spline coefficients such that s(x) = y(i) + b(i)*(x-x(i)) +
     !> c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
     real(dp), allocatable :: bcd(:, :)
     !> Array to quickly get index for a coordinate
     integer, allocatable  :: lookup_index(:)
     !> Inverse spacing in lookup array
     real(dp)              :: lookup_inv_dx
  end type spline_t

  public :: spline_t
  public :: spline_set_coeffs
  public :: spline_evaluate

contains


  !> Calculate coefficients for cubic spline interpolation
  subroutine spline_set_coeffs(x, y, n, spl)
    !> Size of tabulated data
    integer, intent(in)         :: n
    !> Tabulated coordinates (in strictly increasing order)
    real(dp), intent(in)        :: x(n)
    !> Tabulated values
    real(dp), intent(in)        :: y(n)
    type(spline_t), intent(out) :: spl
    real(dp), allocatable       :: b(:), c(:), d(:)
    integer                     :: i, j, lookup_size
    real(dp)                    :: h, x_lookup, min_dx

    if (n < 2) &
         error stop "spline_set_coeffs requires n >= 2"
    if (any(x(2:n) <= x(1:n-1))) &
         error stop "spline_set_coeffs: x(:) not strictly increasing"

    allocate(spl%x(n), spl%y(n), spl%bcd(3, n), b(n), c(n), d(n))
    spl%n = n
    spl%x = x
    spl%y = y

    if (n < 3) then
       ! Handle special case by linear interpolation
       b(1) = (y(2)-y(1))/(x(2)-x(1))
       c(1) = 0
       d(1) = 0
       b(2) = b(1)
       c(2) = 0
       d(2) = 0
       return
    end if

    ! set up tridiagonal system
    ! b = diagonal, d = offdiagonal, c = right hand side.
    d(1) = x(2) - x(1)
    c(2) = (y(2) - y(1))/d(1)
    do i = 2, n-1
       d(i) = x(i+1) - x(i)
       b(i) = 2*(d(i-1) + d(i))
       c(i+1) = (y(i+1) - y(i))/d(i)
       c(i) = c(i+1) - c(i)
    end do

    ! end conditions.  third derivatives at  x(1)  and  x(n)
    ! obtained from divided differences
    b(1) = -d(1)
    b(n) = -d(n-1)
    c(1) = 0
    c(n) = 0
    if(n /= 3) then
       c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
       c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
       c(1) = c(1)*d(1)**2/(x(4)-x(1))
       c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
    end if

    ! forward elimination
    do i = 2, n
       h = d(i-1)/b(i-1)
       b(i) = b(i) - h*d(i-1)
       c(i) = c(i) - h*c(i-1)
    end do

    ! back substitution
    c(n) = c(n)/b(n)
    do j = 1, n-1
       i = n-j
       c(i) = (c(i) - d(i)*c(i+1))/b(i)
    end do

    ! compute spline coefficients
    b(n) = (y(n) - y(n-1))/d(n-1) + d(n-1)*(c(n-1) + 2*c(n))
    do i = 1, n-1
       b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2*c(i))
       d(i) = (c(i+1) - c(i))/d(i)
       c(i) = 3.*c(i)
    end do
    c(n) = 3*c(n)
    d(n) = d(n-1)

    spl%bcd(1, :) = b
    spl%bcd(2, :) = c
    spl%bcd(3, :) = d

    ! Create linear lookup table to find location between data points more
    ! quickly. First determine good size for the lookup table.
    min_dx      = minval(x(2:n) - x(1:n-1))
    lookup_size = min(4 * n, nint(1 + (x(n) - x(1))/min_dx))

    ! The lookup table will have a regular (linear) spacing
    h = (x(n) - x(1)) / (lookup_size - 1)
    spl%lookup_inv_dx = 1/h
    allocate(spl%lookup_index(lookup_size))

    ! At location z, the table index is i = ceiling((z - x(1)) * inv_dx)
    ! The tabulated points should then be x(spl%lookup_index(i))
    spl%lookup_index(1) = 1
    do i = 2, lookup_size
       x_lookup = x(1) + (i-1) * h
       spl%lookup_index(i) = spl%lookup_index(i-1)
       do while (x_lookup > x(spl%lookup_index(i)+1))
          spl%lookup_index(i) = spl%lookup_index(i) + 1
       end do
    end do

  end subroutine spline_set_coeffs

  ! Evaluate the cubic spline interpolation at point u
  ! result = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
  ! where  x(i) <= u <= x(i+1)
  pure elemental function spline_evaluate(u, spl) result(spline_value)
    !> Evaluate at this coordinate
    real(dp), intent(in)       :: u
    !> Spline data
    type(spline_t), intent(in) :: spl
    integer                    :: i, j
    real(dp)                   :: spline_value, dx

    associate (n=>spl%n, x=>spl%x, y=>spl%y, bcd=>spl%bcd)
      if(u <= x(1)) then
         spline_value = y(1)
         return
      else if(u >= x(n)) then
         spline_value = y(n)
         return
      end if

      j = ceiling(spl%lookup_inv_dx * (u - x(1)))

      ! Even with the checks above, we probably have to check whether j < 1 or
      ! j > size(spl%lookup_index) due to numerical round-off error
      if (j < 1) then
         j = 1
      else if (j > size(spl%lookup_index)) then
         j = size(spl%lookup_index)
      end if

      ! TODO: Not sure whether in practical applications it could be useful to
      ! do a binary search here instead of a linear one

      ! Find index i so that x(i) <= u <= x(i+1)
      do i = spl%lookup_index(j), n-1
         if (u <= x(i+1)) exit
      end do

      ! evaluate spline interpolation
      dx = u - x(i)
      spline_value = y(i) + dx*(bcd(1, i) + dx*(bcd(2, i) + dx*bcd(3, i)))
    end associate
  end function spline_evaluate

end module m_spline_interp
