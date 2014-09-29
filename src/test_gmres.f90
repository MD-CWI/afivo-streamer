program test_gmres
  use m_gmres

  implicit none
  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: nn = 1000
  integer, parameter :: max_outer = 1
  integer, parameter :: max_inner = nn
  integer :: i
  real(dp) :: x(nn), ax(nn), rhs(nn)
  real(dp) :: a(nn, nn)

  call random_number(a)
  ! do i = 1, nn
  !    a(i, i) = 0.5_dp * sum(abs(a(:, i)))
  ! end do

  x = 0
  rhs = 1

  call gmr_gmres(x, rhs, my_matrix, max_outer, max_inner, &
       0.0_dp, 0.0_dp)
  call my_matrix(x, ax)
  print *, "norm Ax-rhs", norm2(ax-rhs)

contains

  subroutine my_matrix(x, ax)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: ax(:)
    ! Identity matrix, so A * x = x
    ! ax = 2*x
    ! ax(2) = ax(2) + 0.5_dp * ax(1) - 0.333_dp * ax(3)
    ! ax(3) = 7*x(2) - 5*x(4) + x(3)
    ax = matmul(a, x)
  end subroutine my_matrix

end program test_gmres