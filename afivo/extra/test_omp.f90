program test_omp
  use omp_lib
  use test_mod_omp

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: nn = 10*1000*1000
  integer :: i, n
  real(dp), allocatable :: test_data(:)

  allocate(test_data(nn))
  test_data = 0

  !$omp parallel private(i, n)
  do i = 1, 1000*1000*1000/nn
     call test_2(test_data)
  end do
  !$omp end parallel

  print *, omp_get_max_threads()
contains

  subroutine test_1(rr)
    real(dp), intent(inout) :: rr(:)
    integer                 :: n

    !$omp do schedule(static)
    do n = 1, size(rr)
       ! rr(n) = rr(n) * sin(rr(n)) + cos(rr(n))
       rr(n) = rr(n) * 0.5_dp + 1.0_dp / (rr(n) + abs(rr(n) + 1.0_dp))
    end do
    !$omp end do
  end subroutine test_1
end program test_omp