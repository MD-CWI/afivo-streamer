program test_omp
  use test_mod_omp
  
  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: nn = 100
  integer :: i
  real(dp), allocatable :: test_data(:)

  allocate(test_data(nn))
  test_data = 0
  
  do i = 1, 100*1000*1000/nn
     call test_1(test_data)
  end do

contains

  subroutine test_1(rr)
    real(dp), intent(inout) :: rr(:)
    integer                 :: n

    !$omp parallel do schedule(static)
    do n = 1, size(rr)
       ! rr(n) = rr(n) * sin(rr(n)) + cos(rr(n))
       rr(n) = rr(n) * 0.5_dp + 1.0_dp / (rr(n) + abs(rr(n)))
    end do
    !$omp end parallel do
  end subroutine test_1
end program test_omp