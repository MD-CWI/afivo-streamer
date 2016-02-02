module test_mod_omp
  
  implicit none
  private
  
  integer, parameter :: dp = kind(0.0d0)
  
  ! Public methods
  public :: test_2
  
contains
  
    subroutine test_2(rr)
    real(dp), intent(inout) :: rr(:)
    integer                 :: n

    !$omp do
    do n = 1, size(rr)
       ! rr(n) = rr(n) * sin(rr(n)) + cos(rr(n))
       rr(n) = rr(n) * 0.5_dp + 1.0_dp / (rr(n) + abs(rr(n) + 1.0_dp))
    end do
    !$omp end do
  end subroutine test_2
  
end module test_mod_omp