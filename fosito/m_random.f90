!> Module for pseudo random number generation
module m_random

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: i4 = selected_int_kind(9)

  type, public :: RNG_t
     integer(i4), private :: state(4) = (/521288629, &
          362436069, 16163801, 1131199299/)
   contains
     procedure, non_overridable :: set_seed
     procedure, non_overridable :: int4
     procedure, non_overridable :: int_ab
     procedure, non_overridable :: uni_01
     procedure, non_overridable :: uni_ab
     procedure, non_overridable :: two_normals
     procedure, non_overridable :: poisson
     procedure, non_overridable :: sphere
  end type RNG_t

contains

  subroutine set_seed(self, seed)
    class(RNG_t), intent(out) :: self
    integer(i4), intent(in)   :: seed(4)
    self%state = seed
  end subroutine set_seed

  ! Method mzran (Marsaglia), see http://jblevins.org/log/openmp
  function int4(self) result(rr)
    class(RNG_t), intent(inout) :: self
    integer(i4)                 :: rr

    rr              = self%state(1) - self%state(3)
    if (rr < 0) rr  = rr + 2147483579
    self%state(1:2) = self%state(2:3)
    self%state(3)   = rr
    self%state(4)   = 69069 * self%state(4) + 1013904243
    rr              = rr + self%state(4)
  end function int4

  function int_ab(self, a, b) result(rr)
    class(RNG_t), intent(inout) :: self
    integer, intent(in)         :: a, b
    integer                     :: rr
    rr = a + int(self%uni_01() * (b-a+1))
  end function int_ab

  ! Uniform random number in range [0,1)
  function uni_01(self) result(rr)
    class(RNG_t), intent(inout) :: self
    real(dp)                    :: rr
    real(dp), parameter         :: conv_fac = 2.0_dp**(-32)
    rr = 0.5_dp + conv_fac * self%int4()
  end function uni_01

  ! Uniform random number in range [a,b)
  function uni_ab(self, a, b) result(rr)
    class(RNG_t), intent(inout) :: self
    real(dp), intent(in)        :: a, b
    real(dp)                    :: rr
    rr = a + self%uni_01() * (b-a)
  end function uni_ab

  ! Return two normal random variates with mean 0 and variance 1
  ! Source: http://en.wikipedia.org/wiki/Marsaglia_polar_method
  function two_normals(self) result(rands)
    class(RNG_t), intent(inout) :: self
    real(dp)                    :: rands(2), sum_sq
    do
       rands(1) = self%uni_ab(-1.0_dp, 1.0_dp)
       rands(2) = self%uni_ab(-1.0_dp, 1.0_dp)
       sum_sq = sum(rands**2)
       if (sum_sq < 1.0_dp .and. sum_sq /= 0.0_dp) exit
    end do
    rands = rands * sqrt(-2 * log(sum_sq) / sum_sq)
  end function two_normals

  ! Return Poisson random variate with rate lambda. Works well for lambda < 30
  ! or so. For lambda >> 1 it can produce wrong results due to roundoff error.
  function poisson(self, lambda) result(rr)
    class(RNG_t), intent(inout) :: self
    real(dp), intent(in)        :: lambda
    integer(i4)                 :: rr
    real(dp)                    :: expl, p

    expl = exp(-lambda)
    rr = 0
    p = self%uni_01()

    do while (p > expl)
       rr = rr + 1
       p = p * self%uni_01()
    end do
  end function poisson

  ! Sample point on a sphere with given radius
  function sphere(self, radius) result(xyz)
    class(RNG_t), intent(inout) :: self
    real(dp), intent(in)        :: radius
    real(dp)                    :: rands(2), xyz(3)
    real(dp)                    :: sum_sq, tmp_sqrt

    ! Marsaglia method for uniform sampling on sphere
    do
       rands(1) = self%uni_ab(-1.0_dp, 1.0_dp)
       rands(2) = self%uni_ab(-1.0_dp, 1.0_dp)
       sum_sq   = rands(1)**2 + rands(2)**2
       if (sum_sq <= 1) exit
    end do

    tmp_sqrt = sqrt(1 - sum_sq)
    xyz(1:2) = 2 * rands(1:2) * tmp_sqrt
    xyz(3)   = 1 - 2 * sum_sq
    xyz      = xyz * radius
  end function sphere

end module m_random
