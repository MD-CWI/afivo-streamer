!> Module for pseudo random number generation. The internal pseudo random
!> generator is the xoroshiro128plus method.
module m_random

  implicit none
  private

  ! A 64 bit floating point type
  integer, parameter :: dp = kind(0.0d0)

  ! A 32 bit integer type
  integer, parameter :: i4 = selected_int_kind(9)

  ! A 64 bit integer type
  integer, parameter :: i8 = selected_int_kind(18)

  type rng_t
     !> The rng state (needs initialization)
     integer(i8), private       :: s(2) = [0_i8, 0_i8]
   contains
     procedure, non_overridable :: set_seed    ! Seed the generator
     procedure, non_overridable :: jump        ! Jump function (see below)
     procedure, non_overridable :: int4        ! 4-byte random integer
     procedure, non_overridable :: int8        ! 8-byte random integer
     procedure, non_overridable :: unif_01     ! Uniform (0,1] real
     procedure, non_overridable :: two_normals ! Two normal(0,1) samples
     procedure, non_overridable :: poisson     ! Sample from Poisson-dist.
     procedure, non_overridable :: circle      ! Sample on a circle
     procedure, non_overridable :: sphere      ! Sample on a sphere
     procedure, non_overridable :: next        ! Internal method
  end type rng_t

  public :: rng_t

contains

  !> Set a seed for the rng
  subroutine set_seed(self, the_seed)
    class(rng_t), intent(inout) :: self
    integer(i8), intent(in)     :: the_seed(2)

    self%s = the_seed
  end subroutine set_seed

  ! This is the jump function for the generator. It is equivalent
  ! to 2^64 calls to next(); it can be used to generate 2^64
  ! non-overlapping subsequences for parallel computations.
  subroutine jump(self)
    class(rng_t), intent(inout) :: self
    integer                     :: i, b
    integer(i8)                 :: t(2), dummy

    ! The signed equivalent of the unsigned constants
    integer(i8), parameter      :: jmp_c(2) = &
         (/-4707382666127344949_i8, -2852180941702784734_i8/)

    t = 0
    do i = 1, 2
       do b = 0, 63
          if (iand(jmp_c(i), shiftl(1_i8, b)) /= 0) then
             t = ieor(t, self%s)
          end if
          dummy = self%next()
       end do
    end do

    self%s = t
  end subroutine jump

  !> Return 4-byte integer
  integer(i4) function int4(self)
    class(rng_t), intent(inout) :: self
    int4 = int(self%next(), i4)
  end function int4

  !> Return 8-byte integer
  integer(i8) function int8(self)
    class(rng_t), intent(inout) :: self
    int8 = self%next()
  end function int8

  !> Get a uniform [0,1) random real (double precision)
  real(dp) function unif_01(self)
    class(rng_t), intent(inout) :: self
    integer(i8)                 :: x
    real(dp)                    :: tmp

    x   = self%next()
    x   = ior(shiftl(1023_i8, 52), shiftr(x, 12))
    unif_01 = transfer(x, tmp) - 1.0_dp
  end function unif_01

  !> Return two normal random variates with mean 0 and variance 1.
  !> http://en.wikipedia.org/wiki/Marsaglia_polar_method
  function two_normals(self) result(rands)
    class(rng_t), intent(inout) :: self
    real(dp)                    :: rands(2), sum_sq

    do
       rands(1) = 2 * self%unif_01() - 1
       rands(2) = 2 * self%unif_01() - 1
       sum_sq = sum(rands**2)
       if (sum_sq < 1.0_dp .and. sum_sq > 0.0_dp) exit
    end do
    rands = rands * sqrt(-2 * log(sum_sq) / sum_sq)
  end function two_normals

  !> Return Poisson random variate with rate lambda. Works well for lambda < 30
  !> or so. For lambda >> 1 it can produce wrong results due to roundoff error.
  function poisson(self, lambda) result(rr)
    class(rng_t), intent(inout) :: self
    real(dp), intent(in)        :: lambda
    integer(i4)                 :: rr
    real(dp)                    :: expl, p

    expl = exp(-lambda)
    rr   = 0
    p    = self%unif_01()

    do while (p > expl)
       rr = rr + 1
       p = p * self%unif_01()
    end do
  end function poisson

  !> Sample point on a circle with given radius
  function circle(self, radius) result(xy)
    class(rng_t), intent(inout) :: self
    real(dp), intent(in)        :: radius
    real(dp)                    :: rands(2), xy(2)
    real(dp)                    :: sum_sq

    ! Method for uniform sampling on circle
    do
       rands(1) = 2 * self%unif_01() - 1
       rands(2) = 2 * self%unif_01() - 1
       sum_sq   = sum(rands**2)
       if (sum_sq <= 1) exit
    end do

    xy(1) = (rands(1)**2 - rands(2)**2) / sum_sq
    xy(2) = 2 * rands(1) * rands(2) / sum_sq
    xy    = xy * radius
  end function circle

  !> Sample point on a sphere with given radius
  function sphere(self, radius) result(xyz)
    class(rng_t), intent(inout) :: self
    real(dp), intent(in)        :: radius
    real(dp)                    :: rands(2), xyz(3)
    real(dp)                    :: sum_sq, tmp_sqrt

    ! Marsaglia method for uniform sampling on sphere
    do
       rands(1) = 2 * self%unif_01() - 1
       rands(2) = 2 * self%unif_01() - 1
       sum_sq   = sum(rands**2)
       if (sum_sq <= 1) exit
    end do

    tmp_sqrt = sqrt(1 - sum_sq)
    xyz(1:2) = 2 * rands(1:2) * tmp_sqrt
    xyz(3)   = 1 - 2 * sum_sq
    xyz      = xyz * radius
  end function sphere

    !> Interal routine: get the next value (returned as 64 bit signed integer)
  function next(self) result(res)
    class(rng_t), intent(inout) :: self
    integer(i8)                 :: res
    integer(i8)                 :: t(2)

    t         = self%s
    res       = t(1) + t(2)
    t(2)      = ieor(t(1), t(2))
    self%s(1) = ieor(ieor(rotl(t(1), 55), t(2)), shiftl(t(2), 14))
    self%s(2) = rotl(t(2), 36)
  end function next

  !> Helper function for next()
  pure function rotl(x, k) result(res)
    integer(i8), intent(in) :: x
    integer, intent(in)     :: k
    integer(i8)             :: res

    res = ior(shiftl(x, k), shiftr(x, 64 - k))
  end function rotl

end module m_random