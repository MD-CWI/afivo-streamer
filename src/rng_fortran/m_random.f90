!> Module for pseudo random number generation. The internal pseudo random
!> generator is the xoroshiro128plus method.
module m_random
  use iso_fortran_env, only: int32, int64, real64

  implicit none
  private

  ! A 64 bit floating point type
  integer, parameter :: dp = real64

  !> Random number generator type, which contains the state
  type rng_t
     !> The rng state (always use your own seed)
     integer(int64), private    :: s(2) = [123456789_int64, 987654321_int64]
     integer(int64), private    :: separator(32) ! Separate cache lines (parallel use)
     real(dp), private          :: stored_normal
     logical, private           :: have_stored_normal = .false.
     logical, private           :: initialized = .false.
   contains
     procedure, non_overridable :: set_seed    ! Seed the generator
     procedure, non_overridable :: set_random_seed ! Use a random seed
     procedure, non_overridable :: jump        ! Jump function (see below)
     procedure, non_overridable :: int_32       ! 4-byte random integer
     procedure, non_overridable :: int_64       ! 8-byte random integer
     procedure, non_overridable :: unif_01     ! Uniform (0,1] real
     procedure, non_overridable :: normal      ! One normal(0,1) sample
     procedure, non_overridable :: two_normals ! Two normal(0,1) samples
     procedure, non_overridable :: poisson     ! Sample from Poisson-dist.
     procedure, non_overridable :: poisson_knuth ! Sample from Poisson-dist.
     procedure, non_overridable :: poisson_reject ! Sample from Poisson-dist.
     procedure, non_overridable :: exponential_standard ! Sample from standard exponential dist.
     procedure, non_overridable :: exponential ! Sample from exponential dist.
     procedure, non_overridable :: geometric   ! Sample from geometric distribution
     procedure, non_overridable :: circle      ! Sample on a circle
     procedure, non_overridable :: sphere      ! Sample on a sphere
     procedure, non_overridable :: next        ! Internal method
     procedure, non_overridable, nopass :: log1p ! Internal method
  end type rng_t

  !> Parallel random number generator type
  type prng_t
     type(rng_t), allocatable :: rngs(:)
   contains
     procedure, non_overridable :: init_parallel
     procedure, non_overridable :: update_seed
  end type prng_t

  public :: rng_t
  public :: prng_t

contains

  !> Initialize a collection of rng's for parallel use
  subroutine init_parallel(self, n_proc, rng)
    class(prng_t), intent(inout) :: self
    type(rng_t), intent(inout)   :: rng
    integer, intent(in)          :: n_proc
    integer                      :: n

    if (n_proc < 1) error stop "init_parallel: n_proc < 1"

    allocate(self%rngs(n_proc))
    self%rngs(1) = rng
    call self%rngs(1)%jump()

    do n = 2, n_proc
       self%rngs(n) = self%rngs(n-1)
       call self%rngs(n)%jump()
    end do
  end subroutine init_parallel

  !> Parallel RNG instances are often used temporarily. This routine can
  !> afterwards be used to update the seed of the user's sequential RNG.
  subroutine update_seed(self, rng)
    class(prng_t), intent(inout) :: self
    type(rng_t), intent(inout)   :: rng
    integer                      :: n

    do n = 1, size(self%rngs)
       ! Perform exclusive-or with each parallel rng
       rng%s(1) = ieor(rng%s(1), self%rngs(n)%s(1))
       rng%s(2) = ieor(rng%s(2), self%rngs(n)%s(2))
    end do
  end subroutine update_seed

  !> Set a seed for the rng
  subroutine set_seed(self, the_seed)
    class(rng_t), intent(inout) :: self
    integer(int64), intent(in)  :: the_seed(2)

    self%s = the_seed

    ! Simulate calls to next() to improve randomness of first number
    call self%jump()
  end subroutine set_seed

  subroutine set_random_seed(self)
    class(rng_t), intent(inout) :: self
    integer                     :: i
    real(dp)                    :: rr
    integer(int64)              :: time

    ! Get a random seed from the system (this does not always work)
    call random_seed()

    ! Get some count of the time
    call system_clock(time)

    do i = 1, 2
       call random_number(rr)
       self%s(i) = ieor(transfer(rr, 1_int64), transfer(time, 1_int64))
    end do

    ! Simulate calls to next() to improve randomness of first number
    call self%jump()
  end subroutine set_random_seed

  ! This is the jump function for the generator. It is equivalent
  ! to 2^64 calls to next(); it can be used to generate 2^64
  ! non-overlapping subsequences for parallel computations.
  subroutine jump(self)
    class(rng_t), intent(inout) :: self
    integer                     :: i, b
    integer(int64)              :: t(2), dummy

    ! The signed equivalent of the unsigned constants
    integer(int64), parameter :: jmp_c(2) = &
         (/-4707382666127344949_int64, -2852180941702784734_int64/)

    t = 0
    do i = 1, 2
       do b = 0, 63
          if (iand(jmp_c(i), shiftl(1_int64, b)) /= 0) then
             t = ieor(t, self%s)
          end if
          dummy = self%next()
       end do
    end do

    self%s = t
  end subroutine jump

  !> Return 4-byte integer
  integer(int32) function int_32(self)
    class(rng_t), intent(inout) :: self
    int_32 = int(self%next(), int32)
  end function int_32

  !> Return 8-byte integer
  integer(int64) function int_64(self)
    class(rng_t), intent(inout) :: self
    int_64 = self%next()
  end function int_64

  !> Get a uniform [0,1) random real (double precision)
  real(dp) function unif_01(self)
    class(rng_t), intent(inout) :: self
    integer(int64)              :: x
    real(dp)                    :: tmp

    x   = self%next()
    x   = ior(shiftl(1023_int64, 52), shiftr(x, 12))
    unif_01 = transfer(x, tmp) - 1.0_dp
  end function unif_01

  !> Return normal random variate with mean 0 and variance 1
  real(dp) function normal(self)
    class(rng_t), intent(inout) :: self
    real(dp)                    :: two_normals(2)

    if (self%have_stored_normal) then
       normal = self%stored_normal
       self%have_stored_normal = .false.
    else
       two_normals = self%two_normals()
       normal      = two_normals(1)
       self%stored_normal = two_normals(2)
       self%have_stored_normal = .true.
    end if
  end function normal

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

  !> Compute log(1+x) with good accuracy, see "What Every Computer Scientist
  !> Should Know About Floating-Point Arithmetic"
  real(dp) function log1p(x)
    real(dp), intent(in) :: x

    if (1.0_dp + abs(x) > 1.0_dp) then
       log1p = log(1.0_dp + x) * x / ((1.0_dp + x) - 1.0_dp)
    else
       log1p = x
    endif
  end function log1p

  !> Return exponential random variate with a rate of one
  real(dp) function exponential_standard(self)
    class(rng_t), intent(inout) :: self
    real(dp)                    :: unif_01

    ! It is assumes 1 - unif_01 is in (0, 1], so we avoid log(0.) below
    unif_01 = self%unif_01()

    if (unif_01 < 0.5_dp) then
       exponential_standard = -log1p(-unif_01)
    else
       exponential_standard = -log(1 - unif_01)
    end if
  end function exponential_standard

  !> Return exponential random variate with rate lambda
  real(dp) function exponential(self, lambda)
    class(rng_t), intent(inout) :: self
    real(dp), intent(in)        :: lambda
    exponential = self%exponential_standard()/lambda
  end function exponential

  !> Sample from geometric distribution with Pr(X = k) = (1 - p)^(k-1) * p
  integer(int64) function geometric(self, p)
    class(rng_t), intent(inout) :: self
    real(dp), intent(in)        :: p
    real(dp)                    :: tmp
    real(dp), parameter         :: threshold = real(huge(1_int64) - 1, dp)

    ! Perform inversion sampling X = ceiling(log(U)/log(1-p))
    tmp = -self%exponential_standard() / log1p(-p)

    ! Avoid overflow
    if (tmp < threshold) then
       geometric = ceiling(tmp, int64)
    else
       geometric = huge(1_int64)
    end if

  end function geometric

  !> Return Poisson random variate with rate lambda. Works well for lambda < 30
  !> or so. For lambda >> 1 it can produce wrong results due to roundoff error.
  function poisson_knuth(self, lambda) result(rr)
    class(rng_t), intent(inout) :: self
    real(dp), intent(in)        :: lambda
    integer(int32)                 :: rr
    real(dp)                    :: expl, p

    expl = exp(-lambda)
    rr   = 0
    p    = self%unif_01()

    do while (p > expl)
       rr = rr + 1
       p = p * self%unif_01()
    end do
  end function poisson_knuth

  !> The transformed rejection method for generating Poisson random variables
  !>
  !> Translated from Numpy C code at:
  !> https://github.com/numpy/numpy/blob/main/numpy/random/src/distributions/distributions.c
  !>
  !> W. Hoermann
  !> Insurance: Mathematics and Economics 12, 39-45 (1993)
  function poisson_reject(self, lambda) result(k)
    class(rng_t), intent(inout) :: self
    real(dp), intent(in)        :: lambda
    integer(int32)              :: k
    real(dp)                    :: U, V, sqrt_lambda, log_lambda
    real(dp)                    :: a, b, invalpha, vr, us

    sqrt_lambda = sqrt(lambda)
    log_lambda = log(lambda)

    b = 0.931_dp + 2.53_dp * sqrt_lambda
    a = -0.059_dp + 0.02483_dp * b
    invalpha = 1.1239_dp + 1.1328_dp / (b - 3.4_dp)
    vr = 0.9277_dp - 3.6224_dp / (b - 2)

    do
       U = self%unif_01() - 0.5_dp
       V = 1.0_dp - self%unif_01() ! Avoid 0
       us = 0.5_dp - abs(U);

       k = floor((2 * a / us + b) * U + lambda + 0.43_dp);

       if (us >= 0.07_dp .and. V <= vr) return
       if (k < 0 .or. us < 0.013_dp .and. V > us) cycle

       if ((log(V) + log(invalpha) - log(a / (us * us) + b)) <= &
            (-lambda + k * log_lambda - log_gamma(k + 1.0_dp))) return
    end do
  end function poisson_reject

  !> Return Poisson random variate with rate lambda
  function poisson(self, lambda) result(rr)
    class(rng_t), intent(inout) :: self
    real(dp), intent(in)        :: lambda
    integer(int32)              :: rr

    if (lambda < 10) then
       ! Algorithm for small value of lambda
       rr = self%poisson_knuth(lambda)
    else
       ! Rejection sampling
       rr = self%poisson_reject(lambda)
    end if
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
    integer(int64)              :: res
    integer(int64)              :: t(2)

    t         = self%s
    res       = t(1) + t(2)
    t(2)      = ieor(t(1), t(2))
    self%s(1) = ieor(ieor(rotl(t(1), 55), t(2)), shiftl(t(2), 14))
    self%s(2) = rotl(t(2), 36)
  end function next

  !> Helper function for next()
  pure function rotl(x, k) result(res)
    integer(int64), intent(in) :: x
    integer, intent(in)        :: k
    integer(int64)             :: res

    res = ior(shiftl(x, k), shiftr(x, 64 - k))
  end function rotl

end module m_random
