!> Module for the random number generation
module m_random
   implicit none
   private

   integer, parameter :: dp = kind(0.0d0)

   public :: RNG_normal
   public :: RNG_poisson
   public :: RNG_set_seed

contains

   !> Initialize the RNG seed
   subroutine RNG_set_seed(arg_seed)
      integer, intent(IN)  :: arg_seed
      integer              :: n, seed_size
      integer, allocatable :: used_seed(:)

      call random_seed(size = seed_size) ! Get required size of seed
      allocate(used_seed(seed_size))

      do n = 1, seed_size
         used_seed(n) = arg_seed * n ! Fill with some arbritrary values depending on arg_seed
      end do

      call random_seed(put = used_seed)
      deallocate(used_seed)
   end subroutine RNG_set_seed

   !> Return normal random variate with mean 0 and variance 1
   ! Source: http://en.wikipedia.org/wiki/Marsaglia_polar_method
   real(dp) function RNG_normal()
      logical, save :: have_spare = .false.
      real(dp), save :: spare
      real(dp) :: rands(2), sum_sq, mul

      if (have_spare) then
         have_spare = .false.
         RNG_normal = spare
      else
         do
            call random_number(rands)
            rands = 2 * rands - 1
            sum_sq = sum(rands**2)
            if (sum_sq < 1.0_dp .and. sum_sq /= 0.0_dp) exit
         end do

         mul        = sqrt(-2 * log(sum_sq) / sum_sq)
         RNG_normal = rands(1) * mul
         spare      = rands(2) * mul
         have_spare = .true.
      end if
   end function RNG_normal

   !> Return Poisson random variate with rate labda. Works well for labda < 30 or so.
   !! For labda >> 1 it can produce wrong results due to roundoff error.
   integer function RNG_poisson(labda)
      real(dp), intent(IN) :: labda
      integer :: k
      real(dp) :: expL, p

      expL = exp(-labda)
      k = 0
      call random_number(p)

      do while (p > expL)
         k = k + 1
         call random_number(p)
      end do
      RNG_poisson = k

   end function RNG_poisson

end module m_random
