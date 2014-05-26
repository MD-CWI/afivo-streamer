!> Module for the random number generation
module m_random
   implicit none
   private

   integer, parameter :: dp = kind(0.0d0)
   integer, parameter :: i4 = selected_int_kind(9)

   type, public :: RNG_state_t
      integer(i4), private :: qq(4)
   end type RNG_state_t

   ! Procedures
   public :: RNG_int4
   public :: RNG_U01
   public :: RNG_set_state
   public :: RNG_normal
   public :: RNG_poisson

contains

   subroutine RNG_set_state(state, initializer)
      type(RNG_state_t), intent(out) :: state
      integer(i4), intent(in), optional :: initializer(4)

      if (present(initializer)) then
         state%qq = initializer
      else
         state%qq = (/521288629, 362436069, 16163801, 1131199299/)
      end if
   end subroutine RNG_set_state

   ! Method mzran (Marsaglia), see http://jblevins.org/log/openmp
   function RNG_int4(self) result(r_int)
      type(RNG_state_t), intent(inout) :: self
      integer(i4)                    :: r_int

      r_int = self%qq(1) - self%qq(3)
      if (r_int < 0) r_int = r_int + 2147483579

      self%qq(1:2) = self%qq(2:3)
      self%qq(3)   = r_int
      self%qq(4)   = 69069 * self%qq(4) + 1013904243
      r_int        = r_int + self%qq(4)
   end function RNG_int4

   ! Random [0,1) number
   function RNG_U01(self) result(r_U01)
      type(RNG_state_t), intent(inout) :: self
      real(dp)                       :: r_U01
      real(dp), parameter            :: conv_fac = 0.5_dp / 2.0_dp**31

      r_U01 = 0.5_dp + conv_fac * RNG_int4(self)
   end function RNG_U01

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
   integer(i4) function RNG_poisson(labda)
      real(dp), intent(IN) :: labda
      integer(i4) :: k
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
