#include "cpp_macros.h"
!> Module containing slope limiters
module m_af_limiters
  use m_af_types

  implicit none
  private

  !> Number of limiters
  integer, parameter :: af_num_limiters = 7

  ! Constants to identify limiters
  integer, parameter, public :: af_limiter_none_t = 1
  !> van Leer limiter
  integer, parameter, public :: af_limiter_vanleer_t = 2
  !> Koren limiter
  integer, parameter, public :: af_limiter_koren_t = 3
  !> Minmod limiter
  integer, parameter, public :: af_limiter_minmod_t = 4
  !> MC limiter
  integer, parameter, public :: af_limiter_mc_t = 5
  !> Generalized minmod limiter with theta = 4/3
  integer, parameter, public :: af_limiter_gminmod43_t = 6
  !> All slopes are zero
  integer, parameter, public :: af_limiter_zero_t = 7

  !> Whether limiters are symmetric
  logical, parameter, public :: af_limiter_symmetric(af_num_limiters) = &
       [.true., .true., .false., .true., .true., .true., .true.]

  public :: af_limiter_apply
  public :: af_limiter_koren
  public :: af_limiter_vanleer
  public :: af_limiter_minmod
  public :: af_limiter_gminmod43
  public :: af_limiter_mc

contains

  !> Apply one of the limiters
  elemental function af_limiter_apply(a, b, limiter) result(slope)
    real(dp), intent(in) :: a       !< Slopes from one side
    real(dp), intent(in) :: b       !< Slopes from other side
    integer, intent(in)  :: limiter !< Which limiter to use
    real(dp)             :: slope   !< Limited slope

    select case (limiter)
    case (af_limiter_none_t)
       slope = 0.5_dp * (a + b)
    case (af_limiter_vanleer_t)
       slope = af_limiter_vanleer(a, b)
    case (af_limiter_koren_t)
       slope = af_limiter_koren(a, b)
    case (af_limiter_minmod_t)
       slope = af_limiter_minmod(a, b)
    case (af_limiter_mc_t)
       slope = af_limiter_mc(a, b)
    case (af_limiter_gminmod43_t)
       slope = af_limiter_gminmod43(a, b)
    case (af_limiter_zero_t)
       slope = 0.0_dp
    case default
       ! Cannot (yet) stop in elemental function, so use default limiter
       slope = af_limiter_vanleer(a, b)
    end select
  end function af_limiter_apply

  !> Modified implementation of Koren limiter, to avoid division and the min/max
  !> functions, which can be problematic / expensive. In most literature, you
  !> have r = a / b (ratio of gradients). Then the limiter phi(r) is multiplied
  !> with b. With this implementation, you get phi(r) * b
  elemental function af_limiter_koren(a, b) result(bphi)
    real(dp), intent(in) :: a  !< Density gradient (numerator)
    real(dp), intent(in) :: b  !< Density gradient (denominator)
    real(dp), parameter  :: third = 1/3.0_dp
    real(dp)             :: bphi, aa, ab

    aa = a * a
    ab = a * b

    if (ab <= 0) then
       ! a and b have different sign or one of them is zero, so r is either 0,
       ! inf or negative (special case a == b == 0 is ignored)
       bphi = 0
    else if (aa <= 0.25_dp * ab) then
       ! 0 < a/b <= 1/4, limiter has value 2*a/b
       bphi = 2*a
    else if (aa <= 2.5_dp * ab) then
       ! 1/4 < a/b <= 2.5, limiter has value (1+2*a/b)/3
       bphi = third * (b + 2*a)
    else
       ! (1+2*a/b)/6 >= 1, limiter has value 2
       bphi = 2*b
    end if
  end function af_limiter_koren

  elemental function af_limiter_vanleer(a, b) result(phi)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp)             :: ab, phi

    ab = a * b
    if (ab > 0) then
       ! Jannis: I think roundoff errors can lead to values slightly outside the
       ! desired range
       phi = 2 * ab / (a + b)
    else
       phi = 0
    end if
  end function af_limiter_vanleer

  elemental function af_limiter_minmod(a, b) result(phi)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp)             :: phi
    phi = af_limiter_gminmod(a, b, 1.0_dp)
  end function af_limiter_minmod

  ! Generalized minmod limiter with parameter theta = 4/3
  elemental function af_limiter_gminmod43(a, b) result(phi)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp)             :: phi
    real(dp), parameter  :: theta = 4/3.0_dp
    phi = af_limiter_gminmod(a, b, theta)
  end function af_limiter_gminmod43

  ! MC (monotonized central) limiter
  elemental function af_limiter_mc(a, b) result(phi)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp)             :: phi
    phi = af_limiter_gminmod(a, b, 2.0_dp)
  end function af_limiter_mc

  !> Generalized minmod limiter. The parameter theta controls how dissipative
  !> the limiter is, with 1 corresponding to the minmod limiter and 2 to the MC
  !> limiter.
  elemental function af_limiter_gminmod(a, b, theta) result(phi)
    real(dp), intent(in) :: a, b, theta
    real(dp)             :: phi

    if (a * b > 0) then
       phi = sign(minval(abs([theta * a, theta * b, &
            0.5_dp * (a + b)])), a)
    else
       phi = 0.0_dp
    end if
  end function af_limiter_gminmod

end module m_af_limiters
