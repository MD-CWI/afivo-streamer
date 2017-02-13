#include "cpp_macros_$Dd.h"
!> This module contains routines related to point-based interpolation, which can
!> be useful when you want to have output at given points. The interpolation for
!> meshes is called prolongation, see m_aX_prolong.
module m_a$D_interp
  use m_a$D_types

  implicit none
  private

  public :: a$D_interp1

contains

  !> Using linear interpolation to get a value at r
  function a$D_interp1(tree, r, ivs, n_var) result(vals)
    use m_a$D_utils, only: a$D_get_loc, a$D_r_loc
    type(a$D_t), intent(in) :: tree !< Parent box
    real(dp), intent(in)    :: r($D) !< Where to interpolate
    integer, intent(in)     :: n_var     !< Number of variables
    integer, intent(in)     :: ivs(n_var)   !< Variables to interpolate
    real(dp)                :: vals(n_var)
    integer                 :: i, iv, id, ix($D)
    real(dp)                :: r_loc($D), dvec($D), ovec($D), w(DTIMES(2))
    type(a$D_loc_t)         :: loc  !< Location of cell

    loc   = a$D_get_loc(tree, r)
    id    = loc%id

    if (id <= af_no_box) then
       print *, "a$D_interp1: point outside domain", r
       stop
    end if

    r_loc = a$D_r_loc(tree, loc)
    ix    = loc%ix
    dvec  = r - r_loc

    ! Compute ix such that r lies between ix and ix + 1
    where (r < r_loc)
       ix   = ix - 1
       dvec = dvec + tree%boxes(id)%dr
    end where

    ! Normalize dvec to a value [0, 1]
    dvec = dvec / tree%boxes(id)%dr
    ovec = 1 - dvec

    ! Compute weights of linear interpolation
#if $D == 2
    w(1, 1) = ovec(1) * ovec(2)
    w(2, 1) = dvec(1) * ovec(2)
    w(1, 2) = ovec(1) * dvec(2)
    w(2, 2) = dvec(1) * dvec(2)
#elif $D == 3
    w(1, 1, 1) = ovec(1) * ovec(2) * ovec(3)
    w(2, 1, 1) = dvec(1) * ovec(2) * ovec(3)
    w(1, 2, 1) = ovec(1) * dvec(2) * ovec(3)
    w(2, 2, 1) = dvec(1) * dvec(2) * ovec(3)
    w(1, 1, 2) = ovec(1) * ovec(2) * dvec(3)
    w(2, 1, 2) = dvec(1) * ovec(2) * dvec(3)
    w(1, 2, 2) = ovec(1) * dvec(2) * dvec(3)
    w(2, 2, 2) = dvec(1) * dvec(2) * dvec(3)
#endif

    do i = 1, size(ivs)
       iv = ivs(i)
#if $D == 2
       vals(i) = sum(w * tree%boxes(id)%cc(ix(1):ix(1)+1, &
            ix(2):ix(2)+1, iv))
#elif $D == 3
       vals(i) = sum(w * tree%boxes(id)%cc(ix(1):ix(1)+1, &
            ix(2):ix(2)+1, ix(3):ix(3)+1, iv))
#endif
     end do
  end function a$D_interp1

end module m_a$D_interp
