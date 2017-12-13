#include "cpp_macros_$Dd.h"
!> This module contains routines related to interpolation, which can interpolate
!> 'to' the grid and 'from' the grid (useful for e.g. particle simulations). The
!> interpolation for meshes is called prolongation, see m_aX_prolong.
module m_a$D_interp
  use m_a$D_types

  implicit none
  private

  public :: a$D_interp0
  public :: a$D_interp1
  public :: a$D_interp0_to_grid
  public :: a$D_interp1_to_grid

contains

  !> Using zeroth order interpolation to get a value at r in cell
  function a$D_interp0(tree, rr, ivs, n_var) result(vals)
    use m_a$D_utils, only: a$D_get_loc
    type(a$D_t), intent(in) :: tree !< Parent box
    real(dp), intent(in)    :: rr($D) !< Where to interpolate
    integer, intent(in)     :: n_var     !< Number of variables
    integer, intent(in)     :: ivs(n_var)   !< Variables to interpolate
    real(dp)                :: vals(n_var)
    type(a$D_loc_t)         :: loc

    loc = a$D_get_loc(tree, rr)

    if (loc%id == -1) then
       print *, "a$D_interp0 error, no box at ", rr
       stop
    end if

#if $D == 2
    vals = tree%boxes(loc%id)%cc(loc%ix(1), loc%ix(2), ivs)
#elif $D == 3
    vals = tree%boxes(loc%id)%cc(loc%ix(1), loc%ix(2), loc%ix(3), ivs)
#endif
  end function a$D_interp0

  !> Using linear interpolation to get a value at r
  function a$D_interp1(tree, r, ivs, n_var) result(vals)
    use m_a$D_utils, only: a$D_get_id_at
    type(a$D_t), intent(in) :: tree !< Parent box
    real(dp), intent(in)    :: r($D) !< Where to interpolate
    integer, intent(in)     :: n_var     !< Number of variables
    integer, intent(in)     :: ivs(n_var)   !< Variables to interpolate
    real(dp)                :: vals(n_var)
    integer                 :: i, iv, id, ix($D)
    real(dp)                :: r_loc($D), dvec($D), ovec($D), w(DTIMES(2))

    id = a$D_get_id_at(tree, r)

    if (id <= af_no_box) then
       print *, "a$D_interp1: point outside domain", r
       stop
    end if

    ! Compute ix such that r lies between cell centers at ix and ix + 1
    ix = nint((r - tree%boxes(id)%r_min) / tree%boxes(id)%dr)
    r_loc = a$D_r_cc(tree%boxes(id), ix)
    dvec  = r - r_loc

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

  !> Add 'amount' to the grid cell nearest to rr
  subroutine a$D_interp0_to_grid(tree, rr, iv, amount, to_density)
    use m_a$D_utils, only: a$D_get_loc
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)        :: iv         !< Index of variable
    real(dp), intent(in)       :: rr($D)     !< Location
    real(dp), intent(in)       :: amount     !< How much to add
    logical, intent(in)        :: to_density !< If true, divide by cell volume
    real(dp)                   :: actual_amount
    type(a$D_loc_t)            :: loc
    integer                    :: id, ix($D)

    loc = a$D_get_loc(tree, rr)

    if (loc%id == -1) then
       print *, "a$D_interp0_to_grid error, no box at ", rr
       stop
    end if

    id = loc%id
    ix = loc%ix

    !> @todo Support cylindrical coordinates
    if (to_density) then
       actual_amount = amount / tree%boxes(id)%dr**$D
    else
       actual_amount = amount
    end if

#if $D == 2
    tree%boxes(id)%cc(ix(1), ix(2), iv) = &
         tree%boxes(id)%cc(ix(1), ix(2), iv) + &
         actual_amount
#elif $D == 3
    tree%boxes(id)%cc(ix(1), ix(2), ix(3), iv) = &
         tree%boxes(id)%cc(ix(1), ix(2), ix(3), iv) + &
         actual_amount
#endif
  end subroutine a$D_interp0_to_grid

  !> Add 'amount' to the cell centers around rr, using linear interpolation
  subroutine a$D_interp1_to_grid(tree, rr, iv, amount, to_density)
    use m_a$D_utils, only: a$D_get_id_at
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)        :: iv         !< Index of variable
    real(dp), intent(in)       :: rr($D)     !< Location
    real(dp), intent(in)       :: amount     !< How much to add
    logical, intent(in)        :: to_density !< If true, divide by cell volume
    real(dp)                   :: actual_amount, inv_dr, tmp($D)
    real(dp)                   :: wu($D), wl($D), w(DTIMES(2))
    integer                    :: id, ix($D), nc
    logical                    :: rb_low($D), rb_high($D)

    id = a$D_get_id_at(tree, rr)

    if (id == -1) then
       print *, "a$D_interp1_to_grid error, no box at ", rr
       stop
    end if

    nc     = tree%boxes(id)%n_cell
    inv_dr = 1.0_dp/tree%boxes(id)%dr
    tmp    = (rr - tree%boxes(id)%r_min) * inv_dr + 0.5_dp
    ix     = floor(tmp)
    wu     = tmp - ix
    wl     = 1 - wu

#if $D == 2
    w(:, 1) = [wl(1), wu(1)] * wl(2)
    w(:, 2) = [wl(1), wu(1)] * wu(2)
#elif $D == 3
#endif

    !> @todo Support cylindrical coordinates
    if (to_density) then
       actual_amount = amount / tree%boxes(id)%dr**$D
    else
       actual_amount = amount
    end if

    if (any(ix == 0) .or. any(ix == nc)) then
       ix = ceiling((rr - tree%boxes(id)%r_min) * inv_dr)

       ! Apply special method near refinement boundaries
#if $D == 2
       tree%boxes(id)%cc(ix(1), ix(2), iv) = &
            tree%boxes(id)%cc(ix(1), ix(2), iv) + &
            actual_amount
#elif $D == 3
       tree%boxes(id)%cc(ix(1), ix(2), ix(3), iv) = &
            tree%boxes(id)%cc(ix(1), ix(2), ix(3), iv) + &
            actual_amount
#endif
    else
       ! Linear interpolation
#if $D == 2
       tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, iv) = &
            tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, iv) + &
            w * actual_amount
#elif $D == 3
       ! tree%boxes(id)%cc(ix(1), ix(2), ix(3), iv) = &
       !      tree%boxes(id)%cc(ix(1), ix(2), ix(3), iv) + &
       !      actual_amount
#endif
    end if
  end subroutine a$D_interp1_to_grid

end module m_a$D_interp
