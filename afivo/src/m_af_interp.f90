#include "cpp_macros.h"
!> This module contains routines related to interpolation, which can interpolate
!> 'to' the grid and 'from' the grid (useful for e.g. particle simulations). The
!> interpolation for meshes is called prolongation, see m_aX_prolong.
module m_af_interp
  use m_af_types

  implicit none
  private

  public :: af_interp0
  public :: af_interp1
  public :: af_interp0_to_grid
  public :: af_interp1_fc

contains

  !> Using zeroth order interpolation to get a value at r in cell
  function af_interp0(tree, r, ivs, success, id_guess) result(vals)
    use m_af_utils, only: af_get_loc
    type(af_t), intent(in)        :: tree     !< Parent box
    real(dp), intent(in)          :: r(NDIM) !< Where to interpolate
    integer, intent(in)           :: ivs(:)   !< Variables to interpolate
    logical, intent(out)          :: success  !< Whether the interpolation worked
    integer, intent(inout), optional :: id_guess !< Guess for box id (will be updated)
    real(dp)                      :: vals(size(ivs))
    type(af_loc_t)                :: loc

    loc = af_get_loc(tree, r, guess=id_guess)
    ! Update guess
    if (present(id_guess)) id_guess = loc%id

    if (loc%id == -1) then
       success = .false.
       vals = 0.0_dp
    else
       success = .true.
       vals = tree%boxes(loc%id)%cc(DINDEX(loc%ix), ivs)
    end if
  end function af_interp0

  !> Using linear interpolation to get a value at r
  function af_interp1(tree, r, ivs, success, id_guess) result(vals)
    use m_af_utils, only: af_get_id_at
    type(af_t), intent(in)        :: tree     !< Parent box
    real(dp), intent(in)          :: r(NDIM)  !< Where to interpolate
    integer, intent(in)           :: ivs(:)   !< Variables to interpolate
    logical, intent(out)          :: success  !< Whether the interpolation worked
    integer, intent(inout), optional :: id_guess !< Guess for box id (will be updated)
    real(dp)                      :: vals(size(ivs))
    integer                       :: i, iv, id, ix(NDIM)
    real(dp)                      :: r_loc(NDIM), dvec(NDIM), ovec(NDIM), w(DTIMES(2))

    id = af_get_id_at(tree, r, guess=id_guess)
    ! Update guess
    if (present(id_guess)) id_guess = id

    if (id <= af_no_box) then
       success = .false.
       vals = 0.0_dp
    else
       success = .true.
       ! Compute ix such that r lies between cell centers at ix and ix + 1
       ix = nint((r - tree%boxes(id)%r_min) / tree%boxes(id)%dr)
       r_loc = af_r_cc(tree%boxes(id), ix)
       dvec  = r - r_loc

       ! Normalize dvec to a value [0, 1]
       dvec = dvec / tree%boxes(id)%dr
       ovec = 1 - dvec

       ! Compute weights of linear interpolation
#if NDIM == 1
       w(1) = ovec(1)
       w(2) = dvec(1)
#elif NDIM == 2
       w(1, 1) = ovec(1) * ovec(2)
       w(2, 1) = dvec(1) * ovec(2)
       w(1, 2) = ovec(1) * dvec(2)
       w(2, 2) = dvec(1) * dvec(2)
#elif NDIM == 3
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
#if NDIM == 1
          vals(i) = sum(w * tree%boxes(id)%cc(ix(1):ix(1)+1, iv))
#elif NDIM == 2
          vals(i) = sum(w * tree%boxes(id)%cc(ix(1):ix(1)+1, &
               ix(2):ix(2)+1, iv))
#elif NDIM == 3
          vals(i) = sum(w * tree%boxes(id)%cc(ix(1):ix(1)+1, &
               ix(2):ix(2)+1, ix(3):ix(3)+1, iv))
#endif
       end do
    end if
  end function af_interp1

  !> Add 'amount' to the grid cell nearest to r
  subroutine af_interp0_to_grid(tree, r, iv, amount, to_density)
    use m_af_utils, only: af_get_loc
    type(af_t), intent(inout) :: tree
    integer, intent(in)        :: iv         !< Index of variable
    real(dp), intent(in)       :: r(NDIM)     !< Location
    real(dp), intent(in)       :: amount     !< How much to add
    logical, intent(in)        :: to_density !< If true, divide by cell volume
    real(dp)                   :: actual_amount
    type(af_loc_t)            :: loc
    integer                    :: id, ix(NDIM)

    loc = af_get_loc(tree, r)

    if (loc%id == -1) then
       print *, "af_interp0_to_grid error, no box at ", r
       stop
    end if

    id = loc%id
    ix = loc%ix

    !> @todo Support cylindrical coordinates
    if (to_density) then
       actual_amount = amount / product(tree%boxes(id)%dr)
    else
       actual_amount = amount
    end if

    tree%boxes(id)%cc(DINDEX(ix), iv) = &
         tree%boxes(id)%cc(DINDEX(ix), iv) + &
         actual_amount
  end subroutine af_interp0_to_grid

  !> Linearly interpolate face-centered variable to a position r
  function af_interp1_fc(tree, r, ifc, success, id_guess) result(vals)
    use m_af_utils, only: af_get_id_at
    type(af_t), intent(in)           :: tree     !< Parent box
    real(dp), intent(in)             :: r(NDIM)  !< Where to interpolate
    integer, intent(in)              :: ifc      !< Face-centered variable
    logical, intent(out)             :: success  !< Whether the interpolation worked
    integer, intent(inout), optional :: id_guess !< Guess for box id (will be updated)
    real(dp)                         :: vals(NDIM), inv_dr(NDIM)
    integer                          :: id, ix(NDIM)
    real(dp)                         :: ix_frac(NDIM), r_rel(NDIM)

    id = af_get_id_at(tree, r, guess=id_guess)

    ! Update guess
    if (present(id_guess)) id_guess = id

    if (id <= af_no_box) then
       success = .false.
       vals = 0.0_dp
    else
       success = .true.
       inv_dr  = 1/tree%boxes(id)%dr

       r_rel = r - tree%boxes(id)%r_min
       ix_frac = r_rel * inv_dr + 1
       ix = floor(ix_frac)
       where (ix < 1) ix = 1
       where (ix > tree%n_cell) ix = tree%n_cell
       ix_frac = ix_frac - ix

#if NDIM == 2
       vals(1) = (1 - ix_frac(1)) * tree%boxes(id)%fc(ix(1), ix(2), 1, ifc) + &
            ix_frac(1) * tree%boxes(id)%fc(ix(1)+1, ix(2), 1, ifc)
       vals(2) = (1 - ix_frac(2)) * tree%boxes(id)%fc(ix(1), ix(2), 2, ifc) + &
            ix_frac(2) * tree%boxes(id)%fc(ix(1), ix(2)+1, 2, ifc)
#elif NDIM == 3
       vals(1) = (1 - ix_frac(1)) * tree%boxes(id)%fc(ix(1), ix(2), ix(3), 1, ifc) + &
            ix_frac(1) * tree%boxes(id)%fc(ix(1)+1, ix(2), ix(3), 1, ifc)
       vals(2) = (1 - ix_frac(2)) * tree%boxes(id)%fc(ix(1), ix(2), ix(3), 2, ifc) + &
            ix_frac(2) * tree%boxes(id)%fc(ix(1), ix(2)+1, ix(3), 2, ifc)
       vals(3) = (1 - ix_frac(3)) * tree%boxes(id)%fc(ix(1), ix(2), ix(3), 3, ifc) + &
            ix_frac(3) * tree%boxes(id)%fc(ix(1), ix(2), ix(3)+1, 3, ifc)
#endif
    end if
  end function af_interp1_fc

end module m_af_interp
