! This module contains routines related to interpolation
!
! Author: Jannis Teunissen
! License: GPLv3

module m_a$D_interp
  use m_a$D_t

  implicit none
  private

  public :: a$D_interp1
  public :: a$D_interp2

contains

  !> Using second order interpolation to get a value at r. The result is not
  !> continuous at cell boundaries.
  function a$D_interp1(tree, r, ivs, n_var) result(vals)
    use m_a$D_utils, only: a$D_get_loc, a$D_r_loc
    type(a$D_t), intent(in) :: tree !< Parent box
    real(dp), intent(in)    :: r($D)
    integer, intent(in)     :: n_var     !< Number of variables
    integer, intent(in)     :: ivs(n_var)   !< Variables to interpolate
    real(dp)                :: vals(n_var)
    integer                 :: i, j, id, dix($D)
    real(dp)                :: f0(n_var), grad(n_var, $D)
    real(dp)                :: r_loc($D), dvec($D)
#if $D  == 3
    integer                 :: k
#endif
    type(a$D_loc_t)         :: loc  !< Location of cell

    loc   = a$D_get_loc(tree, r)
    id    = loc%id

    if (id <= af_no_box) then
       print *, "a$D_interp1: point outside domain", r
       stop
    end if

    r_loc = a$D_r_loc(tree, loc)
    dvec  = (r - r_loc) / tree%boxes(id)%dr
    i     = loc%ix(1)
    j     = loc%ix(2)
#if $D == 3
    k     = loc%ix(3)
#endif
    where (dvec > 0)
       dix = 0
    elsewhere
       dix = -1
    end where

#if $D == 2
    f0        = tree%boxes(id)%cc(i, j, ivs)
    grad(:, 1) = tree%boxes(id)%cc(i+dix(1)+1, j, ivs) - &
         tree%boxes(id)%cc(i+dix(1), j, ivs)
    grad(:, 2)   = tree%boxes(id)%cc(i, j+dix(2)+1, ivs) - &
         tree%boxes(id)%cc(i, j+dix(2), ivs)
#elif $D == 3
    f0        = tree%boxes(id)%cc(i, j, k, ivs)
    grad(:, 1)   = tree%boxes(id)%cc(i+dix(1)+1, j, k, ivs) - &
         tree%boxes(id)%cc(i+dix(1), j, k, ivs)
    grad(:, 2)   = tree%boxes(id)%cc(i, j+dix(2)+1, k, ivs) - &
         tree%boxes(id)%cc(i, j+dix(2), k, ivs)
    grad(:, 3)   = tree%boxes(id)%cc(i, j, k+dix(3)+1, ivs) - &
         tree%boxes(id)%cc(i, j, k+dix(3), ivs)
#endif
    vals = f0 + matmul(grad, dvec)
  end function a$D_interp1

  !> Using second order interpolation to get a value at r. The result is not
  !> continuous at cell boundaries.
  function a$D_interp2(tree, r, ivs, n_var) result(vals)
    use m_a$D_utils, only: a$D_get_loc, a$D_r_loc
    type(a$D_t), intent(in) :: tree !< Parent box
    real(dp), intent(in)    :: r($D)
    integer, intent(in)     :: n_var     !< Number of variables
    integer, intent(in)     :: ivs(n_var)   !< Variables to interpolate
    real(dp)                :: vals(n_var)
    integer                 :: i, j, id
    real(dp)                :: f0(n_var), grad(n_var, $D), second(n_var, $D)
    real(dp)                :: r_loc($D), dvec($D)
#if $D  == 3
    integer                 :: k
#endif
    type(a$D_loc_t)         :: loc  !< Location of cell

    loc   = a$D_get_loc(tree, r)
    id    = loc%id

    if (id <= af_no_box) then
       print *, "a$D_interp2: point outside domain", r
       stop
    end if

    r_loc = a$D_r_loc(tree, loc)
    dvec  = (r - r_loc) / tree%boxes(id)%dr
    i     = loc%ix(1)
    j     = loc%ix(2)
#if $D == 3
    k     = loc%ix(3)
#endif

#if $D == 2
    f0           = tree%boxes(id)%cc(i, j, ivs)
    grad(:, 1)   = 0.5_dp * (tree%boxes(id)%cc(i+1, j, ivs) - &
         tree%boxes(id)%cc(i-1, j, ivs))
    grad(:, 2)   = 0.5_dp * (tree%boxes(id)%cc(i, j+1, ivs) - &
         tree%boxes(id)%cc(i, j-1, ivs))
    second(:, 1) = tree%boxes(id)%cc(i-1, j, ivs) - &
         2 * f0 + tree%boxes(id)%cc(i+1, j, ivs)
    second(:, 2) = tree%boxes(id)%cc(i, j-1, ivs) - &
         2 * f0 + tree%boxes(id)%cc(i, j+1, ivs)
#elif $D == 3
    f0           = tree%boxes(id)%cc(i, j, k, ivs)
    grad(:, 1)   = 0.5_dp * (tree%boxes(id)%cc(i+1, j, k, ivs) - &
         tree%boxes(id)%cc(i-1, j, k, ivs))
    grad(:, 2)   = 0.5_dp * (tree%boxes(id)%cc(i, j+1, k, ivs) - &
         tree%boxes(id)%cc(i, j-1, k, ivs))
    grad(:, 3)   = 0.5_dp * (tree%boxes(id)%cc(i, j, k+1, ivs) - &
         tree%boxes(id)%cc(i, j, k-1, ivs))
    second(:, 1) = tree%boxes(id)%cc(i-1, j, k, ivs) - &
         2 * f0 + tree%boxes(id)%cc(i+1, j, k, ivs)
    second(:, 2) = tree%boxes(id)%cc(i, j-1, k, ivs) - &
         2 * f0 + tree%boxes(id)%cc(i, j+1, k, ivs)
    second(:, 3) = tree%boxes(id)%cc(i, j, k-1, ivs) - &
         2 * f0 + tree%boxes(id)%cc(i, j, k+1, ivs)
#endif
    vals = f0 + matmul(grad, dvec) + 0.5_dp * matmul(second, dvec**2)
  end function a$D_interp2

end module m_a$D_interp
