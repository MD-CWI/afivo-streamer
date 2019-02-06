#include "../src/cpp_macros.h"
!> \example test_reduction.f90
!> This example shows the basic reduction functionality of m_af_types.
program test_reduction
  use m_af_types
  use m_af_core
  use m_af_utils

  implicit none

  type(af_t)         :: tree
  integer            :: i
  integer, parameter :: box_size     = 8
  integer            :: i_phi
  type(ref_info_t)   :: ref_info
  real(dp)           :: dr, max_val, min_val
  type(af_loc_t)     :: max_loc, min_loc

  dr = 2 * acos(-1.0_dp) / box_size ! 2 * pi / box_size

  call af_add_cc_variable(tree, "phi", ix=i_phi)

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! Number of cells per coordinate in a box
       dr)             ! Distance between cells on base level

  ! Create the base mesh
  call af_set_coarse_grid(tree, [DTIMES(box_size)], [DTIMES(.true.)])

  ! Set variables on base
  call af_loop_box(tree, set_random_values)

  do i = 1, 16
     print *, "i = ", i, "highest_id", tree%highest_id
     call af_adjust_refinement(tree, ref_func, ref_info, 0)
     call af_loop_box(tree, set_random_values)

     call af_tree_max_cc(tree, i_phi, max_val)
     call af_tree_min_cc(tree, i_phi, min_val)
     print *, "1 - max/min", max_val, min_val
     call af_reduction(tree, box_max, max_ab, -huge(1.0_dp), max_val)
     call af_reduction(tree, box_min, min_ab,  huge(1.0_dp), min_val)
     print *, "2 - max/min", max_val, min_val
     call af_reduction_loc(tree, i_phi, box_max_ix, max_ab, -huge(1.0_dp), &
          max_val, max_loc)
     call af_reduction_loc(tree, i_phi, box_min_ix, min_ab,  huge(1.0_dp), &
          min_val, min_loc)
     print *, "3 - max/min", max_val, min_val
#if NDIM == 2
     print *, "4 - max/min", tree%boxes(max_loc%id)%cc(max_loc%ix(1), &
          max_loc%ix(2), i_phi), tree%boxes(min_loc%id)%cc(min_loc%ix(1), &
          min_loc%ix(2), i_phi)
#elif NDIM == 3
     print *, "4 - max/min", tree%boxes(max_loc%id)%cc(max_loc%ix(1), &
          max_loc%ix(2), max_loc%ix(3), i_phi), &
          tree%boxes(min_loc%id)%cc(min_loc%ix(1), &
          min_loc%ix(2), max_loc%ix(3), i_phi)
#endif
  end do

  call af_destroy(tree)

contains

  real(dp) function box_max(box)
    type(box_t), intent(in) :: box
    box_max = maxval(box%cc(DTIMES(1:box%n_cell), i_phi))
  end function box_max

  subroutine box_max_ix(box, iv, val, ix)
    type(box_t), intent(in) :: box
    integer, intent(in) :: iv
    real(dp), intent(out) :: val
    integer, intent(out) :: ix(NDIM)
    ix = maxloc(box%cc(DTIMES(1:box%n_cell), iv))
#if NDIM == 2
    val = box%cc(ix(1), ix(2), iv)
#elif NDIM == 3
    val = box%cc(ix(1), ix(2), ix(3), iv)
#endif
  end subroutine box_max_ix

  subroutine box_min_ix(box, iv, val, ix)
    type(box_t), intent(in) :: box
    integer, intent(in) :: iv
    real(dp), intent(out) :: val
    integer, intent(out) :: ix(NDIM)
    ix = minloc(box%cc(DTIMES(1:box%n_cell), iv))
#if NDIM == 2
    val = box%cc(ix(1), ix(2), iv)
#elif NDIM == 3
    val = box%cc(ix(1), ix(2), ix(3), iv)
#endif
  end subroutine box_min_ix

  real(dp) function box_min(box)
    type(box_t), intent(in) :: box
    box_min = minval(box%cc(DTIMES(1:box%n_cell), i_phi))
  end function box_min

  real(dp) function max_ab(a,b)
    real(dp), intent(in) :: a, b
    max_ab = max(a,b)
  end function max_ab

  real(dp) function min_ab(a,b)
    real(dp), intent(in) :: a, b
    min_ab = min(a,b)
  end function min_ab

  subroutine ref_func(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))

    if (all(box%r_min < 0.4_dp) .and. box%lvl < 10) then
       cell_flags = af_do_ref
    else
       cell_flags = af_rm_ref
    end if
  end subroutine ref_func

  subroutine set_random_values(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc

    nc = box%n_cell
    box%cc(DTIMES(1:nc), i_phi) = sum(box%ix)
  end subroutine set_random_values

end program test_reduction
