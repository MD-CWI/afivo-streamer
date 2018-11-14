#include "../src/cpp_macros.h"
! Test the refinement procedure
program test_init
  use m_af_types
  use m_af_core
  use m_af_ghostcell
  use m_af_utils

  implicit none

  type(af_t)       :: tree
  type(ref_info_t) :: ref_info
  integer          :: i, n_lvl
  integer          :: ixs(NDIM, 1), nbs(af_num_neighbors, 1)

  call af_add_cc_variable(tree, "phi")
  call af_init(tree, n_cell=8, dr = 1.0_dp)

  ixs = 1                    ! Box at 1,1
  nbs = -1                   ! Boundary condition indicated by -1
  call af_set_base(tree, 1, ixs, nbs)

  n_lvl = 4

  do i = 1, n_lvl
     call af_adjust_refinement(tree, refinement, ref_info)
  end do

  call af_loop_box(tree, init)

  ! Should set all ghost cells to zero
  call af_gc_tree(tree, 1, af_gc_interp, af_bc_neumann_zero)

  call af_loop_box(tree, check_ghostcell)

  print *, "Success"

contains

  subroutine refinement(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(DTIMES(box%n_cell))
    if (all(box%ix == 1)) then
       cell_flags = af_do_ref
    else
       cell_flags = af_keep_ref
    end if
  end subroutine refinement

  subroutine init(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc

    nc                    = box%n_cell
    box%cc(DTIMES(1:nc), 1) = 0
  end subroutine init

  subroutine check_ghostcell(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc
    real(dp)                   :: tmp

    nc  = box%n_cell
#if NDIM == 2
    tmp = sum(box%cc(1:nc, 0, 1)) + &
         sum(box%cc(1:nc, nc+1, 1)) + &
         sum(box%cc(0, 1:nc, 1)) + &
         sum(box%cc(nc+1, 1:nc, 1))
#elif NDIM == 3
    tmp = sum(box%cc(1:nc, 1:nc, 0, 1)) + &
         sum(box%cc(1:nc, 1:nc, nc+1, 1)) + &
         sum(box%cc(1:nc, 0, 1:nc, 1)) + &
         sum(box%cc(1:nc, nc+1, 1:nc, 1)) + &
         sum(box%cc(0, 1:nc, 1:nc, 1)) + &
         sum(box%cc(nc+1, 1:nc, 1:nc, 1))
#endif

    if (tmp > 0 .or. tmp < 0) stop "Wrong ghostcell value"
  end subroutine check_ghostcell

end program
