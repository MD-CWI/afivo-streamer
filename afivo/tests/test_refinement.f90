#include "../src/cpp_macros.h"
! Test the refinement procedure
program test_init
  use m_af_types
  use m_af_core

  implicit none

  type(af_t)       :: tree
  type(ref_info_t) :: ref_info
  integer          :: i, n_lvl
  integer          :: ixs(NDIM, 1), nbs(af_num_neighbors, 1)

  call af_add_cc_variable(tree, "phi")
  call af_add_fc_variable(tree, "flux")

  ! Call init with most options set
  call af_init(tree, n_cell=8, &
       dr = 1.0_dp, r_min=[DTIMES(0.0_dp)], &
       n_boxes=1000, coord=af_xyz, mem_limit_gb=1.0_dp)

  ixs = 1                    ! Box at 1,1
  nbs = 1                    ! Periodic
  call af_set_base(tree, 1, ixs, nbs)

  n_lvl = 4

  do i = 1, n_lvl
     call af_adjust_refinement(tree, refinement, ref_info)
  end do

  print *, af_num_boxes_used(tree), " == ", &
       (1-(2**NDIM)**(n_lvl+1)) / (1-(2**NDIM))

  do i = 1, n_lvl
     call af_adjust_refinement(tree, derefinement, ref_info)
  end do

  print *, af_num_boxes_used(tree), " == ", 1

contains

  subroutine refinement(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))

    if (box%lvl < 20) cell_flags = af_do_ref
  end subroutine refinement

  subroutine derefinement(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))
    cell_flags = af_rm_ref
  end subroutine derefinement

end program
