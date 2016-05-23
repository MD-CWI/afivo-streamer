! Test the refinement procedure
program test_init
  use m_a2_types
  use m_a2_core
  use m_a3_types
  use m_a3_core

  implicit none

  type(a2_t)       :: tree_2d
  type(a3_t)       :: tree_3d
  type(ref_info_t) :: ref_info
  integer          :: i, n_lvl
  integer          :: ixs_2d(2, 1), nbs_2d(a2_num_neighbors, 1)
  integer          :: ixs_3d(3, 1), nbs_3d(a3_num_neighbors, 1)

  ! Call init with most options set
  call a2_init(tree_2d, n_cell=8, n_var_cell=1, n_var_face=1, &
       dr = 1.0_dp, r_min=[0.0_dp, 0.0_dp], lvl_limit=20, &
       n_boxes=1000, coord=af_xyz, cc_names=["phi"], &
       fc_names=["flx"], mem_limit_gb=1.0_dp)

  call a3_init(tree_3d, n_cell=8, n_var_cell=1, n_var_face=1, &
       dr = 1.0_dp, r_min=[0.0_dp, 0.0_dp, 0.0_dp], lvl_limit=20, &
       n_boxes=1000, coord=af_xyz, cc_names=["phi"], &
       fc_names=["flx"], mem_limit_gb=1.0_dp)

  ixs_2d = 1                    ! Box at 1,1
  nbs_2d = 1                    ! Periodic
  ixs_3d = 1                    ! Box at 1,1,1
  nbs_3d = 1                    ! Periodic

  call a2_set_base(tree_2d, ixs_2d, nbs_2d)
  call a3_set_base(tree_3d, ixs_3d, nbs_3d)

  n_lvl = 4

  do i = 1, n_lvl
     call a2_adjust_refinement(tree_2d, refinement_2d, ref_info)
     call a3_adjust_refinement(tree_3d, refinement_3d, ref_info)
  end do

  print *, a2_num_boxes_used(tree_2d), " == ", (1-4**(n_lvl+1)) / (1-4)
  print *, a3_num_boxes_used(tree_3d), " == ", (1-8**(n_lvl+1)) / (1-8)

  do i = 1, n_lvl
     call a2_adjust_refinement(tree_2d, derefinement_2d, ref_info)
     call a3_adjust_refinement(tree_3d, derefinement_3d, ref_info)
  end do

  print *, a2_num_boxes_used(tree_2d), " == ", 1
  print *, a3_num_boxes_used(tree_3d), " == ", 1

contains

  subroutine refinement_2d(boxes, id, ref_flag)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flag
    ref_flag = af_do_ref
  end subroutine refinement_2d

  subroutine derefinement_2d(boxes, id, ref_flag)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flag
    ref_flag = af_rm_ref
  end subroutine derefinement_2d

  subroutine refinement_3d(boxes, id, ref_flag)
    type(box3_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flag
    ref_flag = af_do_ref
  end subroutine refinement_3d

  subroutine derefinement_3d(boxes, id, ref_flag)
    type(box3_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flag
    ref_flag = af_rm_ref
  end subroutine derefinement_3d
end program
