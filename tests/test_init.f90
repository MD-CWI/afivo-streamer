! Test the initialization of 2d and 3d trees
program test_init
  use m_a2_types
  use m_a2_core
  use m_a3_types
  use m_a3_core

  type(a2_t) :: tree_2d
  type(a3_t) :: tree_3d
  integer    :: ixs_2d(2, 1), nbs_2d(a2_num_neighbors, 1)
  integer    :: ixs_3d(3, 1), nbs_3d(a3_num_neighbors, 1)

  ! Call init with most options set
  call a2_init(tree_2d, n_cell=8, n_var_cell=1, n_var_face=1, &
       dr = 1.0_dp, r_min=[0.0_dp, 0.0_dp], lvl_limit=20, &
       n_boxes=100, coord=af_xyz, cc_names=["phi"], &
       fc_names=["flx"], mem_limit_gb=1.0_dp)

  call a3_init(tree_3d, n_cell=8, n_var_cell=1, n_var_face=1, &
       dr = 1.0_dp, r_min=[0.0_dp, 0.0_dp, 0.0_dp], lvl_limit=20, &
       n_boxes=100, coord=af_xyz, cc_names=["phi"], &
       fc_names=["flx"], mem_limit_gb=1.0_dp)

  ixs_2d = 1                    ! Box at 1,1
  nbs_2d = 1                    ! Periodic
  ixs_3d = 1                    ! Box at 1,1,1
  nbs_3d = 1                    ! Periodic

  call a2_set_base(tree_2d, ixs_2d, nbs_2d)
  call a3_set_base(tree_3d, ixs_3d, nbs_3d)

  call a2_print_info(tree_2d)
  call a3_print_info(tree_3d)

end program test_init
