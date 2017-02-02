!> \example computational_domain_$Dd.f90
!> Example showing how create different types of computational domains
program computational_domain_$Dd
  use m_a$D_all

  implicit none

  integer    :: n_cell   = 4    ! Boxes contain n_cell^dim cells
  real(dp)   :: dr       = 1.0_dp ! Coarse grid spacing
  integer    :: n_var_cc = 1      ! Number of cell-centered variables
  integer    :: n_var_fc = 0      ! Number of face-centered variables
  type(a$D_t) :: tree             ! Will contain the quad/octree grid

  integer, parameter :: n_boxes = 2 ! How many boxes to add
  integer :: ix_list($D, n_boxes)                ! Spatial indices of boxes
  integer :: nb_list(a$D_num_neighbors, n_boxes) ! Neighbor information

  ! Create mesh 1: two boxes along x-direction
#if $D == 2
  ix_list(:, 1) = [1, 1]        ! Box 1 at [1, 1]
  ix_list(:, 2) = [2, 1]        ! Box 2 at [2, 1]
#elif $D == 3
  ix_list(:, 1) = [1, 1, 1]     ! Box 1 at [1, 1, 1]
  ix_list(:, 2) = [2, 1, 1]     ! Box 2 at [2, 1, 1]
#endif
  call a$D_init(tree, n_cell, n_var_cc, n_var_fc, dr)
  call a$D_set_base(tree, n_boxes, ix_list)
  call a$D_write_vtk(tree, "computational_domain_$Dd_1", dir="output")
  call a$D_destroy(tree)

  ! Create mesh 2: two boxes along y-direction
#if $D == 2
  ix_list(:, 1) = [1, 1]        ! Box 1 at [1, 1]
  ix_list(:, 2) = [1, 2]        ! Box 2 at [1, 2]
#elif $D == 3
  ix_list(:, 1) = [1, 1, 1]     ! Box 1 at [1, 1, 1]
  ix_list(:, 2) = [1, 2, 1]     ! Box 2 at [1, 2, 1]
#endif
  nb_list(:, :) = af_no_box     ! Afivo resolves the connectivity
  call a$D_init(tree, n_cell, n_var_cc, n_var_fc, dr)
  call a$D_set_base(tree, n_boxes, ix_list)
  call a$D_write_vtk(tree, "computational_domain_$Dd_2", dir="output")
  call a$D_destroy(tree)

  ! Create mesh 3: Two boxes along x-direction that are fully periodic
#if $D == 2
  ix_list(:, 1) = [1, 1]        ! Box 1 at [1, 1]
  ix_list(:, 2) = [2, 1]        ! Box 2 at [2, 1]
#elif $D == 3
  ix_list(:, 1) = [1, 1, 1]     ! Box 1 at [1, 1, 1]
  ix_list(:, 2) = [2, 1, 1]     ! Box 2 at [2, 1, 1]
#endif
  nb_list(:, :) = af_no_box       ! Start with default value
  nb_list(a$D_neighb_lowx, 1) = 2 ! Box 1's lower x-neighbor is box 2
  nb_list(a$D_neighb_lowy, 1) = 1 ! Box 1 is periodic in y-direction
  nb_list(a$D_neighb_lowy, 2) = 2 ! Box 2 is periodic in y-direction
#if $D == 3
  nb_list(a$D_neighb_lowz, 1) = 1
  nb_list(a$D_neighb_lowz, 2) = 2
#endif
  call a$D_init(tree, n_cell, n_var_cc, n_var_fc, dr)
  call a$D_set_base(tree, n_boxes, ix_list, nb_list)
  call a$D_write_vtk(tree, "computational_domain_$Dd_3", dir="output")
  call a$D_destroy(tree)

end program computational_domain_$Dd
