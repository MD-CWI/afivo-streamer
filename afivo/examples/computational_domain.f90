#include "../src/cpp_macros.h"
!> \example computational_domain_Xd.f90
!> Example showing how create different types of computational domains
program computational_domain_Xd
  use m_af_all

  implicit none

  integer    :: n_cell   = 4    ! Boxes contain n_cell^dim cells
  real(dp)   :: dr       = 1.0_dp ! Coarse grid spacing
  type(af_t) :: tree             ! Will contain the quad/octree grid

  integer, parameter :: n_boxes = 2 ! How many boxes to add
  integer :: ix_list(NDIM, n_boxes)                ! Spatial indices of boxes
  integer :: nb_list(af_num_neighbors, n_boxes) ! Neighbor information

  ! Create mesh 1: two boxes along x-direction
#if NDIM == 2
  ix_list(:, 1) = [1, 1]        ! Box 1 at [1, 1]
  ix_list(:, 2) = [2, 1]        ! Box 2 at [2, 1]
#elif NDIM == 3
  ix_list(:, 1) = [1, 1, 1]     ! Box 1 at [1, 1, 1]
  ix_list(:, 2) = [2, 1, 1]     ! Box 2 at [2, 1, 1]
#endif

  call af_add_cc_variable(tree, "phi")
  call af_init(tree, n_cell, dr)
  call af_set_base(tree, n_boxes, ix_list)
  call af_write_vtk(tree, "computational_domain_" // DIMNAME // "_1", dir="output")
  call af_destroy(tree)

  ! Create mesh 2: two boxes along y-direction
#if NDIM == 2
  ix_list(:, 1) = [1, 1]        ! Box 1 at [1, 1]
  ix_list(:, 2) = [1, 2]        ! Box 2 at [1, 2]
#elif NDIM == 3
  ix_list(:, 1) = [1, 1, 1]     ! Box 1 at [1, 1, 1]
  ix_list(:, 2) = [1, 2, 1]     ! Box 2 at [1, 2, 1]
#endif

  call af_add_cc_variable(tree, "phi")
  call af_init(tree, n_cell, dr)
  call af_set_base(tree, n_boxes, ix_list)
  call af_write_vtk(tree, "computational_domain_" // DIMNAME // "_2", dir="output")
  call af_destroy(tree)

  ! Create mesh 3: Two boxes along x-direction that are fully periodic
#if NDIM == 2
  ix_list(:, 1) = [1, 1]        ! Box 1 at [1, 1]
  ix_list(:, 2) = [2, 1]        ! Box 2 at [2, 1]
#elif NDIM == 3
  ix_list(:, 1) = [1, 1, 1]     ! Box 1 at [1, 1, 1]
  ix_list(:, 2) = [2, 1, 1]     ! Box 2 at [2, 1, 1]
#endif

  nb_list(:, :) = af_no_box       ! Start with default value
  nb_list(af_neighb_lowx, 1) = 2 ! Box 1's lower x-neighbor is box 2
  nb_list(af_neighb_lowy, 1) = 1 ! Box 1 is periodic in y-direction
  nb_list(af_neighb_lowy, 2) = 2 ! Box 2 is periodic in y-direction
#if NDIM == 3
  nb_list(af_neighb_lowz, 1) = 1
  nb_list(af_neighb_lowz, 2) = 2
#endif

  call af_add_cc_variable(tree, "phi")
  call af_init(tree, n_cell, dr)
  call af_set_base(tree, n_boxes, ix_list, nb_list)
  call af_write_vtk(tree, "computational_domain_" // DIMNAME // "_3", dir="output")
  call af_destroy(tree)

end program computational_domain_Xd
