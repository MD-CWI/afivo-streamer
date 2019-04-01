#include "../src/cpp_macros.h"
!> \example computational_domain_Xd.f90
!> Example showing how create different types of computational domains
program computational_domain_Xd
  use m_af_all

  implicit none

  integer    :: n_cell = 4 ! Boxes contain n_cell^dim cells
  type(af_t) :: tree       ! Will contain the quad/octree grid
  integer    :: grid_size(NDIM)
  real(dp)   :: domain_size(NDIM)
  logical    :: periodic(NDIM)

  ! Create mesh 1: two boxes along x-direction
  grid_size(:) = n_cell
  grid_size(1) = 2 * n_cell
  domain_size = 1.0_dp * grid_size

  call af_add_cc_variable(tree, "phi")
  call af_init(tree, n_cell, domain_size, grid_size)
  call af_write_vtk(tree, "computational_domain_" // DIMNAME // "_1", dir="output")
  call af_destroy(tree)

  ! Create mesh 2: two boxes along y-direction
  grid_size(:) = n_cell
  grid_size(2) = 2 * n_cell
  domain_size = 1.0_dp * grid_size

  call af_add_cc_variable(tree, "phi")
  call af_init(tree, n_cell, domain_size, grid_size)
  call af_write_vtk(tree, "computational_domain_" // DIMNAME // "_2", dir="output")
  call af_destroy(tree)

  ! Create mesh 3: Two boxes along x-direction that are fully periodic
  grid_size(:) = n_cell
  grid_size(1) = 2 * n_cell
  periodic(:)  = .true.
  domain_size  = 1.0_dp * grid_size

  call af_add_cc_variable(tree, "phi")
  call af_init(tree, n_cell, domain_size, grid_size, periodic=periodic)
  call af_write_vtk(tree, "computational_domain_" // DIMNAME // "_3", dir="output")
  call af_destroy(tree)

end program computational_domain_Xd