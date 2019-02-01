#include "../src/cpp_macros.h"
!> \example hypre_coupling.f90
!>
!> This program tests the coupling with Hypre to solve coarse grid equations
program hypre_coupling
  use m_hypre
  use iso_c_binding

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: nx(2) = [2, 2]

  type(c_ptr) :: grid, A, x, b, solver
  real(dp)    :: vals(product(nx))
  real(dp)    :: matrix_vals(3, nx(1), nx(2))
  real(dp)    :: coeffs(3, product(nx))
  integer     :: i, ierr

  call hypre_create_grid_2d(grid, nx)
  call hypre_create_vector(grid, nx, x)
  call hypre_create_vector(grid, nx, b)

  ! at 1,1
  coeffs(:, 1) = [-2.0_dp, 1.0_dp, 1.0_dp]
  ! at 2,1
  coeffs(:, 2) = [-2.0_dp, 0.0_dp, 1.0_dp]
  ! at 1,2
  coeffs(:, 3) = [-2.0_dp, 1.0_dp, 0.0_dp]
  ! at 2,2
  coeffs(:, 4) = [-2.0_dp, 0.0_dp, 0.0_dp]
  coeffs(:, :) = 0.0_dp
  coeffs(1, :) = 1.0_dp
  ! coeffs(2:3, :) = 1.0_dp

  call hypre_create_matrix_2d(A, grid, nx, coeffs)

  call HYPRE_StructVectorGetBoxValues(x, [0, 0], nx-1, vals, ierr)
  print *, vals

  call HYPRE_StructVectorSetBoxValues(b, [0, 0], nx-1, [(1.0_dp, i=1,product(nx))], ierr)

  call HYPRE_StructVectorGetBoxValues(b, [0, 0], nx-1, vals, ierr)
  print *, vals

  call HYPRE_StructMatrixGetBoxValues(A, [0, 0], nx-1, 3, [0, 1, 2], matrix_vals, ierr)
  print *, matrix_vals
  call hypre_prepare_solve(solver, A, b, x)

  call hypre_solve_smg(solver, A, b, x)

  call HYPRE_StructVectorGetBoxValues(x, [0, 0], nx-1, vals, ierr)
  print *, vals

  call HYPRE_StructGridDestroy(grid, ierr);
  ! call HYPRE_StructStencilDestroy(stencil, ierr);
  call HYPRE_StructMatrixDestroy(A, ierr);
  call HYPRE_StructVectorDestroy(b, ierr);
  call HYPRE_StructVectorDestroy(x, ierr);

end program hypre_coupling

