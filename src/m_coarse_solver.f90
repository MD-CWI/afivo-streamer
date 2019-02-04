#include "../src/cpp_macros.h"
!> Module to solve elliptic PDEs on the coarse grid. This module contains an
!> interface to Hypre, assuming Hypre is compiled with OpenMP and without MPI
module m_coarse_solver
  use m_af_types
  use m_mg_types

  implicit none
  private

  ! Hypre CSR format
  integer, parameter :: HYPRE_PARCSR = 5555

  ! Assume matrices are symmetric
  integer, parameter :: symmetric_matrix = 1

  ! Size of the stencil (only direct neighbors)
  integer, parameter :: stencil_size = NDIM + 1

  ! Order in which the stencil is stored
#if NDIM == 2
  integer, parameter :: stencil_indices(stencil_size) = [0, 1, 2]
#elif NDIM == 3
  integer, parameter :: stencil_indices(stencil_size) = [0, 1, 2, 3]
#endif

  ! Offsets for the stencil elements
#if NDIM == 2
  integer, parameter :: offsets(2, stencil_size) = &
       reshape([0, 0, 1, 0, 0, 1], [2, stencil_size])
#elif NDIM == 3
  integer, parameter :: offsets(3, stencil_size) = &
       reshape([0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1], [3, stencil_size])
#endif

  ! MPI stubs
  integer, parameter :: MPI_COMM_WORLD = 0
  integer, parameter :: num_procs = 1
  integer, parameter :: myid = 0

  public :: coarse_solver_initialize
  ! public :: hypre_create_grid_2d
  ! public :: hypre_create_matrix_2d
  ! public :: hypre_create_vector
  ! public :: hypre_set_vector
  ! public :: hypre_prepare_solve
  public :: hypre_solve_smg

contains

  subroutine coarse_solver_initialize(tree, mg)
    type(af_t), intent(in)    :: tree !< Tree to do multigrid on
    type(mg_t), intent(inout) :: mg
    integer                   :: n, n_boxes, nc, id, ierr
    integer                   :: nx(NDIM), ilo(NDIM), ihi(NDIM)
    real(dp)                  :: symm_coeffs(stencil_size, tree%n_cell**NDIM)
    real(dp)                  :: full_coeffs(2*NDIM+1, DTIMES(tree%n_cell))

    nc                   = tree%n_cell
    nx                   = tree%coarse_grid_size(1:NDIM)
    mg%csolver%grid_size = tree%coarse_grid_size(1:NDIM)
    n_boxes              = size(tree%lvls(1)%ids)

    if (any(nx == -1)) &
         error stop "coarse_solver_initialize: coarse_grid_size not set"

    ! Construct grid and vectors
    print *, "Create grid"
    call hypre_create_grid(mg%csolver%grid, nx, tree%periodic)
    call hypre_create_vector(mg%csolver%grid, nx, mg%csolver%phi)
    call hypre_create_vector(mg%csolver%grid, nx, mg%csolver%rhs)

    ! Construct the matrix
    print *, "Create matrix"
    call hypre_create_matrix(mg%csolver%matrix, mg%csolver%grid)

    if (.not. associated(mg%box_stencil)) &
         error stop "coarse_solver_initialize: mg%box_stencil not set"

    allocate(mg%csolver%bc_to_rhs(DTIMES(nc), n_boxes))

    do n = 1, size(tree%lvls(1)%ids)
       id  = tree%lvls(1)%ids(n)
       call mg%box_stencil(tree%boxes(id), mg, full_coeffs, &
            mg%csolver%bc_to_rhs(DTIMES(:), n))
       ! TODO: check symmetry

       symm_coeffs(1, :) = pack(full_coeffs(1, DTIMES(:)), .true.)
       symm_coeffs(2, :) = pack(full_coeffs(3, DTIMES(:)), .true.)
       symm_coeffs(3, :) = pack(full_coeffs(5, DTIMES(:)), .true.)
#if NDIM == 3
       symm_coeffs(4, :) = pack(full_coeffs(7, DTIMES(:)), .true.)
#endif

       ilo = (tree%boxes(id)%ix - 1) * nc + 1
       ihi = ilo + nc - 1
       call HYPRE_StructMatrixSetBoxValues(mg%csolver%matrix, ilo, ihi, &
            stencil_size, stencil_indices, symm_coeffs, ierr)
    end do

    print *, "Assemble matrix"
    call HYPRE_StructMatrixAssemble(mg%csolver%matrix, ierr)

    print *, "Prepare solve"
    call hypre_prepare_solve(mg%csolver%solver, mg%csolver%matrix, &
         mg%csolver%rhs, mg%csolver%phi)

  end subroutine coarse_solver_initialize

  subroutine hypre_create_grid(grid, nx, periodic)
    type(c_ptr), intent(out) :: grid
    integer, intent(in)      :: nx(NDIM) !< Size of grid
    logical, intent(in)      :: periodic(NDIM) !< Whether the dimension is periodic
    integer                  :: ierr, ilo(NDIM)
    integer                  :: period(NDIM)

    ilo(:) = 1

    ! Create an empty 2D grid object
    call HYPRE_StructGridCreate(MPI_COMM_WORLD, NDIM, grid, ierr)

    ! Add a new box to the grid
    call HYPRE_StructGridSetExtents(grid, ilo, nx, ierr)

    ! Set periodic
    where (periodic)
       period = nx
    elsewhere
       period = 0
    end where

    call HYPRE_StructGridSetPeriodic(grid, period, ierr)

    ! This is a collective call finalizing the grid assembly. The grid is now
    ! ``ready to be used''
    call HYPRE_StructGridAssemble(grid, ierr)
  end subroutine hypre_create_grid

  subroutine hypre_create_vector(grid, nx, vec)
    type(c_ptr), intent(in)  :: grid
    integer, intent(in)      :: nx(:)
    type(c_ptr), intent(out) :: vec
    real(dp), allocatable    :: zeros(:)
    integer                  :: ilo(size(nx))
    integer                  :: ierr

    ilo = 0

    ! Create an empty vector object */
    call HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, vec, ierr)

    ! Indicate that the vector coefficients are ready to be set */
    call HYPRE_StructVectorInitialize(vec, ierr)

    ! Set vector to zero
    allocate(zeros(product(nx)))
    zeros(:) = 0.0_dp
    call HYPRE_StructVectorSetBoxValues(vec, ilo, nx-1, zeros, ierr)

    ! Finalize the construction of the vector before using
    call HYPRE_StructVectorAssemble(vec, ierr)
  end subroutine hypre_create_vector

  subroutine hypre_set_vector(vec, nx, val)
    type(c_ptr), intent(in) :: vec
    integer, intent(in)     :: nx(:)
    real(dp), intent(in)    :: val(product(nx))
    integer                 :: ilo(size(nx))

    ilo = 0

    call HYPRE_StructVectorSetBoxValues(vec, ilo, nx-1, val);
  end subroutine hypre_set_vector

  subroutine hypre_create_matrix(A, grid)
    type(c_ptr), intent(out) :: A
    type(c_ptr), intent(in)  :: grid
    type(c_ptr)              :: stencil
    integer                  :: i, ierr

    ! Create an empty stencil object (symmetric)
    call HYPRE_StructStencilCreate(NDIM, stencil_size, stencil, ierr)

    ! Assign stencil entries */
    do i = 1, stencil_size
       call HYPRE_StructStencilSetElement(stencil, i-1, offsets(:, i), ierr)
    end do

    ! Create an empty matrix object */
    call HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, A, ierr)

    ! Use symmetric storage? */
    call HYPRE_StructMatrixSetSymmetric(A, symmetric_matrix, ierr)

    ! Indicate that the matrix coefficients are ready to be set */
    call HYPRE_StructMatrixInitialize(A, ierr)

  end subroutine hypre_create_matrix

  ! Prepare to solve the system. The coefficient data in b and x is ignored
  ! here, but information about the layout of the data may be used.
  subroutine hypre_prepare_solve(solver, A, b, x)
    type(c_ptr), intent(out) :: solver
    type(c_ptr), intent(in)  :: A, b, x
    integer, parameter       :: n_pre = 1, n_post = 1
    integer                  :: ierr

    call HYPRE_StructSMGCreate(MPI_COMM_WORLD, solver, ierr)
    call HYPRE_StructSMGSetMemoryUse(solver, 0, ierr)
    call HYPRE_StructSMGSetMaxIter(solver, 50, ierr)
    call HYPRE_StructSMGSetTol(solver, 1.0e-09_dp, ierr)
    call HYPRE_StructSMGSetRelChange(solver, 0, ierr)
    call HYPRE_StructSMGSetNumPreRelax(solver, n_pre, ierr)
    call HYPRE_StructSMGSetNumPostRelax(solver, n_post, ierr)
    call HYPRE_StructSMGSetPrintLevel(solver, 1, ierr)
    call HYPRE_StructSMGSetLogging(solver, 1, ierr)
    call HYPRE_StructSMGSetup(solver, A, b, x, ierr)
  end subroutine hypre_prepare_solve

  ! Solve the system A x = b
  subroutine hypre_solve_smg(solver, A, b, x)
    integer :: ierr
    type(c_ptr), intent(in)  :: solver, A, b, x
    call HYPRE_StructSMGSolve(solver, A, b, x, ierr)
  end subroutine hypre_solve_smg

end module m_coarse_solver
