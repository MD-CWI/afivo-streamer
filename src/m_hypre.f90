!> This module contains an interface to Hypre, assuming Hypre is compiled with
!> OpenMP and without MPI
module m_hypre
  use iso_c_binding

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! Hypre CSR format
  integer, parameter :: HYPRE_PARCSR = 5555

  ! MPI stubs
  integer, parameter :: MPI_COMM_WORLD = 0
  integer, parameter :: num_procs = 1
  integer, parameter :: myid = 0


  public :: hypre_create_grid_2d
  public :: hypre_create_matrix_2d
  public :: hypre_create_vector
  public :: hypre_set_vector
  public :: hypre_prepare_solve
  public :: hypre_solve_smg

contains

  subroutine hypre_create_grid_2d(grid, nx)
    type(c_ptr), intent(out) :: grid
    integer, intent(in)         :: nx(2)
    integer, parameter          :: ndim = 2
    integer                     :: ierr
    integer                     :: ilo(ndim), ihi(ndim)

    ilo = 0
    ihi = nx-1

    ! Create an empty 2D grid object
    call HYPRE_StructGridCreate(MPI_COMM_WORLD, ndim, grid, ierr)

    ! Add a new box to the grid
    call HYPRE_StructGridSetExtents(grid, ilo, ihi, ierr)

    ! This is a collective call finalizing the grid assembly. The grid is now
    ! ``ready to be used''
    call HYPRE_StructGridAssemble(grid, ierr)
  end subroutine hypre_create_grid_2d

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

    ! ! Set vector to zero
    ! allocate(zeros(product(nx)))
    ! zeros(:) = 0.0_dp
    ! call HYPRE_StructVectorSetBoxValues(vec, ilo, nx-1, zeros, ierr)

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

  subroutine hypre_create_matrix_2d(A, grid, nx, coeffs)
    integer, parameter          :: ndim                          = 2
    integer, parameter          :: sym                           = 1
    integer, parameter          :: stencil_size                  = ndim + 1
    integer, parameter          :: stencil_indices(stencil_size) = [0, 1, 2]
    integer, parameter          :: offsets(2, stencil_size)      = &
         reshape([0, 0, 1, 0, 0, 1], [2, stencil_size])

    type(c_ptr), intent(out) :: A
    type(c_ptr), intent(in)  :: grid
    integer, intent(in)         :: nx(2)
    real(dp), intent(in)        :: coeffs(stencil_size, product(nx))
    type(c_ptr)              :: stencil
    integer                     :: i, ierr
    integer                     :: ilo(ndim), ihi(ndim)

    ilo = 0
    ihi = nx-1

    ! Create an empty stencil object (symmetric)
    call HYPRE_StructStencilCreate(ndim, stencil_size, stencil, ierr)

    ! Assign stencil entries */
    do i = 1, stencil_size
       call HYPRE_StructStencilSetElement(stencil, i-1, offsets(:, i), ierr)
    end do

    ! Create an empty matrix object */
    call HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, A, ierr)

    ! Use symmetric storage? */
    call HYPRE_StructMatrixSetSymmetric(A, sym, ierr)

    ! Indicate that the matrix coefficients are ready to be set */
    call HYPRE_StructMatrixInitialize(A, ierr)

    call HYPRE_StructMatrixSetBoxValues(A, ilo, ihi, stencil_size, stencil_indices, coeffs, ierr)

    call HYPRE_StructMatrixAssemble(A, ierr)
  end subroutine hypre_create_matrix_2d

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
    call HYPRE_StructSMGSetTol(solver, 1.0e-06_dp, ierr)
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

end module m_hypre
