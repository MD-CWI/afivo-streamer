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

  ! Semi-coarsening multigrid (more robust and expensive)
  integer, parameter, public :: coarse_solver_hypre_smg = 1

  ! Multigrid with point-wise smoother (fast/cheap)
  integer, parameter, public :: coarse_solver_hypre_pfmg = 2

    ! Size of the stencil (only direct neighbors)
  integer, parameter :: max_stencil_size = 2*NDIM + 1

  ! Offsets for the stencil elements
#if NDIM == 2
  integer, parameter :: stencil_offsets(2, max_stencil_size) = reshape([0, 0, &
       -1, 0, 1, 0, &
       0, -1, 0, 1], [2, max_stencil_size])
#elif NDIM == 3
  integer, parameter :: stencil_offsets(3, max_stencil_size) = reshape([0, 0, 0, &
       -1, 0, 0, 1, 0, 0, &
       0, -1, 0, 0, 1, 0, &
       0, 0, -1, 0, 0, 1], [3, max_stencil_size])
#endif

  ! MPI stubs
  integer, parameter :: MPI_COMM_WORLD = 0
  integer, parameter :: num_procs = 1
  integer, parameter :: myid = 0

  interface
     subroutine HYPRE_StructMatrixSetBoxValues(matrix, ilower, iupper, nentries, &
          entries, values, ierr)
       import
       type(c_ptr), intent(in) :: matrix
       integer, intent(in)     :: ilower(*), iupper(*)
       integer, intent(in)     :: nentries, entries(*)
       real(dp), intent(in)    :: values(*)
       integer, intent(out)    :: ierr
     end subroutine HYPRE_StructMatrixSetBoxValues
  end interface

  public :: coarse_solver_initialize
  public :: coarse_solver_destroy
  public :: coarse_solver_set_rhs_phi
  public :: coarse_solver_get_phi
  public :: coarse_solver

contains

  !> Initialize the coarse grid solver
  subroutine coarse_solver_initialize(tree, mg)
    type(af_t), intent(in)    :: tree !< Tree to do multigrid on
    type(mg_t), intent(inout) :: mg
    integer                   :: n, n_boxes, nc, id, IJK
    integer                   :: cnt, stencil_size, ierr
    integer                   :: nx(NDIM), ilo(NDIM), ihi(NDIM)
    real(dp), allocatable     :: coeffs(:, :)
    real(dp)                  :: full_coeffs(2*NDIM+1, DTIMES(tree%n_cell))
    integer, allocatable      :: stencil_ix(:)
    integer, parameter        :: zero_to_n(max_stencil_size) = [(i, i=0, max_stencil_size-1)]

    if (.not. tree%ready) error stop "coarse_solver_initialize: tree not ready"
    if (.not. associated(mg%box_stencil)) &
         error stop "coarse_solver_initialize: mg%box_stencil not set"

    nc                   = tree%n_cell
    nx                   = tree%coarse_grid_size(1:NDIM)
    n_boxes              = size(tree%lvls(1)%ids)

    if (any(nx == -1)) &
         error stop "coarse_solver_initialize: coarse_grid_size not set"

    ! Construct grid and vectors
    call hypre_create_grid(mg%csolver%grid, nx, tree%periodic)
    call hypre_create_vector(mg%csolver%grid, nx, mg%csolver%phi)
    call hypre_create_vector(mg%csolver%grid, nx, mg%csolver%rhs)

    if (tree%coord_t == af_cyl) then
       ! The symmetry option does not seem to work well with axisymmetric problems
       mg%csolver%symmetric = 0
    end if

    if (mg%csolver%symmetric == 1) then
       stencil_size = NDIM + 1
       allocate(stencil_ix(stencil_size))
#if NDIM == 2
       stencil_ix = [1, 3, 5]
#elif NDIM == 3
       stencil_ix = [1, 3, 5, 7]
#endif
    else
       stencil_size = 2*NDIM + 1
       allocate(stencil_ix(stencil_size))
#if NDIM == 2
       stencil_ix = [1, 2, 3, 4, 5]
#elif NDIM == 3
       stencil_ix = [1, 2, 3, 4, 5, 6, 7]
#endif
    end if

    ! Construct the matrix
    call hypre_create_matrix(mg%csolver%matrix, mg%csolver%grid, &
         stencil_size, stencil_offsets(:, stencil_ix), mg%csolver%symmetric)

    allocate(mg%csolver%bc_to_rhs(nc**(NDIM-1), af_num_neighbors, n_boxes))
    allocate(coeffs(stencil_size, tree%n_cell**NDIM))

    do n = 1, size(tree%lvls(1)%ids)
       id  = tree%lvls(1)%ids(n)
       call mg%box_stencil(tree%boxes(id), mg, full_coeffs, &
            mg%csolver%bc_to_rhs(:, :, n))

       cnt = 0
       do KJI_DO(1, nc)
          cnt = cnt + 1
          coeffs(:, cnt) = full_coeffs(stencil_ix, IJK)
       end do; CLOSE_DO

       ilo = (tree%boxes(id)%ix - 1) * nc + 1
       ihi = ilo + nc - 1
       call HYPRE_StructMatrixSetBoxValues(mg%csolver%matrix, ilo, ihi, &
            stencil_size, zero_to_n(1:stencil_size), coeffs, ierr)
       if (ierr /= 0) error stop "HYPRE_StructMatrixSetBoxValues failed"
    end do

    call HYPRE_StructMatrixAssemble(mg%csolver%matrix, ierr)

    if (mg%csolver%solver_type <= 0) then
       mg%csolver%solver_type = coarse_solver_hypre_pfmg
    end if

    call hypre_prepare_solve(mg%csolver)

  end subroutine coarse_solver_initialize

  subroutine coarse_solver_destroy(cs)
    type(coarse_solve_t), intent(inout) :: cs
    integer                             :: ierr

    call HYPRE_StructGridDestroy(cs%grid, ierr)
    call HYPRE_StructMatrixDestroy(cs%matrix, ierr)
    call HYPRE_StructVectorDestroy(cs%rhs, ierr)
    call HYPRE_StructVectorDestroy(cs%phi, ierr)

    select case (cs%solver_type)
    case (coarse_solver_hypre_smg)
       call HYPRE_StructSMGDestroy(cs%solver, ierr)
    case (coarse_solver_hypre_pfmg)
       call HYPRE_StructPFMGDestroy(cs%solver, ierr)
    case default
       error stop "hypre_solver_destroy: unknown solver type"
    end select
  end subroutine coarse_solver_destroy

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

  !> Create a Hypre vector
  subroutine hypre_create_vector(grid, nx, vec)
    type(c_ptr), intent(in)  :: grid
    integer, intent(in)      :: nx(:)
    type(c_ptr), intent(out) :: vec
    integer                  :: ierr

    ! Create an empty vector object */
    call HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, vec, ierr)

    ! Indicate that the vector coefficients are ready to be set */
    call HYPRE_StructVectorInitialize(vec, ierr)

    ! Finalize the construction of the vector before using
    call HYPRE_StructVectorAssemble(vec, ierr)
  end subroutine hypre_create_vector

  !> Set the right-hand side and copy phi from the tree. Also move the boundary
  !> conditions for phi to the rhs.
  subroutine coarse_solver_set_rhs_phi(tree, mg)
    use m_af_ghostcell, only: af_gc_get_boundary_coords
    type(af_t), intent(in) :: tree
    type(mg_t), intent(in) :: mg
    integer                :: n, nb, nc, id, bc_type, ierr
    integer                :: ilo(NDIM), ihi(NDIM)
    real(dp)               :: tmp(DTIMES(tree%n_cell))
    real(dp)               :: coords(NDIM, tree%n_cell**(NDIM-1))
    real(dp)               :: bc_val(tree%n_cell**(NDIM-1))

    nc = tree%n_cell

    do n = 1, size(tree%lvls(1)%ids)
       id  = tree%lvls(1)%ids(n)

       tmp = tree%boxes(id)%cc(DTIMES(1:nc), mg%i_rhs)

       ! Add contribution of boundary conditions to rhs
       do nb = 1, af_num_neighbors
          if (tree%boxes(id)%neighbors(nb) < af_no_box) then
             ! Put the boundary condition into the ghost cells of mg%i_phi
             call af_gc_get_boundary_coords(tree%boxes(id), nb, coords)
             call mg%sides_bc(tree%boxes(id), nb, mg%i_phi, coords, bc_val, bc_type)

             ! Get index range near neighbor
             call af_get_index_bc_inside(nb, nc, ilo, ihi)

             ! Use the stored arrays mg%csolver%bc_to_rhs to convert the value
             ! at the boundary to the rhs
#if NDIM == 2
             tmp(ilo(1):ihi(1), ilo(2):ihi(2)) = &
                  tmp(ilo(1):ihi(1), ilo(2):ihi(2)) + &
                  reshape(mg%csolver%bc_to_rhs(:, nb, n) * bc_val, [ihi - ilo + 1])
#elif NDIM == 3
             tmp(ilo(1):ihi(1), ilo(2):ihi(2), ilo(3):ihi(3)) = &
                  tmp(ilo(1):ihi(1), ilo(2):ihi(2), ilo(3):ihi(3)) + &
                  reshape(mg%csolver%bc_to_rhs(:, nb, n) * pack(bc_val, .true.), &
                  [ihi - ilo + 1])
#endif
          end if
       end do

       ilo = (tree%boxes(id)%ix - 1) * nc + 1
       ihi = ilo + nc - 1
       call HYPRE_StructVectorSetBoxValues(mg%csolver%rhs, ilo, ihi, &
            pack(tmp, .true.), ierr)

       tmp = tree%boxes(id)%cc(DTIMES(1:nc), mg%i_phi)
       call HYPRE_StructVectorSetBoxValues(mg%csolver%phi, ilo, ihi, &
            pack(tmp, .true.), ierr)
    end do
  end subroutine coarse_solver_set_rhs_phi

  !> Copy solution from coarse solver to the tree
  subroutine coarse_solver_get_phi(tree, mg)
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(in)    :: mg
    integer                   :: n, nc, id, ierr
    integer                   :: ilo(NDIM), ihi(NDIM)
    real(dp)                  :: tmp(tree%n_cell**NDIM)

    nc = tree%n_cell

    do n = 1, size(tree%lvls(1)%ids)
       id  = tree%lvls(1)%ids(n)
       ilo = (tree%boxes(id)%ix - 1) * nc + 1
       ihi = ilo + nc - 1

       call HYPRE_StructVectorGetBoxValues(mg%csolver%phi, ilo, ihi, tmp, ierr)
       tree%boxes(id)%cc(DTIMES(1:nc), mg%i_phi) = reshape(tmp, [DTIMES(nc)])
    end do
  end subroutine coarse_solver_get_phi

  !> Create a symmetric matrix on a grid
  subroutine hypre_create_matrix(A, grid, stencil_size, offsets, symmetric)
    type(c_ptr), intent(out) :: A
    type(c_ptr), intent(in)  :: grid
    integer, intent(in)      :: stencil_size
    integer, intent(in)      :: offsets(NDIM, stencil_size)
    integer, intent(in)      :: symmetric
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
    call HYPRE_StructMatrixSetSymmetric(A, symmetric, ierr)

    ! Indicate that the matrix coefficients are ready to be set */
    call HYPRE_StructMatrixInitialize(A, ierr)

    call HYPRE_StructStencilDestroy(stencil, ierr)

  end subroutine hypre_create_matrix

  ! Prepare to solve the system. The coefficient data in b and x is ignored
  ! here, but information about the layout of the data may be used.
  subroutine hypre_prepare_solve(cs)
    type(coarse_solve_t), intent(inout) :: cs
    integer                             :: ierr

    select case (cs%solver_type)
    case (coarse_solver_hypre_smg)
       call HYPRE_StructSMGCreate(MPI_COMM_WORLD, cs%solver, ierr)
       call HYPRE_StructSMGSetMaxIter(cs%solver, cs%max_iterations, ierr)
       call HYPRE_StructSMGSetTol(cs%solver, cs%tolerance, ierr)
       call HYPRE_StructSMGSetNumPreRelax(cs%solver, cs%n_cycle_down, ierr)
       call HYPRE_StructSMGSetNumPostRelax(cs%solver, cs%n_cycle_up, ierr)
       call HYPRE_StructSMGSetup(cs%solver, cs%matrix, cs%rhs, cs%phi, ierr)
    case (coarse_solver_hypre_pfmg)
       call HYPRE_StructPFMGCreate(MPI_COMM_WORLD, cs%solver, ierr)
       call HYPRE_StructPFMGSetMaxIter(cs%solver, cs%max_iterations, ierr)
       call HYPRE_StructPFMGSetTol(cs%solver, cs%tolerance, ierr)
       call HYPRE_StructPFMGSetNumPreRelax(cs%solver, cs%n_cycle_down, ierr)
       call HYPRE_StructPFMGSetNumPostRelax(cs%solver, cs%n_cycle_up, ierr)
       call HYPRE_StructPFMGSetup(cs%solver, cs%matrix, cs%rhs, cs%phi, ierr)
    case default
       error stop "hypre_prepare_solve: unknown solver type"
    end select
  end subroutine hypre_prepare_solve

  ! Solve the system A x = b
  subroutine coarse_solver(cs)
    type(coarse_solve_t), intent(in) :: cs
    integer                          :: ierr

    select case (cs%solver_type)
    case (coarse_solver_hypre_smg)
       call HYPRE_StructSMGSolve(cs%solver, cs%matrix, cs%rhs, cs%phi, ierr)
    case (coarse_solver_hypre_pfmg)
       call HYPRE_StructPFMGSolve(cs%solver, cs%matrix, cs%rhs, cs%phi, ierr)
    case default
       error stop "coarse_solver: unknown solver type"
    end select
  end subroutine coarse_solver

end module m_coarse_solver
