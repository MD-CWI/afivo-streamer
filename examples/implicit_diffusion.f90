#include "../src/cpp_macros.h"
!> \example implicit_diffusion_Xd.f90
!>
!> An implicit diffusion example, showing how the multigrid methods can be used
!> to solve the diffusion equation with a backward Euler scheme.
program implicit_diffusion_Xd
  use m_af_all

  implicit none

  integer, parameter :: box_size = 8
  integer            :: i_phi
  integer            :: i_rhs
  integer            :: i_tmp
  integer            :: i_err

  real(dp), parameter :: domain_len = 2 * acos(-1.0_dp)
  real(dp), parameter :: dr = domain_len / box_size
  real(dp), parameter :: diffusion_coeff = 1.0_dp
#if NDIM == 2
  real(dp), parameter :: solution_modes(NDIM) = [1, 1]
#elif NDIM == 3
  real(dp), parameter :: solution_modes(NDIM) = [1, 1, 0]
#endif

  type(af_t)         :: tree
  type(mg_t)         :: mg
  type(ref_info_t)   :: refine_info
  integer            :: time_steps, output_cnt
  real(dp)           :: time, end_time
  real(dp), parameter :: dt = 0.1_dp
  character(len=100) :: fname

  print *, "Running implicit_diffusion_" // DIMNAME // ""
  print *, "Number of threads", af_get_max_threads()

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)
  call af_add_cc_variable(tree, "err", ix=i_err)

  ! Initialize tree
  call af_init(tree, box_size, dr=dr)

  call af_set_coarse_grid(tree, [DTIMES(box_size)], [DTIMES(.true.)])
  call af_print_info(tree)

  ! Set up the initial conditions
  do

     ! For each box, set the initial conditions
     call af_loop_box(tree, set_initial_condition)

     ! Fill ghost cells for variables i_phi on the sides of all boxes, using
     ! af_gc_interp on refinement boundaries: Interpolation between fine
     ! points and coarse neighbors to fill ghost cells near refinement
     ! boundaries.
     ! Fill ghost cells near physical boundaries using Dirichlet zero

     call af_gc_tree(tree, i_phi, af_gc_interp, af_bc_dirichlet_zero)

     ! Adjust the refinement of a tree using refine_routine (see below) for grid
     ! refinement.
     ! Routine af_adjust_refinement sets the bit af_bit_new_children for each box
     ! that is refined.  On input, the tree should be balanced. On output,
     ! the tree is still balanced, and its refinement is updated (with at most
     ! one level per call).
     call af_adjust_refinement(tree, &           ! tree
                               refine_routine, & ! Refinement function
                               refine_info)      ! Information about refinement

     ! If no new boxes have been added, exit the loop
     if (refine_info%n_add == 0) exit
  end do

  mg%i_phi    = i_phi                 ! Solution variable
  mg%i_rhs    = i_rhs                 ! Right-hand side variable
  mg%i_tmp    = i_tmp                 ! Variable for temporary space
  mg%sides_bc => af_bc_dirichlet_zero ! Method for boundary conditions

  ! The methods defined below implement a backward Euler method for the heat
  ! equation, by changing the elliptic operator for the multigrid procedure.
  mg%box_op   => box_op_diff
  mg%box_gsrb => box_gsrb_diff
  mg%box_stencil => box_op_stencil

  ! This routine does not initialize the multigrid fields boxes%i_phi,
  ! boxes%i_rhs and boxes%i_tmp. These fields will be initialized at the
  ! first call of mg_fas_fmg
  call mg_init(tree, mg)

  output_cnt = 0
  time       = 0
  end_time   = 2.0_dp
  time_steps = 0
  time       = 0

  ! Starting simulation
  do
     time_steps = time_steps + 1

     output_cnt = output_cnt + 1
     write(fname, "(A,I0)") "implicit_diffusion_" // DIMNAME // "_", output_cnt

     ! Call procedure set_error (see below) for each box in tree, with argument time
     call af_loop_box_arg(tree, set_error, [time])

     ! Write the cell centered data of tree to a vtk unstructured file fname.
     ! Only the leaves of the tree are used
     call af_write_silo(tree, trim(fname), output_cnt, time, dir="output")

     if (time > end_time) exit

     ! Call set_rhs (see below) for each box in tree
     call af_loop_box(tree, set_rhs)

     ! Perform a FAS-FMG cycle (full approximation scheme, full multigrid).
     call mg_fas_fmg(tree, &     ! Tree to do multigrid on
          mg, &                    ! Multigrid options
          set_residual = .true., & ! If true, store residual in i_tmp
          have_guess   = .true.)  ! If false, start from phi = 0

     time = time + dt
  end do

contains

  !> Return the refinement flag for box
  subroutine refine_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(DTIMES(box%n_cell))

    if (box%dr > 2e-2_dp * domain_len) then
       cell_flags = af_do_ref
    else
       cell_flags = af_keep_ref
    end if
  end subroutine refine_routine

  !> This routine sets the initial conditions for each box
  subroutine set_initial_condition(box)
    type(box_t), intent(inout) :: box
    integer                     :: IJK, nc
    real(dp)                    :: rr(NDIM)

    nc = box%n_cell

    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, i_phi) = solution(rr, 0.0_dp)
    end do; CLOSE_DO
  end subroutine set_initial_condition

  !> This routine computes the error in i_phi
  subroutine set_error(box, time)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)        :: time(:)
    integer                     :: IJK, nc
    real(dp)                    :: rr(NDIM)

    nc = box%n_cell
    do KJI_DO(1,nc)
          rr = af_r_cc(box, [IJK])
          box%cc(IJK, i_err) = &
               box%cc(IJK, i_phi) - solution(rr, time(1))
    end do; CLOSE_DO
  end subroutine set_error

  function solution(rr, t) result(sol)
    real(dp), intent(in) :: rr(NDIM), t
    real(dp)             :: sol, tmp(NDIM)

    tmp = solution_modes * rr
    sol = 1 + product(cos(tmp)) * &
         exp(-sum(solution_modes**2) * diffusion_coeff * t)
  end function solution

  !> This routine computes the right hand side per box
  subroutine set_rhs(box)
    type(box_t), intent(inout) :: box
    integer                     :: nc

    nc = box%n_cell
    box%cc(DTIMES(1:nc), i_rhs) = box%cc(DTIMES(1:nc), i_phi)
  end subroutine set_rhs

  ! Compute L * phi, where L corresponds to (D * dt * nabla^2 - 1)
  subroutine box_op_diff(box, i_out, mg)
    type(box_t), intent(inout) :: box    !< Box to operate on
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    type(mg_t), intent(in)     :: mg
    integer                     :: IJK, nc, i_phi
    real(dp)                    :: tmp

    nc    = box%n_cell
    tmp   = diffusion_coeff * dt / box%dr**2
    i_phi = mg%i_phi

    do KJI_DO(1,nc)
#if NDIM == 2
       box%cc(i, j, i_out) = tmp * ((4 + 1/tmp) * &
            box%cc(i, j, i_phi) &
            - box%cc(i-1, j, i_phi) &
            - box%cc(i+1, j, i_phi) &
            - box%cc(i, j-1, i_phi) &
            - box%cc(i, j+1, i_phi))
#elif NDIM == 3
       box%cc(i, j, k, i_out) = tmp * ((6 + 1/tmp) * &
            box%cc(i, j, k, i_phi) &
            - box%cc(i-1, j, k, i_phi) &
            - box%cc(i+1, j, k, i_phi) &
            - box%cc(i, j-1, k, i_phi) &
            - box%cc(i, j+1, k, i_phi) &
            - box%cc(i, j, k-1, i_phi) &
            - box%cc(i, j, k+1, i_phi))
#endif
    end do; CLOSE_DO
  end subroutine box_op_diff

  subroutine box_op_stencil(box, mg, stencil, bc_to_rhs)
    type(box_t), intent(in) :: box
    type(mg_t), intent(in)  :: mg
    real(dp), intent(inout) :: stencil(2*NDIM+1, DTIMES(box%n_cell))
    real(dp), intent(inout) :: bc_to_rhs(box%n_cell**(NDIM-1), af_num_neighbors)
    real(dp)                :: tmp

    tmp = diffusion_coeff * dt / box%dr**2

    stencil(1, DTIMES(:)) = tmp * (2.0_dp * NDIM + 1/tmp)
    stencil(2:, DTIMES(:)) = -tmp

    call mg_stencil_handle_boundaries(box, mg, stencil, bc_to_rhs)
  end subroutine box_op_stencil

  ! Locally solve L * phi = rhs, where L corresponds to (D * dt * nabla^2 - 1)
  subroutine box_gsrb_diff(box, redblack_cntr, mg)
    type(box_t), intent(inout) :: box            !< Box to operate on
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    type(mg_t), intent(in)     :: mg
    integer                     :: IJK, i0, nc, i_phi, i_rhs
    real(dp)                    :: tmp

    tmp   = diffusion_coeff * dt / box%dr**2
    nc    = box%n_cell
    i_phi = mg%i_phi
    i_rhs = mg%i_rhs

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
#if NDIM == 2
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          box%cc(i, j, i_phi) = 1 / (4 + 1/tmp) * ( &
               box%cc(i+1, j, i_phi) + box%cc(i-1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi) + &
               1/tmp * box%cc(i, j, i_rhs))
       end do
    end do
#elif NDIM == 3
    do k = 1, nc
       do j = 1, nc
          i0 = 2 - iand(ieor(redblack_cntr, k+j), 1)
          do i = i0, nc, 2
             box%cc(i, j, k, i_phi) = 1 / (6 + 1/tmp) * ( &
                  box%cc(i+1, j, k, i_phi) + box%cc(i-1, j, k, i_phi) + &
                  box%cc(i, j+1, k, i_phi) + box%cc(i, j-1, k, i_phi) + &
                  box%cc(i, j, k+1, i_phi) + box%cc(i, j, k-1, i_phi) + &
                  1/tmp * box%cc(i, j, k, i_rhs))
          end do
       end do
    end do
#endif
  end subroutine box_gsrb_diff

end program implicit_diffusion_Xd
