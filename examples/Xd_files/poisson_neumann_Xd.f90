#include "../src/cpp_macros_$Dd.h"
!> \example poisson_neumann_$Dd.f90
!>
!> Example showing how to use multigrid in combination with Neumann boundary
!> conditions.
program poisson_neumann_$Dd
  use m_a$D_all
  use m_gaussians

  implicit none

  integer, parameter :: box_size = 8
  integer, parameter :: n_boxes_base = 1
  integer, parameter :: n_iterations = 10
  integer, parameter :: n_var_cell = 4
  integer, parameter :: i_phi = 1
  integer, parameter :: i_rhs = 2
  integer, parameter :: i_tmp = 3
  integer, parameter :: i_err = 4

  type(a$D_t)        :: tree
  type(ref_info_t)   :: refine_info
  integer            :: mg_iter
  integer            :: ix_list($D, n_boxes_base)
  real(dp)           :: dr
  character(len=100) :: fname
  type(mg$D_t)       :: mg
  integer            :: count_rate,t_start,t_end

  print *, "Running poisson_neumann_$Dd"
  print *, "Number of threads", af_get_max_threads()

  ! The cell spacing at the coarsest grid level
  dr = 1.0_dp / box_size

  ! Initialize tree
  call a$D_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       n_var_cell, &   ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       coarsen_to=2, & ! Add coarsened levels for multigrid
       cc_names=["phi", "rhs", "tmp", "err"]) ! Variable names

  ! Set up geometry. These indices are used to define the coordinates of a box,
  ! by default the box at [1,1] touches the origin (x,y) = (0,0)
  ix_list(:, 1) = [DTIMES(1)]         ! Set index of box 1

  ! Create the base mesh, using the box indices and their neighbor information
  call a$D_set_base(tree, 1, ix_list)

  call system_clock(t_start, count_rate)
  do
     ! For each box, set the initial conditions
     call a$D_loop_box(tree, set_initial_condition)

     ! This updates the refinement of the tree, by at most one level per call.
     ! The second argument is a subroutine that is called for each box that can
     ! be refined or derefined, and it should set refinement flags. Information
     ! about the changes in refinement are returned in the third argument.
     call a$D_adjust_refinement(tree, refine_routine, refine_info)

     ! If no new boxes have been added, exit the loop
     if (refine_info%n_add == 0) exit
  end do
  call system_clock(t_end, count_rate)

  write(*,"(A,Es10.3,A)") " Wall-clock time generating AMR grid: ", &
       (t_end-t_start) / real(count_rate,dp), " seconds"

  call a$D_print_info(tree)

  ! Initialize the multigrid options. This performs some basics checks and sets
  ! default values where necessary.
  mg%i_phi        = i_phi       ! Solution variable
  mg%i_rhs        = i_rhs       ! Right-hand side variable
  mg%i_tmp        = i_tmp       ! Variable for temporary space
  mg%sides_bc     => sides_bc   ! Method for boundary conditions

  ! This routine does not initialize the multigrid fields boxes%i_phi,
  ! boxes%i_rhs and boxes%i_tmp. These fileds will be initialized at the
  ! first call of mg$D_fas_fmg
  call mg$D_init_mg(mg)

  print *, "Multigrid iteration | max residual | max error"
  call system_clock(t_start, count_rate)

  do mg_iter = 1, n_iterations
     call mg$D_fas_fmg(tree, mg, set_residual=.true., have_guess=(mg_iter>1))
     ! call mg$D_fas_vcycle(tree, mg, tree%highest_lvl, set_residual=.true.)

     call a$D_loop_box(tree, set_error)

     write(fname, "(A,I0)") "poisson_neumann_$Dd_", mg_iter
     call a$D_write_silo(tree, trim(fname), dir="output")
  end do
  call system_clock(t_end, count_rate)

  write(*, "(A,I0,A,E10.3,A)") &
       " Wall-clock time after ", n_iterations, &
       " iterations: ", (t_end-t_start) / real(count_rate, dp), &
       " seconds"

contains

  ! Return the refinement flags for box
  subroutine refine_routine(box, cell_flags)
    type(box$D_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(DTIMES(box%n_cell))

    if (box%lvl < 3 .and. all(box%r_min < 0.5_dp)) then
       cell_flags(DTIMES(:)) = af_do_ref
    else
       cell_flags(DTIMES(:)) = af_keep_ref
    end if
  end subroutine refine_routine

  ! This routine sets the initial conditions for each box
  subroutine set_initial_condition(box)
    type(box$D_t), intent(inout) :: box

    box%cc(DTIMES(:), i_rhs) = 0.0_dp
  end subroutine set_initial_condition

  ! Set the error compared to the analytic solution
  subroutine set_error(box)
    type(box$D_t), intent(inout) :: box
    integer                     :: IJK, nc
    real(dp)                    :: rr($D)

    nc = box%n_cell
    do KJI_DO(1,nc)
       rr = a$D_r_cc(box, [IJK])
       box%cc(IJK, i_err) = box%cc(IJK, i_phi) - rr(1)
    end do; CLOSE_DO
  end subroutine set_error

  ! This routine sets boundary conditions for a box, by filling its ghost cells
  ! with approriate values.
  subroutine sides_bc(box, nb, iv, bc_type)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)          :: nb      ! Direction for the boundary condition
    integer, intent(in)          :: iv      ! Index of variable
    integer, intent(out)         :: bc_type ! Type of boundary condition
    integer                      :: nc

    nc = box%n_cell

    ! Below the solution is specified in the approriate ghost cells
    select case (nb)
    case (a$D_neighb_lowx)             ! Lower-x direction
       call a$D_bc_dirichlet_zero(box, nb, iv, bc_type)
    case (a$D_neighb_highx)             ! Higher-x direction
#if $D == 2
       bc_type = af_bc_dirichlet
       box%cc(nc+1, 1:nc, iv) = 1
#elif $D == 3
       bc_type = af_bc_dirichlet
       box%cc(nc+1, 1:nc, 1:nc, iv) = 1
#endif
    case default
       call a$D_bc_neumann_zero(box, nb, iv, bc_type)
    end select
  end subroutine sides_bc

end program
