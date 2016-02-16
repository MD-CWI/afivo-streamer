!> \example test_mg_benchmark_2d.f90

! This program can be used to benchmark the multigrid routines. For simplicity,
! it does not compare results with known solution.
program test_mg2_2d
  use m_a2_t
  use m_a2_core
  use m_a2_mg
  use m_a2_utils
  use m_a2_gc
  use m_a2_io

  implicit none

  integer, parameter :: n_boxes_base = 1
  integer, parameter :: n_var_cell = 3
  integer, parameter :: i_phi = 1
  integer, parameter :: i_rhs = 2
  integer, parameter :: i_tmp = 3

  type(a2_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: i, n_args, n_cell, max_ref_lvl
  integer            :: n_iterations, t_start, t_end
  integer            :: ix_list(2, n_boxes_base)
  integer            :: nb_list(4, n_boxes_base)
  real(dp)           :: dr, count_rate
  character(len=40)  :: arg_string
  type(mg2_t)        :: mg

  ! Get box size and mesh size from command line argument
  n_args = command_argument_count()

  if (n_args >= 1) then
     call get_command_argument(1, arg_string)
     read(arg_string, *) n_cell
  else
     print *, "No arguments specified, using default values"
     print *, "Usage: ./poisson_benchmark_2d n_cell max_ref_lvl n_iterations"
     n_cell = 16
  end if

  if (n_args >= 2) then
     call get_command_argument(2, arg_string)
     read(arg_string, *) max_ref_lvl
  else
     max_ref_lvl = 4
  end if

  if (n_args >= 3) then
     call get_command_argument(3, arg_string)
     read(arg_string, *) n_iterations
  else
     n_iterations = 100
  end if

  print *, "Box size:           ", n_cell
  print *, "Max refinement lvl: ", max_ref_lvl
  print *, "Num iterations:     ", n_iterations
  print *, ""

  dr = 1.0_dp / n_cell

  ! Initialize tree
  call a2_init(tree, & ! Tree to initialize
       n_cell, &       ! A box contains n_cell**DIM cells
       n_var_cell, &   ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       coarsen_to=2, & ! Add coarsened levels for multigrid
       cc_names=["phi", "rhs", "tmp"]) ! Variable names

  ! Set up geometry. These indices are used to define the coordinates of a box,
  ! by default the box at [1,1] touches the origin (x,y) = (0,0)
  ix_list(:, 1) = [1,1]         ! Set index of boxnn
  nb_list(:, 1) = -1            ! Dirichlet zero -> -1

  ! Create the base mesh, using the box indices and their neighbor information
  call a2_set_base(tree, ix_list, nb_list)

  do
     ! For each box, set the initial conditions
     call a2_loop_box(tree, set_init_cond)

     ! This updates the refinement of the tree, by at most one level per call.
     call a2_adjust_refinement(tree, ref_routine, ref_info)

     ! If no new boxes have been added, exit the loop
     if (ref_info%n_add == 0) exit
  end do

  call a2_print_info(tree)

  ! Set the multigrid options.
  mg%i_phi        = i_phi       ! Solution variable
  mg%i_rhs        = i_rhs       ! Right-hand side variable
  mg%i_tmp        = i_tmp       ! Variable for temporary space
  mg%sides_bc     => a2_bc_dirichlet_zero ! Method for boundary conditions

  ! Initialize the multigrid options. This performs some basics checks and sets
  ! default values where necessary.
  call mg2_init_mg(mg)

  ! Do the actual benchmarking
  call system_clock(t_start, count_rate)
  do i = 1, n_iterations
     ! Perform a FAS-FMG cycle (full approximation scheme, full multigrid). The
     ! third argument controls whether the residual is stored in i_tmp. The
     ! fourth argument controls whether to improve the current solution.
     call mg2_fas_fmg(tree, mg, .true., i>1)
  end do
  call system_clock(t_end, count_rate)

  write(*, "(A,E12.4)") " Wall-clock time (s): ", (t_end-t_start) / count_rate

  ! This writes a Silo output file containing the cell-centered values of the
  ! leaves of the tree (the boxes not covered by refinement).
  call a2_write_silo(tree, "poisson_benchmark_2d", dir="output")

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call a2_destroy(tree)

contains

  ! Return the refinement flag for boxes(id)
  subroutine ref_routine(boxes, id, ref_flag)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flag

    ! Fully refine up to max_ref_lvl
    if (boxes(id)%lvl < max_ref_lvl) ref_flag = a5_do_ref
  end subroutine ref_routine

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: nc

    nc = box%n_cell
    box%cc(1:nc, 1:nc, i_rhs) = 1.0_dp
  end subroutine set_init_cond

end program test_mg2_2d
