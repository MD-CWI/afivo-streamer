!> \example poisson_benchmark_3d.f90

! This program can be used to benchmark the multigrid routines. For simplicity,
! it does not compare results with known solution.
program poisson_benchmark_3d
  use m_a3_types
  use m_a3_core
  use m_a3_multigrid
  use m_a3_utils
  use m_a3_ghostcell
  use m_a3_output

  implicit none

  integer, parameter :: n_boxes_base = 1
  integer, parameter :: n_var_cell = 3
  integer, parameter :: i_phi = 1
  integer, parameter :: i_rhs = 2
  integer, parameter :: i_tmp = 3

  type(a3_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: mg_iter, n_args
  integer            :: n_cell, n_iterations, max_ref_lvl
  integer            :: ix_list(3, n_boxes_base)
  integer            :: nb_list(6, n_boxes_base)
  real(dp)           :: dr, time
  character(len=100) :: fname, arg_string
  type(mg3_t)        :: mg
  integer            :: count_rate,t_start, t_end

  print *, "Running poisson_benchmark_3d"
  print *, "Number of threads", af_get_max_threads()

  ! Get box size and mesh size from command line argument
  n_args = command_argument_count()

  if (n_args >= 1) then
     call get_command_argument(1, arg_string)
     read(arg_string, *) n_cell
  else
     print *, "No arguments specified, using default values"
     print *, "Usage: ./poisson_benchmark_3d n_cell max_ref_lvl n_iterations"
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

  dr = 1.0_dp / n_cell

  ! Initialize tree
  call a3_init(tree, & ! Tree to initialize
       n_cell, &       ! A box contains box_size**DIM cells
       n_var_cell, &   ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       coarsen_to=2, & ! Add coarsened levels for multigrid
       cc_names=["phi", "rhs", "tmp"]) ! Variable names

  ! Set up geometry. These indices are used to define the coordinates of a box,
  ! by default the box at [1,1] touches the origin (x,y) = (0,0)
  ix_list(:, 1) = [1,1,1]       ! Set index of box 1

  ! Set neighbors for box one, negative values indicate a physical boundary
  nb_list(:, 1) = -1            ! Dirichlet zero -> -1

  ! Create the base mesh, using the box indices and their neighbor information
  call a3_set_base(tree, ix_list, nb_list)
  call a3_print_info(tree)

  call system_clock(t_start, count_rate)
  do
     ! For each box, set the initial conditions
     call a3_loop_box(tree, set_init_cond)

     ! This updates the refinement of the tree, by at most one level per call.
     ! The second argument is a subroutine that is called for each box that can
     ! be refined or derefined, and it should set refinement flags. Information
     ! about the changes in refinement are returned in the third argument.
     call a3_adjust_refinement(tree, ref_routine, ref_info)

     ! If no new boxes have been added, exit the loop
     if (ref_info%n_add == 0) exit
  end do
  call system_clock(t_end, count_rate)

  write(*,"(A,Es10.3,A)") " Wall-clock time generating AMR grid: ", &
       (t_end-t_start) / real(count_rate,dp), " seconds"

  call a3_print_info(tree)

  ! Set the multigrid options.
  mg%i_phi        = i_phi       ! Solution variable
  mg%i_rhs        = i_rhs       ! Right-hand side variable
  mg%i_tmp        = i_tmp       ! Variable for temporary space
  mg%sides_bc     => a3_bc_dirichlet_zero ! Method for boundary conditions

  ! Initialize the multigrid options. This performs some basics checks and sets
  ! default values where necessary.
  ! This routine does not initialize the multigrid variables i_phi, i_rhs
  ! and i_tmp. These variables will be initialized at the first call of mg3_fas_fmg
  call mg3_init_mg(mg)

  ! Do the actual benchmarking
  call system_clock(t_start, count_rate)
  do mg_iter = 1, n_iterations
     ! Perform a FAS-FMG cycle (full approximation scheme, full multigrid). The
     ! third argument controls whether the residual is stored in i_tmp. The
     ! fourth argument controls whether to improve the current solution.
     call mg3_fas_fmg(tree, mg, .true., mg_iter>1)

     ! If uncommented, this writes Silo output files containing the
     ! cell-centered values of the leaves of the tree
     ! write(fname, "(A,I0)") "poisson_benchmark_3d_", mg_iter
     ! call a3_write_silo(tree, trim(fname), dir="output")
  end do
  call system_clock(t_end, count_rate)

  time = (t_end-t_start) / real(count_rate, dp)
  write(*, "(A,I0,A,E10.3,A)") &
       " Wall-clock time after ", n_iterations, &
       " iterations: ", time, " seconds"
  write(*, "(A,E10.3,A)") " Per iteration: ", time/n_iterations, " seconds"

  ! This writes a Silo output file containing the cell-centered values of the
  ! leaves of the tree (the boxes not covered by refinement).
  fname = "poisson_benchmark_3d"
  call a3_write_silo(tree, trim(fname), dir="output")

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call a3_destroy(tree)

contains

  ! Set refinement flags for box
  subroutine ref_routine(box, cell_flags)
    type(box3_t), intent(in) :: box ! A list of all boxes in the tree
    integer, intent(out)     :: cell_flags(box%n_cell, box%n_cell, box%n_cell)

    ! Fully refine up to max_ref_lvl
    if (box%lvl < max_ref_lvl) then
       cell_flags = af_do_ref
    else
       cell_flags = af_keep_ref
    end if
  end subroutine ref_routine

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box3_t), intent(inout) :: box
    integer                     :: nc

    nc = box%n_cell
    box%cc(1:nc, 1:nc, 1:nc, i_rhs) = 1.0_dp
  end subroutine set_init_cond

end program poisson_benchmark_3d
