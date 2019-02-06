#include "../src/cpp_macros.h"
!> \example poisson_benchmark_Xd.f90
!>
!> This program can be used to benchmark the multigrid routines. For simplicity,
!> it does not compare results with a known solution.
program poisson_benchmark_Xd
  use m_af_all

  implicit none

  integer            :: i_phi
  integer            :: i_rhs
  integer            :: i_tmp

  type(af_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: mg_iter, n_args
  integer            :: n_cell, n_iterations, max_ref_lvl
  real(dp)           :: dr, time, runtime
  character(len=100) :: arg_string !, fname
  type(mg_t)         :: mg
  integer            :: count_rate,t_start, t_end

  print *, "Running poisson_benchmark_" // DIMNAME // ""
  print *, "Number of threads", af_get_max_threads()

  ! Get box size and mesh size from command line argument
  n_args = command_argument_count()

  if (n_args >= 1) then
     call get_command_argument(1, arg_string)
     read(arg_string, *) n_cell
  else
     print *, "No arguments specified, using default values"
     print *, "Usage: ./poisson_benchmark_" // DIMNAME // " n_cell max_ref_lvl runtime(s)"
     n_cell = 16
  end if

  if (n_args >= 2) then
     call get_command_argument(2, arg_string)
     read(arg_string, *) max_ref_lvl
  else
     max_ref_lvl = 2
  end if

  if (n_args >= 3) then
     call get_command_argument(3, arg_string)
     read(arg_string, *) runtime
     if (runtime < 1e-3_dp) stop "Run time should be > 1e-3 seconds"
  else
     runtime = 0.2_dp
  end if

  print *, "Box size:           ", n_cell
  print *, "Max refinement lvl: ", max_ref_lvl
  print *, "Run time (s):       ", runtime

  dr = 1.0_dp / n_cell

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       n_cell, &       ! A box contains box_size**DIM cells
       dr)             ! Distance between cells on base level

  call af_set_coarse_grid(tree, [DTIMES(n_cell)])

  call system_clock(t_start, count_rate)
  do
     ! For each box, set the initial conditions
     call af_loop_box(tree, set_init_cond)

     ! This updates the refinement of the tree, by at most one level per call.
     ! The second argument is a subroutine that is called for each box that can
     ! be refined or derefined, and it should set refinement flags. Information
     ! about the changes in refinement are returned in the third argument.
     call af_adjust_refinement(tree, ref_routine, ref_info)

     ! If no new boxes have been added, exit the loop
     if (ref_info%n_add == 0) exit
  end do
  call system_clock(t_end, count_rate)

  write(*,"(A,Es10.3,A)") " Wall-clock time generating AMR grid: ", &
       (t_end-t_start) / real(count_rate,dp), " seconds"

  call af_print_info(tree)

  ! Set the multigrid options.
  mg%i_phi        = i_phi       ! Solution variable
  mg%i_rhs        = i_rhs       ! Right-hand side variable
  mg%i_tmp        = i_tmp       ! Variable for temporary space
  mg%sides_bc     => af_bc_dirichlet_zero ! Method for boundary conditions

  ! Initialize the multigrid options. This performs some basics checks and sets
  ! default values where necessary.
  ! This routine does not initialize the multigrid variables i_phi, i_rhs
  ! and i_tmp. These variables will be initialized at the first call of mg_fas_fmg
  call mg_init(tree, mg)

  ! Warm-up call
  call mg_fas_fmg(tree, mg, .false., .false.)

  ! Test how long cycles take to determine the number of cycles
  n_iterations = 1000
  call system_clock(t_start, count_rate)
  do mg_iter = 1, n_iterations
     call mg_fas_fmg(tree, mg, .false., mg_iter>1)
     call system_clock(t_end, count_rate)
     time = (t_end-t_start) / real(count_rate, dp)
     if (time > 0.2_dp * runtime) exit
  end do

  n_iterations = ceiling((runtime/time) * mg_iter)

  ! Do the actual benchmarking
  call system_clock(t_start, count_rate)
  do mg_iter = 1, n_iterations
     ! Perform a FAS-FMG cycle (full approximation scheme, full multigrid). The
     ! third argument controls whether the residual is stored in i_tmp. The
     ! fourth argument controls whether to improve the current solution.
     call mg_fas_fmg(tree, mg, .false., mg_iter>1)

     ! This uses a V-cycle instead of an FMG-cycle
     ! call mg_fas_vcycle(tree, mg, .false.)

     ! If uncommented, this writes Silo output files containing the
     ! cell-centered values of the leaves of the tree
     ! write(fname, "(A,I0)") "poisson_benchmark_" // DIMNAME // "_", mg_iter
     ! call af_write_silo(tree, trim(fname), dir="output")
  end do
  call system_clock(t_end, count_rate)

  time = (t_end-t_start) / real(count_rate, dp)
  write(*, "(A,I0,A,E10.3,A)") &
       " Wall-clock time after ", n_iterations, &
       " iterations: ", time, " seconds"
  write(*, "(A,E10.3,A)") " Per iteration: ", time/n_iterations, " seconds"

  ! This writes a Silo output file containing the cell-centered values of the
  ! leaves of the tree (the boxes not covered by refinement).
  ! fname = "poisson_benchmark_" // DIMNAME // ""
  ! call af_write_silo(tree, trim(fname), dir="output")

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call af_destroy(tree)

contains

  ! Set refinement flags for box
  subroutine ref_routine(box, cell_flags)
    type(box_t), intent(in) :: box ! A list of all boxes in the tree
    integer, intent(out) :: cell_flags(DTIMES(box%n_cell))

    ! Fully refine up to max_ref_lvl
    if (box%lvl < max_ref_lvl) then
       cell_flags = af_do_ref
    else
       cell_flags = af_keep_ref
    end if
  end subroutine ref_routine

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box_t), intent(inout) :: box
    integer                     :: nc

    nc = box%n_cell
    box%cc(DTIMES(1:nc), i_rhs) = 1.0_dp
  end subroutine set_init_cond

end program poisson_benchmark_Xd
