#include "../src/cpp_macros.h"
!> \example ghostcell_benchmark_Xd.f90
!>
!> This program can be used to benchmark the ghostcells routines
program ghostcell_benchmark_Xd
  use m_af_all

  implicit none

  integer, parameter :: n_boxes_base = 1
  integer, parameter :: n_var_cell = 1
  integer, parameter :: i_phi = 1

  type(af_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: n_args
  integer            :: n_cell, it, n_iterations, max_ref_lvl
  integer            :: ix_list(NDIM, n_boxes_base)
  real(dp)           :: dr, time
  character(len=100) :: arg_string
  integer            :: count_rate, t_start, t_end

  print *, "Running ghostcell_benchmark_" // DIMNAME // ""
  print *, "Number of threads", af_get_max_threads()

  ! Get box size and mesh size from command line argument
  n_args = command_argument_count()

  if (n_args >= 1) then
     call get_command_argument(1, arg_string)
     read(arg_string, *) n_cell
  else
     print *, "No arguments specified, using default values"
     print *, "Usage: ./ghostcell_benchmark_" // DIMNAME // " n_cell max_ref_lvl n_iterations"
     print *, ""
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
  call af_init(tree, & ! Tree to initialize
       n_cell, &       ! A box contains box_size**DIM cells
       n_var_cell, &   ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       cc_names=["phi"]) ! Variable names

  ! Set up geometry. These indices are used to define the coordinates of a box,
  ! by default the box at [1,1] touches the origin (x,y) = (0,0)
  ix_list(:, 1) = [DTIMES(1)]       ! Set index of box 1

  ! Create the base mesh, using the box indices and their neighbor information
  call af_set_base(tree, 1, ix_list)

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

  ! Do the actual benchmarking
  call system_clock(t_start, count_rate)
  do it = 1, n_iterations
     call af_gc_tree(tree, i_phi, af_gc_interp, af_bc_dirichlet_zero, .false.)
  end do
  call system_clock(t_end, count_rate)

  time = (t_end-t_start) / real(count_rate, dp)
  write(*, "(A,I0,A,E10.3,A)") &
       " Wall-clock time after ", n_iterations, &
       " iterations: ", time, " seconds"
  write(*, "(A,E10.3,A)") " Per iteration: ", time/n_iterations, " seconds"

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
    box%cc(DTIMES(1:nc), i_phi) = 1.0_dp
  end subroutine set_init_cond

end program ghostcell_benchmark_Xd
