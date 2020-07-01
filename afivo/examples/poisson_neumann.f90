#include "../src/cpp_macros.h"
!> \example poisson_neumann.f90
!>
!> Example showing how to use multigrid in combination with Neumann boundary
!> conditions.
program poisson_neumann
  use m_af_all
  use m_gaussians

  implicit none

  integer, parameter :: box_size     = 8
  integer, parameter :: n_iterations = 10
  integer            :: i_phi
  integer            :: i_rhs
  integer            :: i_tmp
  integer            :: i_err

  type(af_t)        :: tree
  type(ref_info_t)   :: refine_info
  integer            :: mg_iter
  character(len=100) :: fname
  type(mg_t)       :: mg
  integer            :: count_rate,t_start,t_end

  print *, "Running poisson_neumann_" // DIMNAME // ""
  print *, "Number of threads", af_get_max_threads()

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)
  call af_add_cc_variable(tree, "err", ix=i_err)

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [DTIMES(1.0_dp)], &
       [DTIMES(box_size)])

  call system_clock(t_start, count_rate)
  do
     ! For each box, set the initial conditions
     call af_loop_box(tree, set_initial_condition)

     ! This updates the refinement of the tree, by at most one level per call.
     ! The second argument is a subroutine that is called for each box that can
     ! be refined or derefined, and it should set refinement flags. Information
     ! about the changes in refinement are returned in the third argument.
     call af_adjust_refinement(tree, refine_routine, refine_info, 0)

     ! If no new boxes have been added, exit the loop
     if (refine_info%n_add == 0) exit
  end do
  call system_clock(t_end, count_rate)

  write(*,"(A,Es10.3,A)") " Wall-clock time generating AMR grid: ", &
       (t_end-t_start) / real(count_rate,dp), " seconds"

  call af_print_info(tree)

  ! Initialize the multigrid options. This performs some basics checks and sets
  ! default values where necessary.
  mg%i_phi        = i_phi       ! Solution variable
  mg%i_rhs        = i_rhs       ! Right-hand side variable
  mg%i_tmp        = i_tmp       ! Variable for temporary space
  mg%sides_bc     => sides_bc   ! Method for boundary conditions

  ! This routine does not initialize the multigrid fields boxes%i_phi,
  ! boxes%i_rhs and boxes%i_tmp. These fileds will be initialized at the
  ! first call of mg_fas_fmg
  call mg_init(tree, mg)

  print *, "Multigrid iteration | max residual | max error"
  call system_clock(t_start, count_rate)

  do mg_iter = 1, n_iterations
     call mg_fas_fmg(tree, mg, set_residual=.true., have_guess=(mg_iter>1))
     ! call mg_fas_vcycle(tree, mg, tree%highest_lvl, set_residual=.true.)

     call af_loop_box(tree, set_error)

     write(fname, "(A,I0)") "poisson_neumann_" // DIMNAME // "_", mg_iter
     call af_write_silo(tree, trim(fname), dir="output")
  end do
  call system_clock(t_end, count_rate)

  write(*, "(A,I0,A,E10.3,A)") &
       " Wall-clock time after ", n_iterations, &
       " iterations: ", (t_end-t_start) / real(count_rate, dp), &
       " seconds"

contains

  ! Return the refinement flags for box
  subroutine refine_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(DTIMES(box%n_cell))

    ! Refine around one corner
    if (box%lvl <= 4 .and. all(box%r_min < 0.25_dp)) then
       cell_flags(DTIMES(:)) = af_do_ref
    else
       cell_flags(DTIMES(:)) = af_keep_ref
    end if
  end subroutine refine_routine

  ! This routine sets the initial conditions for each box
  subroutine set_initial_condition(box)
    type(box_t), intent(inout) :: box

    box%cc(DTIMES(:), i_rhs) = 0.0_dp
  end subroutine set_initial_condition

  ! Set the error compared to the analytic solution
  subroutine set_error(box)
    type(box_t), intent(inout) :: box
    integer                     :: IJK, nc
    real(dp)                    :: rr(NDIM)

    nc = box%n_cell
    do KJI_DO(1,nc)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, i_err) = box%cc(IJK, i_phi) - rr(1)
    end do; CLOSE_DO
  end subroutine set_error

  ! This routine sets boundary conditions for a box
  subroutine sides_bc(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type
    integer                 :: nc

    nc = box%n_cell

    ! Below the solution is specified in the approriate ghost cells
    select case (nb)
    case (af_neighb_lowx)             ! Lower-x direction
       bc_type = af_bc_dirichlet
       bc_val = 0.0_dp
    case (af_neighb_highx)             ! Higher-x direction
       bc_type = af_bc_dirichlet
       bc_val = 1.0_dp
    case default
       bc_type = af_bc_neumann
       bc_val = 0.0_dp
    end select
  end subroutine sides_bc

end program
