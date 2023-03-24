#include "../src/cpp_macros.h"
!> \example time_integration.f90
!>
!> An example to test the accuracy of time integration, solving a simple
!> equation: d/dt y = -y
program time_integration
  use m_af_all

  implicit none

  integer             :: i_phi, i_err
  integer             :: n, n_steps, integrator
  real(dp), parameter :: initial_value = exp(1.0_dp)
  real(dp)            :: dt, time, end_time
  real(dp)            :: max_error
  character(len=100)  :: argument
  type(af_t)          :: tree

  end_time = 1.0_dp
  dt       = 0.1_dp

  ! Read dt from command line if given
  if (command_argument_count() > 0) then
     call get_command_argument(1, argument)
     read(argument, *) dt
  end if

  n_steps  = nint(end_time/dt)
  dt       = end_time / n_steps

  write(*, "(A20,2A12)") "integrator          ", "dt", "max_error"

  do integrator = 1, af_num_integrators
     time = 0.0_dp
     call integrate(tree, integrator, time, dt, n_steps)
     call af_loop_box_arg(tree, set_error, [time])
     call af_tree_maxabs_cc(tree, i_err, max_error)
     write(*, "(A20,2E12.4)") af_integrator_names(integrator), dt, max_error
     call af_destroy(tree)
  end do

contains

  subroutine set_initial_condition(box)
    type(box_t), intent(inout) :: box
    box%cc(DTIMES(:), i_phi) = initial_value
  end subroutine set_initial_condition

  elemental real(dp) function solution(t)
    real(dp), intent(in) :: t
    solution = initial_value * exp(-t)
  end function solution

  subroutine set_error(box, time)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: time(:)
    integer                    :: IJK, nc

    nc = box%n_cell
    do KJI_DO(1,nc)
       box%cc(IJK, i_err) = box%cc(IJK, i_phi) - solution(time(1))
    end do; CLOSE_DO
  end subroutine set_error

  subroutine forward_euler(tree, dt, dt_stiff, dt_lim, time, s_deriv, n_prev, s_prev, &
       w_prev, s_out, i_step, n_steps)
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt             !< Time step
    real(dp), intent(in)      :: dt_stiff       !< Time step for stiff terms
    real(dp), intent(inout)   :: dt_lim         !< Computed time step limit
    real(dp), intent(in)      :: time           !< Current time
    integer, intent(in)       :: s_deriv        !< State to compute derivatives from
    integer, intent(in)       :: n_prev         !< Number of previous states
    integer, intent(in)       :: s_prev(n_prev) !< Previous states
    real(dp), intent(in)      :: w_prev(n_prev) !< Weights of previous states
    integer, intent(in)       :: s_out          !< Output state
    integer, intent(in)       :: i_step         !< Step of the integrator
    integer, intent(in)       :: n_steps        !< Total number of steps
    integer                   :: lvl, IJK, n, id, nc

    nc = tree%n_cell

    !$omp parallel private(lvl, n, id, IJK)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do n = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(n)
          associate(cc => tree%boxes(id)%cc)
            do KJI_DO(1,nc)
               cc(IJK, i_phi+s_out) = sum(w_prev * cc(IJK, i_phi+s_prev)) - &
                    dt_stiff * cc(IJK, i_phi+s_deriv)
            end do; CLOSE_DO
          end associate
       end do
       !$omp end do
    end do
    !$omp end parallel

    call af_gc_tree(tree, [i_phi])
  end subroutine forward_euler

  subroutine integrate(tree, integrator, time, dt, n_steps)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: integrator
    real(dp), intent(inout)   :: time
    real(dp), intent(in)      :: dt
    integer, intent(in)       :: n_steps
    real(dp)                  :: dt_lim
    integer                   :: n_copies
    integer, parameter        :: box_size = 8

    n_copies = af_advance_num_steps(integrator)
    call af_add_cc_variable(tree, "phi", ix=i_phi, n_copies=n_copies)
    call af_add_cc_variable(tree, "err", ix=i_err)
    call af_set_cc_methods(tree, i_phi, af_bc_neumann_zero)
    call af_init(tree, box_size, [DTIMES(1.0_dp)], [DTIMES(box_size)], &
         periodic=[DTIMES(.true.)], mem_limit_gb=1e-3_dp)

    call af_loop_box(tree, set_initial_condition)

    do n = 1, n_steps
       call af_advance(tree, dt, dt_lim, time, [i_phi], integrator, &
            forward_euler, implicit_solver)
    end do
  end subroutine integrate

  !> Implicit solver for equation d/dt y = -y, which has the solution
  !> y_n+1 = y_n / (1 + dt)
  subroutine implicit_solver(tree, dt_stiff, time, n_prev, s_prev, w_prev, s_out)
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt_stiff       !< Time step for stiff terms
    real(dp), intent(in)      :: time           !< Current time
    integer, intent(in)       :: n_prev         !< Number of previous states
    integer, intent(in)       :: s_prev(n_prev) !< Previous states
    real(dp), intent(in)      :: w_prev(n_prev) !< Weights of previous states
    integer, intent(in)       :: s_out          !< Output state

    integer                   :: lvl, IJK, n, id, nc

    nc = tree%n_cell

    !$omp parallel private(lvl, n, id, IJK)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do n = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(n)
          associate(cc => tree%boxes(id)%cc)
            do KJI_DO(1,nc)
               cc(IJK, i_phi+s_out) = sum(w_prev * cc(IJK, i_phi+s_prev)) / &
                    (1 + dt_stiff)
            end do; CLOSE_DO
          end associate
       end do
       !$omp end do
    end do
    !$omp end parallel

    call af_gc_tree(tree, [i_phi])
  end subroutine implicit_solver

end program time_integration
