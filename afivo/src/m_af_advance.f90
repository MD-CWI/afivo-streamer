#include "cpp_macros.h"
!> Module with methods to perform time integration
module m_af_advance
  use m_af_types

  implicit none
  private

  integer, parameter         :: n_integrators      = 3
  integer, parameter, public :: af_forward_euler   = 1
  integer, parameter, public :: af_heuns_method    = 2
  integer, parameter, public :: af_midpoint_method = 3

  !> How many steps the time integrators take
  integer, parameter, public :: &
       af_advance_num_steps(n_integrators) = [1, 2, 2]

  !> How many variable copies are required for the time integrators
  integer, parameter :: req_copies(n_integrators) = af_advance_num_steps

  interface
     !> Interface for a generic forward Euler scheme for time integration
     !>
     !> This method should advance the solution over a time dt. The method
     !> assumes that several copies are stored for the variables to be
     !> integrated. It should then operate on these different copies, which
     !> correspond to temporal states. In this way, higher-order time
     !> integration schemes can be constructed.
     !>
     !> The meaning of the temporal states is as follows. For an equation y' =
     !> f(y), the method should perform: y_out = y_prev + dt * f(y_deriv).
     !>
     !> If the index of the variable `y` is `i`, then the index of `y_out` is
     !> `i+s_out`, the index of `y_prev` is `i+s_prev` etc.
     subroutine subr_feuler(tree, dt, dt_lim, time, s_deriv, s_prev, s_out, &
          i_step, n_steps)
       import
       type(af_t), intent(inout) :: tree
       real(dp), intent(in)      :: dt      !< Time step
       real(dp), intent(inout)   :: dt_lim  !< Computed time step limit
       real(dp), intent(in)      :: time    !< Current time
       integer, intent(in)       :: s_deriv !< State to compute derivatives from
       integer, intent(in)       :: s_prev  !< Previous state
       integer, intent(in)       :: s_out   !< Output state
       integer, intent(in)       :: i_step   !< Step of the integrator
       integer, intent(in)       :: n_steps   !< Total number of steps
     end subroutine subr_feuler
  end interface

  public :: af_advance

contains

  !> Advance solution over dt using time_integrator
  !>
  !> The user should supply a forward Euler method as documented in subr_feuler.
  !> The indices of the cell-centered variables that will be operated on should
  !> also be provided, so that higher-order schemes can be constructed
  !> automatically from the forward Euler method.
  subroutine af_advance(tree, dt, dt_lim, time, i_cc, time_integrator, forward_euler)
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt      !< Current time step
    real(dp), intent(out)     :: dt_lim  !< Time step limit
    real(dp), intent(inout)   :: time    !< Current time
    integer, intent(in)       :: i_cc(:) !< Index of cell-centered variables
    !> One of the pre-defined time integrators (e.g. af_heuns_method)
    integer, intent(in)       :: time_integrator
    !> Forward Euler method provided by the user
    procedure(subr_feuler)    :: forward_euler
    integer                   :: n_steps

    if (time_integrator < 1 .or. time_integrator > n_integrators) &
         error stop "Invalid time integrator"

    if (any(tree%cc_num_copies(i_cc) < req_copies(time_integrator))) &
         error stop "Not enough copies available"

    n_steps = af_advance_num_steps(time_integrator)
    dt_lim = 1e100_dp

    select case (time_integrator)
    case (af_forward_euler)
       call forward_euler(tree, dt, dt_lim, time, 0, 0, 0, 1, n_steps)
       time = time + dt
    case (af_midpoint_method)
       call forward_euler(tree, 0.5_dp * dt, dt_lim, time, 0, 0, 1, 1, n_steps)
       call forward_euler(tree, dt, dt_lim, time, 1, 0, 0, 2, n_steps)
       time = time + dt
    case (af_heuns_method)
       call forward_euler(tree, dt, dt_lim, time, 0, 0, 1, 1, n_steps)
       time = time + dt
       call forward_euler(tree, dt, dt_lim, time, 1, 1, 1, 2, n_steps)
       call combine_substeps(tree, i_cc, 2, [0, 1], [0.5_dp, 0.5_dp], 0)
    end select

  end subroutine af_advance

  subroutine combine_substeps(tree, ivs, n_in, s_in, coeffs, s_out)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: ivs(:)
    integer, intent(in)       :: n_in
    integer, intent(in)       :: s_in(n_in)
    real(dp), intent(in)      :: coeffs(n_in)
    integer, intent(in)       :: s_out
    integer                   :: lvl, i, id, n, nc
    real(dp), allocatable     :: tmp(DTIMES(:), :)

    nc = tree%n_cell

    !$omp parallel private(lvl, i, id, tmp, n)
    allocate(tmp(DTIMES(0:nc+1), size(ivs)))

    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)

          tmp = 0.0_dp
          do n = 1, n_in
             tmp = tmp + coeffs(n) * &
                  tree%boxes(id)%cc(DTIMES(:), ivs+s_in(n))
          end do
          tree%boxes(id)%cc(DTIMES(:), ivs+s_out) = tmp
       end do
       !$omp end do nowait
    end do
    !$omp end parallel
  end subroutine combine_substeps

end module m_af_advance
