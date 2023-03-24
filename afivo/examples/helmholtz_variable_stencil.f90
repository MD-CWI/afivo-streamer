#include "../src/cpp_macros.h"
!> \example helmholtz_variable_stencil.f90
!>
!> Example of solving a Helmholtz equation with a stencil that varies in time.
!> Useful ingredient for the implicit solution of diffusion equations.
program helmholtz_variable_stencil
  use m_af_all

  implicit none

  integer, parameter :: box_size = 8
  integer            :: i_err

  real(dp), parameter :: domain_len = 2 * acos(-1.0_dp)
  real(dp), parameter :: diffusion_coeff = 1.0_dp
#if NDIM == 1
  real(dp), parameter :: solution_modes(NDIM) = [1]
#elif NDIM == 2
  real(dp), parameter :: solution_modes(NDIM) = [1, 1]
#elif NDIM == 3
  real(dp), parameter :: solution_modes(NDIM) = [1, 1, 0]
#endif

  type(af_t)         :: tree
  type(mg_t)         :: mg
  integer            :: i, n, k_factor, n_steps, n_steps_initial
  real(dp)           :: time, dt, dt_initial
  character(len=100) :: fname

  print *, "Running helmholtz_variable_stencil_" // DIMNAME // ""
  print *, "Number of threads", af_get_max_threads()

  dt_initial = 0.1_dp
  n_steps_initial = 10

  call af_add_cc_variable(tree, "phi", ix=mg%i_phi)
  call af_add_cc_variable(tree, "rhs", ix=mg%i_rhs)
  call af_add_cc_variable(tree, "tmp", ix=mg%i_tmp)
  call af_add_cc_variable(tree, "err", ix=i_err)

  call af_init(tree, box_size, [DTIMES(domain_len)], &
       [DTIMES(box_size)], periodic=[DTIMES(.true.)])

  mg%sides_bc => af_bc_dirichlet_zero ! Method for boundary conditions
  mg%helmholtz_lambda = 0.            ! Updated below
  mg%helmholtz_lambda = 1/(diffusion_coeff * dt_initial)

  call mg_init(tree, mg)
  call af_print_info(tree)
  call af_refine_up_to_lvl(tree, 3)

  i = 0

  do k_factor = 1, 10
     time                = 0
     dt                  = dt_initial / k_factor
     n_steps             = n_steps_initial * k_factor
     mg%helmholtz_lambda = 1/(diffusion_coeff * dt)
     print *, k_factor, dt, n_steps

     if (k_factor > 1) call mg_update_operator_stencil(tree, mg)

     call af_loop_box(tree, set_initial_condition)
     call af_gc_tree(tree, [mg%i_phi])

     do n = 1, n_steps
        call af_loop_box(tree, set_rhs)
        call mg_fas_fmg(tree, mg, .true., .true.)
        time = time + dt
        call af_loop_box_arg(tree, set_error, [time])
     end do

     write(fname, "(A,I0)") "output/helmholtz_variable_stencil_" // DIMNAME // "_", k_factor
     call af_write_silo(tree, trim(fname), k_factor, time)
  end do

contains

  !> This routine sets the initial conditions for each box
  subroutine set_initial_condition(box)
    type(box_t), intent(inout) :: box
    integer                     :: IJK, nc
    real(dp)                    :: rr(NDIM)

    nc = box%n_cell

    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, mg%i_phi) = solution(rr, 0.0_dp)
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
            box%cc(IJK, mg%i_phi) - solution(rr, time(1))
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
    integer                    :: nc

    nc = box%n_cell
    box%cc(DTIMES(1:nc), mg%i_rhs) = -1/(dt * diffusion_coeff) * &
         box%cc(DTIMES(1:nc), mg%i_phi)
  end subroutine set_rhs

end program helmholtz_variable_stencil
