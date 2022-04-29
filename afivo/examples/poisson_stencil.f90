#include "../src/cpp_macros.h"
!> \example poisson_stencil.f90
!>
!> Example showing how to use multigrid with a numerical stencil
program poisson_stencil
  use m_af_all
  use m_gaussians

  implicit none

  integer, parameter :: box_size          = 16
  integer, parameter :: n_iterations      = 10
  integer            :: domain_size(NDIM)
  real(dp)           :: domain_len(NDIM)
  integer            :: i_phi
  integer            :: i_sol
  integer            :: i_rhs
  integer            :: i_err
  integer            :: i_tmp
  type(af_t)         :: tree
  type(ref_info_t)   :: refine_info
  integer            :: mg_iter
  real(dp)           :: residu(2), anal_err(2)
  character(len=100) :: fname
  type(mg_t)         :: mg
  type(gauss_t)      :: gs
  integer            :: count_rate,t_start,t_end

  print *, "Running poisson_stencil_" // DIMNAME
  print *, "Number of threads", af_get_max_threads()

  call gauss_init(gs, [1.0_dp, 1.0_dp], [0.04_dp, 0.04_dp], &
       reshape([DTIMES(0.1_dp), &
       DTIMES(0.75_dp)], [NDIM,2]))

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "sol", ix=i_sol)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)
  call af_add_cc_variable(tree, "err", ix=i_err)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)

  domain_len(1)   = 3.0_dp
  domain_len(2:)  = 1.0_dp
  domain_size(1)  = 3 * box_size
  domain_size(2:) = box_size

  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       domain_len, &
       domain_size)

  call system_clock(t_start, count_rate)
  do
     call af_loop_box(tree, set_initial_condition)
     call af_adjust_refinement(tree, refine_routine, refine_info)
     if (refine_info%n_add == 0) exit
  end do
  call system_clock(t_end, count_rate)

  write(*,"(A,Es10.3,A)") " Wall-clock time generating AMR grid: ", &
       (t_end-t_start) / real(count_rate,dp), " seconds"

  call af_print_info(tree)

  mg%i_phi        = i_phi       ! Solution variable
  mg%i_rhs        = i_rhs       ! Right-hand side variable
  mg%i_tmp        = i_tmp       ! Variable for temporary space
  mg%sides_bc     => sides_bc   ! Method for boundary conditions

  mg%operator_key = mg_normal_box

  call mg_init(tree, mg)

  print *, "Multigrid iteration | max residual | max error"
  call system_clock(t_start, count_rate)

  do mg_iter = 1, n_iterations
     call mg_fas_fmg(tree, mg, set_residual=.true., have_guess=(mg_iter>1))

     call af_loop_box(tree, set_error)

     call af_tree_min_cc(tree, i_tmp, residu(1))
     call af_tree_max_cc(tree, i_tmp, residu(2))
     call af_tree_min_cc(tree, i_err, anal_err(1))
     call af_tree_max_cc(tree, i_err, anal_err(2))
     write(*,"(I8,13x,2(Es14.5))") mg_iter, maxval(abs(residu)), &
          maxval(abs(anal_err))

     write(fname, "(A,I0)") "output/poisson_stencil_" // DIMNAME // "_", mg_iter
     call af_write_silo(tree, trim(fname))
  end do
  call system_clock(t_end, count_rate)

  write(*, "(A,I0,A,E10.3,A)") &
       " Wall-clock time after ", n_iterations, &
       " iterations: ", (t_end-t_start) / real(count_rate, dp), &
       " seconds"

  call af_destroy(tree)

contains

  subroutine refine_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(DTIMES(box%n_cell))
    integer                  :: IJK, nc
    real(dp)                 :: rr(NDIM), dr2, drhs

    nc = box%n_cell
    dr2 = maxval(box%dr)**2

    do KJI_DO(1,nc)
       rr = af_r_cc(box, [IJK])

       drhs = dr2 * box%cc(IJK, i_rhs)

       if (abs(drhs) > 1e-3_dp .and. box%lvl < 7 - NDIM) then
          cell_flags(IJK) = af_do_ref
       else
          cell_flags(IJK) = af_keep_ref
       end if
    end do; CLOSE_DO
  end subroutine refine_routine

  subroutine set_initial_condition(box)
    type(box_t), intent(inout) :: box
    integer                     :: IJK, nc
    real(dp)                    :: rr(NDIM)

    nc = box%n_cell
    do KJI_DO(1,nc)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, i_rhs) = gauss_laplacian(gs, rr)
       box%cc(IJK, i_sol) = gauss_value(gs, rr)
    end do; CLOSE_DO
  end subroutine set_initial_condition

  subroutine set_error(box)
    type(box_t), intent(inout) :: box
    integer                     :: IJK, nc
    real(dp)                    :: rr(NDIM)

    nc = box%n_cell

    do KJI_DO(1,nc)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, i_err) = box%cc(IJK, i_phi) - gauss_value(gs, rr)
    end do; CLOSE_DO
  end subroutine set_error

  subroutine sides_bc(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type
    integer                 :: n

    bc_type = af_bc_dirichlet

    do n = 1, box%n_cell**(NDIM-1)
       bc_val(n) = gauss_value(gs, coords(:, n))
    end do
  end subroutine sides_bc

end program poisson_stencil
