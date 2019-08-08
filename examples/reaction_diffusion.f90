#include "../src/cpp_macros.h"
!> \example reaction_diffusion.f90
program reaction_diffusion_Xd
  use m_af_all

  implicit none

  integer, parameter :: box_size = 8
  integer            :: i_u
  integer            :: i_v
  integer            :: i_phi1, i_phi2
  integer            :: i_rhs1, i_rhs2
  integer            :: i_tmp

  real(dp), parameter :: domain_len = 1.0_dp

  type(af_t)         :: tree
  type(mg_t)         :: mg1, mg2
  type(ref_info_t)   :: refine_info
  integer            :: time_steps, output_cnt
  real(dp)           :: time, end_time
  real(dp)           :: dt
  real(dp)           :: dt_output
  character(len=100) :: fname

  real(dp), parameter :: alpha = 0.1305_dp
  real(dp), parameter :: beta = 0.7695_dp
  real(dp), parameter :: D1 = 0.05_dp
  real(dp), parameter :: D2 = 1.0_dp
  real(dp), parameter :: kappa = 100.0_dp

  print *, "Running reaction_diffusion_" // DIMNAME // ""
  print *, "Number of threads", af_get_max_threads()

  call af_add_cc_variable(tree, "u", ix=i_u, n_copies=3)
  call af_add_cc_variable(tree, "v", ix=i_v, n_copies=3)
  call af_add_cc_variable(tree, "phi1", ix=i_phi1)
  call af_add_cc_variable(tree, "phi2", ix=i_phi2)
  call af_add_cc_variable(tree, "rhs1", ix=i_rhs1)
  call af_add_cc_variable(tree, "rhs2", ix=i_rhs2)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)

  call af_set_cc_methods(tree, i_u, af_bc_neumann_zero)
  call af_set_cc_methods(tree, i_v, af_bc_neumann_zero)

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [DTIMES(domain_len)], &
       [DTIMES(box_size)], &
       periodic=[DTIMES(.false.)])

  call af_print_info(tree)

  ! Set up the initial conditions
  do
     call af_loop_box(tree, set_initial_condition)
     call af_gc_tree(tree, [i_u, i_v])
     call af_adjust_refinement(tree, &           ! tree
          refine_routine, & ! Refinement function
          refine_info)      ! Information about refinement
     if (refine_info%n_add == 0) exit
  end do

  mg1%i_phi            = i_phi1
  mg1%i_rhs            = i_rhs1
  mg1%i_tmp            = i_tmp

  mg1%sides_bc    => af_bc_neumann_zero
  mg1%box_op      => mg_box_lpl
  mg1%box_gsrb    => mg_box_gsrb_lpl
  mg1%box_stencil => mg_box_lpl_stencil

  mg2%i_phi            = i_phi2
  mg2%i_rhs            = i_rhs2
  mg2%i_tmp            = i_tmp

  mg2%sides_bc    => af_bc_neumann_zero
  mg2%box_op      => mg_box_lpl
  mg2%box_gsrb    => mg_box_gsrb_lpl
  mg2%box_stencil => mg_box_lpl_stencil


  output_cnt = 0
  time       = 0
  end_time   = 2.0_dp
  time_steps = 0
  time       = 0
  dt         = 0.5 * af_min_dr(tree)**2/(2*NDIM*max(D1, D2))
  dt         = 1e-3_dp
  dt_output  = 1e-1_dp

  mg1%helmholtz_lambda = 1/(0.5_dp * dt * D1)
  mg2%helmholtz_lambda = 1/(0.5_dp * dt * D2)

  call mg_init(tree, mg1)
  call mg_init(tree, mg2)

  ! Starting simulation
  do
     time_steps = time_steps + 1

     if (output_cnt * dt_output <= time) then
        output_cnt = output_cnt + 1
        write(fname, "(A,I0)") "reaction_diffusion_" // DIMNAME // "_", output_cnt
        call af_write_silo(tree, trim(fname), output_cnt, time, dir="output")
     end if

     ! call af_loop_box(tree, gs_explicit)
     call gs_imex(tree)
     call af_gc_tree(tree, [i_u, i_v])

     if (time > end_time) exit
     time = time + dt
  end do

contains

  subroutine gs_explicit(box)
    type(box_t), intent(inout) :: box
    integer, parameter         :: explicit_order = 1

    select case (explicit_order)
    case (1)
       call step_F(box, dt, [i_u, i_v], [i_u, i_v], [i_u, i_v])
    case (2)
       ! RK2 method
       call step_F(box, 0.5_dp * dt, [i_u, i_v], [i_u, i_v], [i_u+1, i_v+1])
       call step_F(box, dt, [i_u+1, i_v+1], [i_u, i_v], [i_u, i_v])
    end select
  end subroutine gs_explicit

  subroutine gs_imex(tree)
    type(af_t), intent(inout) :: tree
    integer                   :: n
    real(dp)                  :: max_res, max_rhs
    real(dp), parameter       :: max_rel_residual = 1e-6_dp

    call af_loop_box(tree, gs_imex_step1)

    do n = 1, 10
       call af_tree_maxabs_cc(tree, mg1%i_rhs, max_rhs)
       call mg_fas_fmg(tree, mg1, .true., .true.)
       ! call mg_fas_vcycle(tree, mg1, .true.)
       call af_tree_maxabs_cc(tree, mg1%i_tmp, max_res)
       if (max_res/max_rhs < max_rel_residual) exit
    end do

    do n = 1, 10
       call af_tree_maxabs_cc(tree, mg2%i_rhs, max_rhs)
       call mg_fas_fmg(tree, mg2, .true., .true.)
       ! call mg_fas_vcycle(tree, mg2, .true.)
       call af_tree_maxabs_cc(tree, mg2%i_tmp, max_res)
       if (max_res/max_rhs < max_rel_residual) exit
    end do

    call af_loop_box(tree, gs_imex_step2)
  end subroutine gs_imex

  subroutine gs_imex_step1(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc

    nc = box%n_cell

    ! Set rhs for Helmholtz equation
    call step_F0(box, dt, [i_u, i_v], [i_u, i_v], [mg1%i_rhs, mg2%i_rhs])
    call step_F1(box, 0.5_dp * dt, [i_u, i_v], [mg1%i_rhs, mg2%i_rhs], &
         [mg1%i_rhs, mg2%i_rhs])

    box%cc(DTIMES(1:nc), mg1%i_rhs) = -box%cc(DTIMES(1:nc), mg1%i_rhs) * &
         mg1%helmholtz_lambda
    box%cc(DTIMES(1:nc), mg2%i_rhs) = -box%cc(DTIMES(1:nc), mg2%i_rhs) * &
         mg2%helmholtz_lambda

  end subroutine gs_imex_step1

  subroutine gs_imex_step2(box)
    type(box_t), intent(inout) :: box

    call step_F(box, 0.5_dp * dt, [i_u, i_v], [i_u, i_v], [i_u, i_v])
    call step_F(box, 0.5_dp * dt, [mg1%i_phi, mg2%i_phi], &
         [i_u, i_v], [i_u, i_v])
  end subroutine gs_imex_step2

  subroutine step_F0(box, dt, i_deriv, i_prev, i_out)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: dt
    integer, intent(in)        :: i_deriv(2)
    integer, intent(in)        :: i_prev(2)
    integer, intent(in)        :: i_out(2)
    real(dp)                   :: tmp(DTIMES(box%n_cell), 2)
    integer                    :: nc

    nc = box%n_cell

    ! Store the derivatives in case i_out overlaps with i_deriv
    tmp(DTIMES(:), 1) = kappa * (alpha - box%cc(DTIMES(1:nc), i_deriv(1)) + &
         box%cc(DTIMES(1:nc), i_deriv(1))**2 * box%cc(DTIMES(1:nc), i_deriv(2)))
    tmp(DTIMES(:), 2) = kappa * (beta - &
         box%cc(DTIMES(1:nc), i_deriv(1))**2 * box%cc(DTIMES(1:nc), i_deriv(2)))

    box%cc(DTIMES(1:nc), i_out) = &
         box%cc(DTIMES(1:nc), i_prev) + dt * tmp
  end subroutine step_F0

  subroutine step_F1(box, dt, i_deriv, i_prev, i_out)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: dt
    integer, intent(in)        :: i_deriv(2)
    integer, intent(in)        :: i_prev(2)
    integer, intent(in)        :: i_out(2)
    real(dp)                   :: tmp(DTIMES(box%n_cell), 2)
    integer                    :: nc

    nc = box%n_cell
    call laplacian(box, i_deriv(1), tmp(DTIMES(:), 1))
    call laplacian(box, i_deriv(2), tmp(DTIMES(:), 2))
    tmp(DTIMES(:), 1) = tmp(DTIMES(:), 1) * D1
    tmp(DTIMES(:), 2) = tmp(DTIMES(:), 2) * D2

    box%cc(DTIMES(1:nc), i_out) = &
         box%cc(DTIMES(1:nc), i_prev) + dt * tmp
  end subroutine step_F1

  subroutine step_F(box, dt, i_deriv, i_prev, i_out)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: dt
    integer, intent(in)        :: i_deriv(2)
    integer, intent(in)        :: i_prev(2)
    integer, intent(in)        :: i_out(2)
    integer                    :: nc
    real(dp)                   :: tmp(DTIMES(box%n_cell), 2)

    nc = box%n_cell

    call laplacian(box, i_deriv(1), tmp(DTIMES(:), 1))
    call laplacian(box, i_deriv(2), tmp(DTIMES(:), 2))
    tmp(DTIMES(:), 1) = tmp(DTIMES(:), 1) * D1
    tmp(DTIMES(:), 2) = tmp(DTIMES(:), 2) * D2

    ! Add the other derivatives
    tmp(DTIMES(:), 1) = tmp(DTIMES(:), 1) + &
         kappa * (alpha - box%cc(DTIMES(1:nc), i_deriv(1)) + &
         box%cc(DTIMES(1:nc), i_deriv(1))**2 * box%cc(DTIMES(1:nc), i_deriv(2)))
    tmp(DTIMES(:), 2) = tmp(DTIMES(:), 2) + &
         kappa * (beta - &
         box%cc(DTIMES(1:nc), i_deriv(1))**2 * box%cc(DTIMES(1:nc), i_deriv(2)))

    box%cc(DTIMES(1:nc), i_out) = box%cc(DTIMES(1:nc), i_prev) + dt * tmp
  end subroutine step_F

  !> Return the refinement flag for box
  subroutine refine_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))

    if (maxval(box%dr) > domain_len/50) then
       cell_flags = af_do_ref
    else
       cell_flags = af_keep_ref
    end if
  end subroutine refine_routine

  !> This routine sets the initial conditions for each box
  subroutine set_initial_condition(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM), r0(NDIM)

    nc = box%n_cell
    r0(1) = 1.0_dp/3
    r0(2:) = 0.5_dp

    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, i_u) = alpha + beta + 1e-3_dp * &
            exp(-100_dp * sum((rr - r0)**2))
       box%cc(IJK, i_v) = beta / (alpha + beta)**2
    end do; CLOSE_DO
  end subroutine set_initial_condition

  ! !> This routine computes the right hand side per box
  ! subroutine set_rhs(box)
  !   type(box_t), intent(inout) :: box
  !   integer                     :: nc

  !   nc = box%n_cell
  !   box%cc(DTIMES(1:nc), i_rhs) = -1/(dt * diffusion_coeff) * &
  !        box%cc(DTIMES(1:nc), i_phi)
  ! end subroutine set_rhs

  !> Perform Laplacian operator
  subroutine laplacian(box, i_in, lpl)
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: i_in
    real(dp), intent(out)      :: lpl(DTIMES(1:box%n_cell))
    integer                    :: IJK, nc
    real(dp)                   :: idr2(NDIM)

    nc   = box%n_cell
    idr2 = 1 / box%dr**2

    associate (cc => box%cc, n => i_in)
      do KJI_DO(1, nc)
#if NDIM == 2
         lpl(i, j) = &
              idr2(1) * (cc(i-1, j, n) + cc(i+1, j, n) - 2 * cc(i, j, n)) + &
              idr2(2) * (cc(i, j-1, n) + cc(i, j+1, n) - 2 * cc(i, j, n))
#elif NDIM == 3
         lpl(i, j, k) = &
              idr2(1) * (cc(i-1, j, k, n) + cc(i+1, j, k, n) &
              - 2 * cc(i, j, k, n)) &
              + idr2(2) * (cc(i, j-1, k, n) + cc(i, j+1, k, n) &
              - 2 * cc(i, j, k, n)) &
              + idr2(3) * (cc(i, j, k-1, n) + cc(i, j, k+1, n) &
              - 2 * cc(i, j, k, n))
#endif
      end do; CLOSE_DO
    end associate
  end subroutine laplacian

end program reaction_diffusion_Xd
