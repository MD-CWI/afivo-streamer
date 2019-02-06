#include "../src/cpp_macros.h"
!> \example poisson_helmholtz_Xd.f90
!>
!> Example showing how to use multigrid for Helmholtz equation and compare with
!> an analytic solution. A standard 5 or 7 point stencil is used.
program poisson_helmholtz_Xd
  use m_af_all
  use m_gaussians

  implicit none

  integer, parameter  :: box_size     = 8
  integer, parameter  :: n_iterations = 10
  integer             :: i_phi
  integer             :: i_rhs
  integer             :: i_err
  integer             :: i_tmp
  real(dp), parameter :: lambda       = 1.0e3_dp

  type(af_t)        :: tree
  type(ref_info_t)   :: refine_info
  integer            :: mg_iter
  real(dp)           :: dr, residu(2), anal_err(2)
  character(len=100) :: fname
  type(mg_t)       :: mg
  type(gauss_t)      :: gs
  integer            :: count_rate,t_start,t_end

  print *, "Running poisson_helmholtz_" // DIMNAME // ""
  print *, "Number of threads", af_get_max_threads()

  call gauss_init(gs, [1.0_dp, 1.0_dp], [0.04_dp, 0.04_dp], &
       reshape([DTIMES(0.25_dp), &
       DTIMES(0.75_dp)], [NDIM,2]))

  ! The cell spacing at the coarsest grid level
  dr = 1.0_dp / box_size

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)
  call af_add_cc_variable(tree, "err", ix=i_err)

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       dr)             ! Distance between cells on base level

  call af_set_coarse_grid(tree, [DTIMES(box_size)])

  do
     call af_loop_box(tree, set_initial_condition)
     call af_adjust_refinement(tree, refine_routine, refine_info)
     if (refine_info%n_add == 0) exit
  end do

  call af_print_info(tree)

  mg%i_phi    =  i_phi     ! Solution variable
  mg%i_rhs    =  i_rhs     ! Right-hand side variable
  mg%i_tmp    =  i_tmp     ! Variable for temporary space
  mg%sides_bc => sides_bc ! Method for boundary conditions
  mg%box_op   => helmholtz_operator
  mg%box_gsrb => helmholtz_gsrb

  call mg_init_mg(tree, mg)

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

     write(fname, "(A,I0)") "poisson_helmholtz_" // DIMNAME // "_", mg_iter
     call af_write_vtk(tree, trim(fname), dir="output")
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
    integer                  :: IJK, nc
    real(dp)                 :: rr(NDIM), dr2, drhs

    nc = box%n_cell
    dr2 = box%dr**2

    do KJI_DO(1,nc)
       rr = af_r_cc(box, [IJK])

       ! This is an estimate of the truncation error in the right-hand side,
       ! which is related to the fourth derivative of the solution.
       drhs = dr2 * box%cc(IJK, i_rhs)

       if (abs(drhs) > 1e-3_dp .and. box%lvl < 5) then
          cell_flags(IJK) = af_do_ref
       else
          cell_flags(IJK) = af_keep_ref
       end if
    end do; CLOSE_DO
  end subroutine refine_routine

  ! This routine sets the initial conditions for each box
  subroutine set_initial_condition(box)
    type(box_t), intent(inout) :: box
    integer                     :: IJK, nc
    real(dp)                    :: rr(NDIM)

    nc = box%n_cell

    do KJI_DO(1,nc)
       ! Get the coordinate of the cell center at IJK
       rr = af_r_cc(box, [IJK])

       ! And set the rhs values
       box%cc(IJK, i_rhs) = gauss_laplacian(gs, rr) - &
            lambda * gauss_value(gs, rr)
    end do; CLOSE_DO
  end subroutine set_initial_condition

  ! Set the error compared to the analytic solution
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

  ! This routine sets boundary conditions for a box, by filling its ghost cells
  ! with approriate values.
  subroutine sides_bc(box, nb, iv, bc_type)
    type(box_t), intent(inout) :: box
    integer, intent(in)          :: nb      ! Direction for the boundary condition
    integer, intent(in)          :: iv      ! Index of variable
    integer, intent(out)         :: bc_type ! Type of boundary condition
    real(dp)                     :: rr(NDIM)
#if NDIM == 2
    integer                      :: n, nc
#elif NDIM == 3
    integer                      :: IJK, ix, nc
    real(dp)                     :: loc
#endif

    nc = box%n_cell

    ! We use dirichlet boundary conditions
    bc_type = af_bc_dirichlet

    ! Below the solution is specified in the approriate ghost cells
#if NDIM == 2
    select case (nb)
    case (af_neighb_lowx)             ! Lower-x direction
       do n = 1, nc
          rr = af_rr_cc(box, [0.5_dp, real(n, dp)])
          box%cc(0, n, iv) = gauss_value(gs, rr)
       end do
    case (af_neighb_highx)             ! Higher-x direction
       do n = 1, nc
          rr = af_rr_cc(box, [nc+0.5_dp, real(n, dp)])
          box%cc(nc+1, n, iv) = gauss_value(gs, rr)
       end do
    case (af_neighb_lowy)             ! Lower-y direction
       do n = 1, nc
          rr = af_rr_cc(box, [real(n, dp), 0.5_dp])
          box%cc(n, 0, iv) = gauss_value(gs, rr)
       end do
    case (af_neighb_highy)             ! Higher-y direction
       do n = 1, nc
          rr = af_rr_cc(box, [real(n, dp), nc+0.5_dp])
          box%cc(n, nc+1, iv) = gauss_value(gs, rr)
       end do
    end select
#elif NDIM == 3
    ! Determine whether the direction nb is to "lower" or "higher" neighbors
    if (af_neighb_low(nb)) then
       ix = 0
       loc = 0.5_dp
    else
       ix = nc+1
       loc = nc + 0.5_dp
    end if

    ! Below the solution is specified in the approriate ghost cells
    select case (af_neighb_dim(nb))
    case (1)
       do k = 1, nc
          do j = 1, nc
             rr = af_rr_cc(box, [loc, real(j, dp), real(k, dp)])
             box%cc(ix, j, k, iv) = gauss_value(gs, rr)
          end do
       end do
    case (2)
       do k = 1, nc
          do i = 1, nc
             rr = af_rr_cc(box, [real(i, dp), loc, real(k, dp)])
             box%cc(i, ix, k, iv) = gauss_value(gs, rr)
          end do
       end do
    case (3)
       do j = 1, nc
          do i = 1, nc
             rr = af_rr_cc(box, [real(i, dp), real(j, dp), loc])
             box%cc(i, j, ix, iv) = gauss_value(gs, rr)
          end do
       end do
    end select
#endif
  end subroutine sides_bc

  !> Perform Helmholtz operator on a box
  subroutine helmholtz_operator(box, i_out, mg)
    type(box_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    type(mg_t), intent(in)     :: mg !< Multigrid options
    integer                     :: i, j, nc, i_phi
    real(dp)                    :: inv_dr_sq
#if NDIM == 3
    integer                     :: k
#endif

    nc = box%n_cell
    inv_dr_sq = 1 / box%dr**2
    i_phi = mg%i_phi

#if NDIM == 2
    do j = 1, nc
       do i = 1, nc
          box%cc(i, j, i_out) = inv_dr_sq * (box%cc(i-1, j, i_phi) + &
               box%cc(i+1, j, i_phi) + box%cc(i, j-1, i_phi) + &
               box%cc(i, j+1, i_phi) - 4 * box%cc(i, j, i_phi)) - &
               lambda * box%cc(i, j, i_phi)
       end do
    end do
#elif NDIM == 3
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             box%cc(i, j, k, i_out) = inv_dr_sq * (box%cc(i-1, j, k, i_phi) + &
                  box%cc(i+1, j, k, i_phi) + box%cc(i, j-1, k, i_phi) + &
                  box%cc(i, j+1, k, i_phi) + box%cc(i, j, k-1, i_phi) + &
                  box%cc(i, j, k+1, i_phi) - 6 * box%cc(i, j, k, i_phi)) - &
                  lambda * box%cc(i, j, k, i_phi)
          end do
       end do
    end do
#endif
  end subroutine helmholtz_operator

  !> Perform Gauss-Seidel relaxation on box for a Helmholtz operator
  subroutine helmholtz_gsrb(box, redblack_cntr, mg)
    type(box_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    type(mg_t), intent(in)     :: mg !< Multigrid options
    integer                     :: i, i0, j, nc, i_phi, i_rhs
    real(dp)                    :: dx2
#if NDIM == 3
    integer                     :: k
#endif

    dx2   = box%dr**2
    nc    = box%n_cell
    i_phi = mg%i_phi
    i_rhs = mg%i_rhs

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
#if NDIM == 2
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          box%cc(i, j, i_phi) = 1/(4 + lambda * dx2) * ( &
               box%cc(i+1, j, i_phi) + box%cc(i-1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi) - &
               dx2 * box%cc(i, j, i_rhs))
       end do
    end do
#elif NDIM == 3
    do k = 1, nc
       do j = 1, nc
          i0 = 2 - iand(ieor(redblack_cntr, k+j), 1)
          do i = i0, nc, 2
             box%cc(i, j, k, i_phi) = 1/(6 + lambda * dx2) * ( &
                  box%cc(i+1, j, k, i_phi) + box%cc(i-1, j, k, i_phi) + &
                  box%cc(i, j+1, k, i_phi) + box%cc(i, j-1, k, i_phi) + &
                  box%cc(i, j, k+1, i_phi) + box%cc(i, j, k-1, i_phi) - &
                  dx2 * box%cc(i, j, k, i_rhs))
          end do
       end do
    end do
#endif
  end subroutine helmholtz_gsrb

end program poisson_helmholtz_Xd


