#include "../src/cpp_macros.h"
!> \example poisson_basic_Xd.f90
!>
!> Example showing how to use multigrid and compare with an analytic solution. A
!> standard 5-point Laplacian is used.
program poisson_basic
  use m_af_all
  use m_gaussians

  implicit none

  integer, parameter :: box_size = 8
  integer, parameter :: n_iterations = 10
  integer :: i_phi
  integer :: i_rhs
  integer :: i_err
  integer :: i_tmp
  integer :: i_gradx
  integer :: i_egradx

  type(af_t)        :: tree
  type(ref_info_t)   :: refine_info
  integer            :: mg_iter
  real(dp)           :: dr, residu(2), anal_err(2)
  character(len=100) :: fname
  type(mg_t)       :: mg
  type(gauss_t)      :: gs
  integer            :: count_rate,t_start,t_end

  print *, "Running poisson_basic_" // DIMNAME
  print *, "Number of threads", af_get_max_threads()

  !> [Gauss_init]
  ! The manufactured solution exists of Gaussians
  ! Amplitudes:  [1.0_dp, 1.0_dp]
  ! Sigmas    :  [0.04_dp, 0.04_dp]
  ! Locations :  x, y, z = 0.25 or x, y, z = 0.75
  call gauss_init(gs, [1.0_dp, 1.0_dp], [0.04_dp, 0.04_dp], &
       reshape([DTIMES(0.25_dp), &
       DTIMES(0.75_dp)], [NDIM,2]))
  !> [Gauss_init]

  ! The cell spacing at the coarsest grid level
  dr = 1.0_dp / box_size

  !> [af_init]
  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)
  call af_add_cc_variable(tree, "err", ix=i_err)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)
  call af_add_cc_variable(tree, "Dx",  ix=i_gradx)
  call af_add_cc_variable(tree, "eDx", ix=i_egradx)

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [DTIMES(1.0_dp)], &
       [DTIMES(box_size)])
  !> [af_init]

  call system_clock(t_start, count_rate)
  !> [set_refinement]
  do
     ! For each box, set the initial conditions
     call af_loop_box(tree, set_initial_condition)

     ! This updates the refinement of the tree, by at most one level per call.
     ! The second argument is a subroutine that is called for each box that can
     ! be refined or derefined, and it should set refinement flags. Information
     ! about the changes in refinement are returned in the third argument.
     call af_adjust_refinement(tree, refine_routine, refine_info)

     ! If no new boxes have been added, exit the loop
     if (refine_info%n_add == 0) exit
  end do
  !> [set_refinement]
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
     ! Perform a FAS-FMG cycle (full approximation scheme, full multigrid). The
     ! third argument controls whether the residual is stored in i_tmp. The
     ! fourth argument controls whether to improve the current solution.
     call mg_fas_fmg(tree, mg, set_residual=.true., have_guess=(mg_iter>1))

     ! Compute the error compared to the analytic solution
     call af_loop_box(tree, set_error)

     ! Determine the minimum and maximum residual and error
     call af_tree_min_cc(tree, i_tmp, residu(1))
     call af_tree_max_cc(tree, i_tmp, residu(2))
     call af_tree_min_cc(tree, i_err, anal_err(1))
     call af_tree_max_cc(tree, i_err, anal_err(2))
     write(*,"(I8,13x,2(Es14.5))") mg_iter, maxval(abs(residu)), &
          maxval(abs(anal_err))

     !> [write_output]
     ! This writes a Silo output file containing the cell-centered values of the
     ! leaves of the tree (the boxes not covered by refinement).
     write(fname, "(A,I0)") "poisson_basic_" // DIMNAME // "_", mg_iter
     call af_write_silo(tree, trim(fname), dir="output")
     !> [write_output]
  end do
  call system_clock(t_end, count_rate)

  write(*, "(A,I0,A,E10.3,A)") &
       " Wall-clock time after ", n_iterations, &
       " iterations: ", (t_end-t_start) / real(count_rate, dp), &
       " seconds"

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call af_destroy(tree)

contains

  !> [refine_routine]
  ! Return the refinement flags for box
  subroutine refine_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(DTIMES(box%n_cell))
    integer                  :: IJK, nc
    real(dp)                 :: rr(NDIM), dr2, drhs

    nc = box%n_cell
    dr2 = maxval(box%dr)**2

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
  !> [refine_routine]

  !> [set_initial_condition]
  ! This routine sets the initial conditions for each box
  subroutine set_initial_condition(box)
    type(box_t), intent(inout) :: box
    integer                     :: IJK, nc
    real(dp)                    :: rr(NDIM), grad(NDIM)

    nc = box%n_cell

    do KJI_DO(1,nc)
       ! Get the coordinate of the cell center at IJK
       rr = af_r_cc(box, [IJK])

       ! And set the rhs values
       box%cc(IJK, i_rhs) = gauss_laplacian(gs, rr)
       call gauss_gradient(gs, rr, grad)
       box%cc(IJK, i_gradx) = grad(1)
    end do; CLOSE_DO
  end subroutine set_initial_condition
  !> [set_initial_condition]

  !> [set_error]
  ! Set the error compared to the analytic solution
  subroutine set_error(box)
    type(box_t), intent(inout) :: box
    integer                     :: IJK, nc
    real(dp)                    :: rr(NDIM), gradx

    nc = box%n_cell

    do KJI_DO(1,nc)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, i_err) = box%cc(IJK, i_phi) - gauss_value(gs, rr)
#if NDIM == 2
       gradx = 0.5_dp * (box%cc(i+1, j, i_phi) - &
            box%cc(i-1, j, i_phi)) / box%dr(1)
#elif NDIM == 3
       gradx = 0.5_dp * (box%cc(i+1, j, k, i_phi) - &
            box%cc(i-1, j, k, i_phi)) / box%dr(1)
#endif
       box%cc(IJK, i_egradx) = gradx - box%cc(IJK, i_gradx)

    end do; CLOSE_DO
  end subroutine set_error
  !> [set_error]

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

end program poisson_basic
