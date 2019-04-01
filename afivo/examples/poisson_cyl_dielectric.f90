!> \example poisson_cyl_dielectric.f90
!>
!> Example showing how to use m_af_multigrid in cylindrical coordinates with an abrubt
!> change in "eps", and compare with an analytic solution.
program poisson_cyl_dielectric
  use m_af_all
  use m_gaussians

  implicit none

  integer, parameter :: box_size     = 8
  integer, parameter :: n_iterations = 10
  integer            :: i_phi
  integer            :: i_rhs
  integer            :: i_err
  integer            :: i_tmp
  integer            :: i_eps

  type(af_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: mg_iter
  real(dp)           :: residu(2), anal_err(2)
  character(len=100) :: fname
  type(mg_t)        :: mg
  type(gauss_t)      :: gs
  integer            :: count_rate,t_start, t_end

  print *, "Running poisson_cyl_dielectric"
  print *, "Number of threads", af_get_max_threads()

  ! The manufactured solution exists of two Gaussians, which are stored in gs
  call gauss_init(gs, [1.0_dp], [0.1_dp], &
       reshape([0.0_dp, 0.25_dp], [2,1]))

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)
  call af_add_cc_variable(tree, "err", ix=i_err)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)
  call af_add_cc_variable(tree, "eps", ix=i_eps)

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [1.0_dp, 1.0_dp], &
       [box_size, box_size], &
       coord=af_cyl)   ! Cylindrical coordinates

  call af_print_info(tree)

  call system_clock(t_start, count_rate)
  do
     ! For each box, set the initial conditions
     call af_loop_box(tree, set_init_cond)

     ! This updates the refinement of the tree, by at most one level per call.
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
  mg%i_eps        = i_eps       ! Variable for epsilon coefficient
  mg%sides_bc     => sides_bc   ! Method for boundary conditions Because we use

  ! Initialize the multigrid options. This performs some basics checks and sets
  ! default values where necessary.
  ! This routine does not initialize the multigrid variables i_phi, i_rhs
  ! and i_tmp. These variables will be initialized at the first call of mg_fas_fmg
  call mg_init(tree, mg)

  print *, "Multigrid iteration | max residual | max error"
  call system_clock(t_start, count_rate)

  do mg_iter = 1, n_iterations
     ! Perform a FAS-FMG cycle (full approximation scheme, full multigrid). The
     ! third argument controls whether the residual is stored in i_tmp. The
     ! fourth argument controls whether to improve the current solution.
     call mg_fas_fmg(tree, mg, .true., mg_iter>1)

     ! Compute the error compared to the analytic solution
     call af_loop_box(tree, set_err)

     ! Determine the minimum and maximum residual and error
     call af_tree_min_cc(tree, i_tmp, residu(1))
     call af_tree_max_cc(tree, i_tmp, residu(2))
     call af_tree_min_cc(tree, i_err, anal_err(1))
     call af_tree_max_cc(tree, i_err, anal_err(2))
     write(*,"(I8,2Es14.5)") mg_iter, maxval(abs(residu)), &
          maxval(abs(anal_err))

     write(fname, "(A,I0)") "poisson_cyl_dielectric_", mg_iter
     call af_write_silo(tree, trim(fname), dir="output")
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

  ! Set refinement flags for box
  subroutine ref_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(box%n_cell, box%n_cell)
    integer                  :: i, j, nc
    real(dp)                 :: crv, dr2

    nc = box%n_cell
    dr2 = maxval(box%dr)**2

    ! Compute the "curvature" in phi
    do j = 1, nc
       do i = 1, nc
          crv = dr2 * abs(box%cc(i, j, i_rhs)) / box%cc(i, j, i_eps)

          ! And refine if it exceeds a threshold
          if (crv > 1.0e-3_dp) then
             cell_flags(i, j) = af_do_ref
          else
             cell_flags(i, j) = af_keep_ref
          end if
       end do
    end do
  end subroutine ref_routine

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: rz(2), grad(2), qbnd, tmp

    nc                  = box%n_cell
    box%cc(:, :, i_phi) = 0

    do j = 0, nc+1
       do i = 0, nc+1
          rz = af_r_cc(box, [i,j])

          ! Change epsilon in part of the domain
          if (rz(1) < 0.5_dp .and. rz(2) < 0.5_dp) then
             box%cc(i, j, i_eps) = 100.0_dp
          else
             box%cc(i, j, i_eps) = 1.0_dp
          end if

          ! Partially compute the right-hand side (see below)
          box%cc(i, j, i_rhs) = gauss_laplacian_cyl(gs, rz) * box%cc(i, j, i_eps)
       end do
    end do

    ! We have to place surface charges where epsilon has a jump, this is first
    ! done in the r-direction
    do j = 1, nc
       do i = 0, nc
          rz = af_rr_cc(box, [i + 0.5_dp, real(j, dp)])

          ! Determine amount of charge
          call gauss_gradient(gs, rz, grad)
          qbnd = (box%cc(i+1, j, i_eps) - box%cc(i, j, i_eps)) * &
               grad(1) / box%dr(1)

          ! Place surface charge weighted with eps
          tmp = box%cc(i+1, j, i_eps) / &
               (box%cc(i, j, i_eps) + box%cc(i+1, j, i_eps))
          box%cc(i+1, j, i_rhs) = box%cc(i+1, j, i_rhs) + tmp * qbnd
          box%cc(i, j, i_rhs) = box%cc(i, j, i_rhs) + (1-tmp) * qbnd
       end do
    end do

    ! Set surface charge in z-direction
    do j = 0, nc
       do i = 1, nc
          rz = af_rr_cc(box, [real(i, dp), j + 0.5_dp])

          ! Determine amount of charge
          call gauss_gradient(gs, rz, grad)
          qbnd = (box%cc(i, j+1, i_eps) - box%cc(i, j, i_eps)) * &
               grad(2) / box%dr(2)

          ! Place surface charge weighted with eps
          tmp = box%cc(i, j+1, i_eps) / &
               (box%cc(i, j, i_eps) + box%cc(i, j+1, i_eps))
          box%cc(i, j+1, i_rhs) = box%cc(i, j+1, i_rhs) + tmp * qbnd
          box%cc(i, j, i_rhs) = box%cc(i, j, i_rhs) + (1-tmp) * qbnd
       end do
    end do
  end subroutine set_init_cond

  ! Compute error compared to the analytic solution
  subroutine set_err(box)
    type(box_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: rz(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          rz = af_r_cc(box, [i,j])
          box%cc(i, j, i_err) = box%cc(i, j, i_phi) - gauss_value(gs, rz)
       end do
    end do
  end subroutine set_err

  ! This routine sets boundary conditions for a box
  subroutine sides_bc(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type
    integer                 :: n

    if (nb == af_neighb_lowx) then
       ! On the axis, apply Neumann zero conditions
       bc_type = af_bc_neumann
       bc_val = 0.0_dp
    else
       ! We use dirichlet boundary conditions
       bc_type = af_bc_dirichlet

       ! Below the solution is specified in the approriate ghost cells
       do n = 1, box%n_cell**(NDIM-1)
          bc_val(n) = gauss_value(gs, coords(:, n))
       end do
    end if
  end subroutine sides_bc

end program poisson_cyl_dielectric
