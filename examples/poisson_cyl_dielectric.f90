!> \example poisson_cyl_dielectric.f90
!>
!> Example showing how to use m_a2_multigrid in cylindrical coordinates with an abrubt
!> change in "eps", and compare with an analytic solution.
program poisson_cyl_dielectric
  use m_a2_all
  use m_gaussians

  implicit none

  integer, parameter :: box_size = 8
  integer, parameter :: n_boxes_base = 1
  integer, parameter :: n_iterations = 10
  integer, parameter :: n_var_cell = 5
  integer, parameter :: i_phi = 1
  integer, parameter :: i_rhs = 2
  integer, parameter :: i_err = 3
  integer, parameter :: i_tmp = 4
  integer, parameter :: i_eps = 5

  type(a2_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: mg_iter
  integer            :: ix_list(2, n_boxes_base)
  real(dp)           :: dr, residu(2), anal_err(2)
  character(len=100) :: fname
  type(mg2_t)        :: mg
  type(gauss_t)      :: gs
  integer            :: count_rate,t_start, t_end

  print *, "Running poisson_cyl_dielectric"
  print *, "Number of threads", af_get_max_threads()

  ! The manufactured solution exists of two Gaussians, which are stored in gs
  call gauss_init(gs, [1.0_dp, 1.0_dp], [0.04_dp, 0.04_dp], &
       reshape([0.25_dp, 0.25_dp, 0.75_dp, 0.75_dp], [2,2]))

  ! The cell spacing at the coarsest grid level
  dr = 1.0_dp / box_size

  ! Initialize tree
  call a2_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       n_var_cell, &   ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       coarsen_to=2, & ! Add coarsened levels for multigrid
       coord=af_cyl, & ! Cylindrical coordinates
       cc_names=["phi", "rhs", "err", "tmp", "eps"]) ! Variable names

  ! Set up geometry. These indices are used to define the coordinates of a box,
  ! by default the box at [1,1] touches the origin (x,y) = (0,0)
  ix_list(:, 1) = [1,1]         ! Set index of box 1

  ! Create the base mesh, using the box indices and their neighbor information
  call a2_set_base(tree, 1, ix_list)
  call a2_print_info(tree)

  call system_clock(t_start, count_rate)
  do
     ! For each box, set the initial conditions
     call a2_loop_box(tree, set_init_cond)

     ! This updates the refinement of the tree, by at most one level per call.
     call a2_adjust_refinement(tree, ref_routine, ref_info)

     ! If no new boxes have been added, exit the loop
     if (ref_info%n_add == 0) exit
  end do
  call system_clock(t_end, count_rate)

  write(*,"(A,Es10.3,A)") " Wall-clock time generating AMR grid: ", &
       (t_end-t_start) / real(count_rate,dp), " seconds"

  call a2_print_info(tree)

  ! Set the multigrid options.
  mg%i_phi        = i_phi       ! Solution variable
  mg%i_rhs        = i_rhs       ! Right-hand side variable
  mg%i_tmp        = i_tmp       ! Variable for temporary space
  mg%i_eps        = i_eps       ! Variable for epsilon coefficient
  mg%sides_bc     => sides_bc   ! Method for boundary conditions Because we use

  ! Automatically detect the right methods
  mg%box_op       => mg2_auto_op
  mg%box_gsrb     => mg2_auto_gsrb
  mg%box_corr     => mg2_auto_corr

  ! Initialize the multigrid options. This performs some basics checks and sets
  ! default values where necessary.
  ! This routine does not initialize the multigrid variables i_phi, i_rhs
  ! and i_tmp. These variables will be initialized at the first call of mg2_fas_fmg
  call mg2_init_mg(mg)

  print *, "Multigrid iteration | max residual | max error"
  call system_clock(t_start, count_rate)

  do mg_iter = 1, n_iterations
     ! Perform a FAS-FMG cycle (full approximation scheme, full multigrid). The
     ! third argument controls whether the residual is stored in i_tmp. The
     ! fourth argument controls whether to improve the current solution.
     call mg2_fas_fmg(tree, mg, .true., mg_iter>1)

     ! Compute the error compared to the analytic solution
     call a2_loop_box(tree, set_err)

     ! Determine the minimum and maximum residual and error
     call a2_tree_min_cc(tree, i_tmp, residu(1))
     call a2_tree_max_cc(tree, i_tmp, residu(2))
     call a2_tree_min_cc(tree, i_err, anal_err(1))
     call a2_tree_max_cc(tree, i_err, anal_err(2))
     write(*,"(I8,2Es14.5)") mg_iter, maxval(abs(residu)), &
          maxval(abs(anal_err))

     write(fname, "(A,I0)") "poisson_cyl_dielectric_", mg_iter
     call a2_write_vtk(tree, trim(fname), dir="output")
  end do
  call system_clock(t_end, count_rate)

  write(*, "(A,I0,A,E10.3,A)") &
       " Wall-clock time after ", n_iterations, &
       " iterations: ", (t_end-t_start) / real(count_rate, dp), &
       " seconds"

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call a2_destroy(tree)

contains

  ! Set refinement flags for box
  subroutine ref_routine(box, cell_flags)
    type(box2_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(box%n_cell, box%n_cell)
    integer                  :: i, j, nc
    real(dp)                 :: crv

    nc = box%n_cell

    ! Compute the "curvature" in phi
    do j = 1, nc
       do i = 1, nc
          crv = box%dr**2 * abs(box%cc(i, j, i_rhs)) / box%cc(i, j, i_eps)

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
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: rz(2), grad(2), dr, qbnd, tmp

    nc                  = box%n_cell
    box%cc(:, :, i_phi) = 0
    dr                  = box%dr

    do j = 0, nc+1
       do i = 0, nc+1
          rz = a2_r_cc(box, [i,j])

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
          rz = a2_rr_cc(box, [i + 0.5_dp, real(j, dp)])

          ! Determine amount of charge
          call gauss_gradient(gs, rz, grad)
          qbnd = (box%cc(i+1, j, i_eps) - box%cc(i, j, i_eps)) * &
               grad(1) / dr

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
          rz = a2_rr_cc(box, [real(i, dp), j + 0.5_dp])

          ! Determine amount of charge
          call gauss_gradient(gs, rz, grad)
          qbnd = (box%cc(i, j+1, i_eps) - box%cc(i, j, i_eps)) * &
               grad(2) / dr

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
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: rz(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          rz = a2_r_cc(box, [i,j])
          box%cc(i, j, i_err) = box%cc(i, j, i_phi) - gauss_value(gs, rz)
       end do
    end do
  end subroutine set_err

  ! This routine sets boundary conditions for a box, by filling its ghost cells
  ! with approriate values. Note that on the axis (a boundary in the lower-x
  ! direction) we should use a Neumann zero condition in cylindrical
  ! coordinates.
  subroutine sides_bc(box, nb, iv, bc_type)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: nb ! Direction for the boundary condition
    integer, intent(in)         :: iv ! Index of variable
    integer, intent(out)        :: bc_type ! Type of boundary condition
    real(dp)                    :: rz(2)
    integer                     :: n, nc

    nc = box%n_cell

    select case (nb)
    case (a2_neighb_lowx)             ! Neumann zero on axis
       bc_type = af_bc_neumann
       box%cc(0, 1:nc, iv) = 0
    case (a2_neighb_highx)             ! Use solution on other boundaries
       bc_type = af_bc_dirichlet
       do n = 1, nc
          rz = a2_rr_cc(box, [nc+0.5_dp, real(n, dp)])
          box%cc(nc+1, n, iv) = gauss_value(gs, rz)
       end do
    case (a2_neighb_lowy)
       bc_type = af_bc_dirichlet
       do n = 1, nc
          rz = a2_rr_cc(box, [real(n, dp), 0.5_dp])
          box%cc(n, 0, iv) = gauss_value(gs, rz)
       end do
    case (a2_neighb_highy)
       bc_type = af_bc_dirichlet
       do n = 1, nc
          rz = a2_rr_cc(box, [real(n, dp), nc+0.5_dp])
          box%cc(n, nc+1, iv) = gauss_value(gs, rz)
       end do
    end select
  end subroutine sides_bc

end program poisson_cyl_dielectric
