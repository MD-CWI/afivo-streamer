!> \example poisson_basic_2d.f90
!>
!> Example showing how to use multigrid and compare with an analytic solution. A
!> standard 5-point Laplacian is used.
program poisson_basic_2d
  use m_a2_types
  use m_a2_core
  use m_a2_multigrid
  use m_a2_utils
  use m_a2_output
  use m_gaussians

  implicit none

  integer, parameter :: box_size = 8
  integer, parameter :: n_boxes_base = 1
  integer, parameter :: n_iterations = 10
  integer, parameter :: n_var_cell = 4
  integer, parameter :: i_phi = 1
  integer, parameter :: i_rhs = 2
  integer, parameter :: i_err = 3
  integer, parameter :: i_tmp = 4

  type(a2_t)         :: tree
  type(ref_info_t)   :: refine_info
  integer            :: n_gaussian=2, mg_iter
  integer            :: ix_list(2, n_boxes_base)
  integer            :: nb_list(4, n_boxes_base)
  real(dp)           :: dr, residu(2), anal_err(2)
  character(len=100) :: fname
  type(mg2_t)        :: mg
  type(gauss_t)      :: gs
  integer            :: count_rate,t_start,t_end

  print *, "Running poisson_basic_2d"
  print *, "Number of threads", af_get_max_threads()

  ! The manufactured solution exists Gaussians, which are stored in gs
  if (n_gaussian == 1) then
     ! Amplitudes:  [1.0_dp]
     ! Sigmas    :  [0.04_dp]
     ! Locations :  [0.5_dp, 0.5_dp]
     call gauss_init(gs, [1.0_dp], [0.04_dp], &
          reshape([0.5_dp, 0.5_dp], [2,1]))
     write(*,"(2(A17,2x,i2,/),2(A17,1x,Es10.2,/),A17,1x,2(Es10.2,1x))") &
       "gs%n_gauss:",gs%n_gauss, &
       "gs%n_dim  :",gs%n_dim, &
       "gs%ampl   :",gs%ampl(:), &
       "gs%sigma  :",gs%sigma(:), &
       "gs%r0     :",gs%r0(:,:)
  else if (n_gaussian == 2) then
     ! Amplitudes:  [1.0_dp, 1.0_dp]
     ! Sigmas    :  [0.04_dp, 0.04_dp]
     ! Locations :  [[0.25_dp, 0.25_dp], [0.75_dp, 0.75_dp]]
     call gauss_init(gs, [1.0_dp, 1.0_dp], [0.04_dp, 0.04_dp], &
          reshape([0.25_dp, 0.25_dp, 0.75_dp, 0.75_dp], [2,2]))
  end if

  ! The cell spacing at the coarsest grid level
  dr = 1.0_dp / box_size

  ! Initialize tree
  call a2_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       n_var_cell, &   ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       coarsen_to=2, & ! Add coarsened levels for multigrid
       cc_names=["phi", "rhs", "err", "tmp"]) ! Variable names

  ! Set up geometry. These indices are used to define the coordinates of a box,
  ! by default the box at [1,1] touches the origin (x,y) = (0,0)
  ix_list(:, 1) = [1,1]         ! Set index of box 1

  ! Set neighbors for box one, negative values indicate a physical boundary
  nb_list(:, 1) = -1            ! Dirichlet zero -> -1

  ! Create the base mesh, using the box indices and their neighbor information
  call a2_set_base(tree, ix_list, nb_list)
  call a2_print_info(tree)

  call system_clock(t_start, count_rate)
  do
     ! For each box, set the initial conditions
     call a2_loop_box(tree, set_initial_condition)

     ! This updates the refinement of the tree, by at most one level per call.
     ! The second argument is a subroutine that is called for each box that can
     ! be refined or derefined, and it should set refinement flags. Information
     ! about the changes in refinement are returned in the third argument.
     call a2_adjust_refinement(tree, refine_routine, refine_info)

     ! If no new boxes have been added, exit the loop
     if (refine_info%n_add == 0) exit
  end do
  call system_clock(t_end, count_rate)

  write(*,"(A,Es10.3,A)") " Wall-clock time generating AMR grid: ", &
       (t_end-t_start) / real(count_rate,dp), " seconds"

  call a2_print_info(tree)

  ! This writes a VTK output file containing the cell-centered values of the
  ! leaves of the tree (the boxes not covered by refinement).
  write(fname, "(A,I0)") "poisson_basic_2d_", 0
  call a2_write_vtk(tree, trim(fname), dir="output")

  ! Initialize the multigrid options. This performs some basics checks and sets
  ! default values where necessary.
  mg%i_phi        = i_phi       ! Solution variable
  mg%i_rhs        = i_rhs       ! Right-hand side variable
  mg%i_tmp        = i_tmp       ! Variable for temporary space
  mg%sides_bc     => sides_bc   ! Method for boundary conditions
  ! This routine does not initialize the multigrid fields boxes%i_phi,
  ! boxes%i_rhs and boxes%i_tmp. These fileds will be initialized at the
  ! first call of mg2_fas_fmg
  call mg2_init_mg(mg)

  print *, "Multigrid iteration | max residual | max error"
  call system_clock(t_start, count_rate)

  do mg_iter = 1, n_iterations
     ! Perform a FAS-FMG cycle (full approximation scheme, full multigrid). The
     ! third argument controls whether the residual is stored in i_tmp. The
     ! fourth argument controls whether to improve the current solution.
     call mg2_fas_fmg(tree, mg, set_residual=.true., have_guess=(mg_iter>1))

     ! Compute the error compared to the analytic solution
     call a2_loop_box(tree, set_error)

     ! Determine the minimum and maximum residual and error
     call a2_tree_min_cc(tree, i_tmp, residu(1))
     call a2_tree_max_cc(tree, i_tmp, residu(2))
     call a2_tree_min_cc(tree, i_err, anal_err(1))
     call a2_tree_max_cc(tree, i_err, anal_err(2))
     write(*,"(I8,13x,2(Es14.5))") mg_iter, maxval(abs(residu)), &
          maxval(abs(anal_err))

     ! This writes a VTK output file containing the cell-centered values of the
     ! leaves of the tree (the boxes not covered by refinement).
     write(fname, "(A,I0)") "poisson_basic_2d_", mg_iter
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

  ! Return the refinement flags for box
  subroutine refine_routine(box, cell_flags)
    type(box2_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(box%n_cell, box%n_cell)
    integer                  :: i, j, nc
    real(dp)                 :: xy(2), dr2, drhs

    nc = box%n_cell
    dr2 = box%dr**2

    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i, j])

          ! This is an estimate of the truncation error in the right-hand side,
          ! which is related to the fourth derivative of the solution.
          drhs = dr2 * box%cc(i, j, i_rhs)

          if (abs(drhs) > 1e-3_dp) then
             cell_flags(i, j) = af_do_ref
          else
             cell_flags(i, j) = af_keep_ref
          end if
       end do
    end do
  end subroutine refine_routine

  ! This routine sets the initial conditions for each box
  subroutine set_initial_condition(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell

    do j = 1, nc
       do i = 1, nc
          ! Get the coordinate of the cell center at i,j
          xy = a2_r_cc(box, [i,j])

          ! And set the rhs values
          box%cc(i, j, i_rhs) = gauss_laplacian(gs, xy)
       end do
    end do
  end subroutine set_initial_condition

  ! Set the error compared to the analytic solution
  subroutine set_error(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, i_err) = box%cc(i, j, i_phi) - gauss_value(gs, xy)
       end do
    end do
  end subroutine set_error

  ! This routine sets boundary conditions for a box, by filling its ghost cells
  ! with approriate values.
  subroutine sides_bc(box, nb, iv, bc_type)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: nb ! Direction for the boundary condition
    integer, intent(in)         :: iv ! Index of variable
    integer, intent(out)        :: bc_type ! Type of boundary condition
    real(dp)                    :: xy(2)
    integer                     :: n, nc

    nc = box%n_cell

    ! We use dirichlet boundary conditions
    bc_type = af_bc_dirichlet

    ! Below the solution is specified in the approriate ghost cells
    select case (nb)
    case (a2_neighb_lowx)             ! Lower-x direction
       do n = 1, nc
          xy = a2_rr_cc(box, [0.5_dp, real(n, dp)])
          box%cc(0, n, iv) = gauss_value(gs, xy)
       end do
    case (a2_neighb_highx)             ! Higher-x direction
       do n = 1, nc
          xy = a2_rr_cc(box, [nc+0.5_dp, real(n, dp)])
          box%cc(nc+1, n, iv) = gauss_value(gs, xy)
       end do
    case (a2_neighb_lowy)             ! Lower-y direction
       do n = 1, nc
          xy = a2_rr_cc(box, [real(n, dp), 0.5_dp])
          box%cc(n, 0, iv) = gauss_value(gs, xy)
       end do
    case (a2_neighb_highy)             ! Higher-y direction
       do n = 1, nc
          xy = a2_rr_cc(box, [real(n, dp), nc+0.5_dp])
          box%cc(n, nc+1, iv) = gauss_value(gs, xy)
       end do
    end select
  end subroutine sides_bc

end program poisson_basic_2d
