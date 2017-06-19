program casper_test
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

  print *, "Number of threads", af_get_max_threads()

  ! The cell spacing at the coarsest grid level
  dr = 1.0_dp / box_size

  ! Initialize tree
  call a2_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       n_var_cell, &   ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       coarsen_to=2, & ! Add coarsened levels for multigrid
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
  mg%sides_bc     => sides_bc

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

     ! Determine the minimum and maximum residual and error
     call a2_tree_min_cc(tree, i_tmp, residu(1))
     call a2_tree_max_cc(tree, i_tmp, residu(2))
     write(*,"(I8,Es14.5)") mg_iter, maxval(abs(residu))

     write(fname, "(A,I0)") "casper_test_", mg_iter
     call a2_write_silo(tree, trim(fname), dir="output")
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

    if (box%lvl < 7 .and. box%r_min(1) < 0.5_dp) then
       cell_flags(:,:) = af_do_ref
    else
       cell_flags(:, :) = af_keep_ref
    end if
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
          if (norm2(rz - 0.5_dp) < 0.25_dp) then
             box%cc(i, j, i_eps) = 100.0_dp
          else
             box%cc(i, j, i_eps) = 1.0_dp
          end if

          box%cc(i, j, i_rhs) = 0.0d0
          box%cc(i, j, i_phi) = 0.0d0
       end do
    end do

  end subroutine set_init_cond

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
       bc_type = af_bc_neumann
       box%cc(nc+1, 1:nc, iv) = 0
    case (a2_neighb_lowy)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, 0, iv) = 0
    case (a2_neighb_highy)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, nc+1, iv) = 1
    end select
  end subroutine sides_bc

end program casper_test
