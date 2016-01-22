!> \example poisson_basic_2d.f90
!>
!> Example showing how to use multigrid and compare with an analytic solution. A
!> standard 5-point Laplacian is used.
program poisson_basic_2d
  use m_a2_t
  use m_a2_core
  use m_a2_mg
  use m_a2_utils
  use m_a2_io

  implicit none

  integer, parameter :: n_cell = 8
  integer, parameter :: n_boxes_base = 1
  integer, parameter :: n_var_cell = 4
  integer, parameter :: i_phi = 1
  integer, parameter :: i_rhs = 2
  integer, parameter :: i_err = 3
  integer, parameter :: i_tmp = 4

  ! The manufactured solution exists of two Gaussians here. For each Gaussian, 4
  ! constants are used: pre-factor, x0, y0, sigma.
  integer, parameter :: n_gaussians = 2
  real(dp), parameter :: g_params(4, n_gaussians) = reshape(&
       [1.0_dp, 0.25_dp, 0.25_dp, 0.04_dp, &
       1.0_dp, 0.75_dp, 0.75_dp, 0.04_dp], [4,2])

  type(a2_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: i
  integer            :: ix_list(2, n_boxes_base)
  integer            :: nb_list(4, n_boxes_base)
  real(dp)           :: dr, min_res, max_res
  character(len=40)  :: fname
  type(mg2_t)        :: mg

  ! The cell spacing at the coarsest grid level
  dr = 1.0_dp / n_cell

  ! Initialize tree
  call a2_init(tree, & ! Tree to initialize
       n_cell, &       ! A box contains n_cell**DIM cells
       n_var_cell, &   ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       coarsen_to=2, & ! Add coarsened levels for multigrid
       cc_names=["phi", "rhs", "err", "tmp"]) ! Variable names

  ! Set up geometry. These indices are used to define the coordinates of a box,
  ! by default the box at [1,1] touches the origin (x,y) = (0,0)
  ix_list(:, 1) = [1,1]         ! Set index of box 1

  ! Set neighbors for box one, negative values indicate a physical boundary
  nb_list(:, 1) = -1

  ! Create the base mesh, using the box indices and their neighbor information
  call a2_set_base(tree, ix_list, nb_list)

  do
     ! For each box, set the initial conditions
     call a2_loop_box(tree, set_init_cond)

     ! This updates the refinement of the tree, by at most one level per call.
     call a2_adjust_refinement(tree, set_ref_flags, ref_info)

     ! If no new boxes have been added, exit the loop
     if (ref_info%n_add == 0) exit
  end do

  ! Set the multigrid options.
  mg%i_phi        = i_phi       ! Solution variable
  mg%i_rhs        = i_rhs       ! Right-hand side variable
  mg%i_tmp        = i_tmp       ! Variable for temporary space
  mg%sides_bc     => sides_bc   ! Method for boundary conditions

  ! Initialize the multigrid options. This performs some basics checks and sets
  ! default values where necessary.
  call mg2_init_mg(mg)

  do i = 1, 12
     ! Perform a FAS-FMG cycle (full approximation scheme, full multigrid). The
     ! third argument controls whether the residual is stored in i_tmp. The
     ! fourth argument controls whether to improve the current solution.
     call mg2_fas_fmg(tree, mg, set_residual=.true., have_guess=(i>1))

     ! Compute the error compared to the analytic solution
     call a2_loop_box(tree, set_err)

     ! Determine the minimum and maximum residual
     call a2_tree_min_cc(tree, i_tmp, min_res)
     call a2_tree_max_cc(tree, i_tmp, max_res)
     print *, "Iteration ", i, "max residual: ", max(abs(min_res), abs(max_res))

     ! This writes a Silo output file containing the cell-centered values of the
     ! leaves of the tree (the boxes not covered by refinement).
     write(fname, "(A,I0)") "poisson_basic_2d_", i
     call a2_write_silo(tree, trim(fname), dir="output")
  end do

  call a2_destroy(tree)

contains

  ! This routine sets refinement flags
  subroutine set_ref_flags(boxes, id, ref_flags)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flags(:)
    integer                  :: i, j, nc
    real(dp)                 :: max_crv, xy(2), dr2, derror

    nc = boxes(id)%n_cell
    ! dr2 = boxes(id)%dr**2

    ! do j = 1, nc
    !    do i = 1, nc
    !       xy = a2_r_cc(boxes(id), [i, j])
    !       derror = dr2 * discr_error(xy)
    !    end do
    ! end do

    ! Compute the "curvature" in phi
    max_crv = boxes(id)%dr**2 * &
         maxval(abs(boxes(id)%cc(1:nc, 1:nc, i_rhs)))

    ! ! And refine if it exceeds a threshold
    ! if (max_crv > 10.0e-4_dp) then
    !    ref_flags(id) = a5_do_ref
    ! end if

    if (boxes(id)%lvl < 5 .and. max_crv > 0.0e-4_dp) ref_flags(id) = a5_do_ref
  end subroutine set_ref_flags

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2), dr2

    nc = box%n_cell
    dr2 = box%dr**2

    do j = 1, nc
       do i = 1, nc
          ! Get the coordinate of the cell center at i,j
          xy = a2_r_cc(box, [i,j])

          ! And set the rhs values
          box%cc(i, j, i_rhs) = analytic_rhs(xy)! + &
               ! dr2 * analytic_fourth(xy) / 12
       end do
    end do
  end subroutine set_init_cond

  ! Set the error compared to the analytic solution
  subroutine set_err(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, i_err) = box%cc(i, j, i_phi) - analytic_solution(xy)
       end do
    end do
  end subroutine set_err

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
    bc_type = a5_bc_dirichlet

    ! Below the solution is specified in the approriate ghost cells
    select case (nb)
    case (a2_nb_lx)             ! Lower-x direction
       do n = 1, nc
          xy = a2_rr_cc(box, [0.5_dp, real(n, dp)])
          box%cc(0, n, iv) = analytic_solution(xy)
       end do
    case (a2_nb_hx)             ! Higher-x direction
       do n = 1, nc
          xy = a2_rr_cc(box, [nc+0.5_dp, real(n, dp)])
          box%cc(nc+1, n, iv) = analytic_solution(xy)
       end do
    case (a2_nb_ly)             ! Lower-y direction
       do n = 1, nc
          xy = a2_rr_cc(box, [real(n, dp), 0.5_dp])
          box%cc(n, 0, iv) = analytic_solution(xy)
       end do
    case (a2_nb_hy)             ! Higher-y direction
       do n = 1, nc
          xy = a2_rr_cc(box, [real(n, dp), nc+0.5_dp])
          box%cc(n, nc+1, iv) = analytic_solution(xy)
       end do
    end select
  end subroutine sides_bc

  ! Analytic solution to the Poisson problem
  real(dp) function analytic_solution(x)
    real(dp), intent(in) :: x(2)
    integer :: n

    analytic_solution = 0
    do n = 1, n_gaussians
       analytic_solution = analytic_solution + g_params(1, n) * &
            gaussian_2d(x, g_params(2:3, n), g_params(4, n))
    end do
  end function analytic_solution

  ! Analytic right-hand side to the Poisson problem
  real(dp) function analytic_rhs(x)
    real(dp), intent(in) :: x(2)
    integer              :: n

    analytic_rhs = 0
    do n = 1, n_gaussians
       analytic_rhs = analytic_rhs + g_params(1, n) * &
            lpl_gaussian_2d(x, g_params(2:3, n), g_params(4, n))
    end do
  end function analytic_rhs

  ! A Gaussian in xy coordinates
  real(dp) function gaussian_2d(x, x0, sigma)
    real(dp), intent(in) :: x(2), x0(2), sigma
    real(dp)             :: xrel(2)

    xrel = (x-x0)/sigma
    gaussian_2d = exp(-sum(xrel**2))
  end function gaussian_2d

  ! Laplacian of a Gaussian in xy coordinates
  real(dp) function lpl_gaussian_2d(x, x0, sigma)
    real(dp), intent(in) :: x(2), x0(2), sigma
    real(dp)             :: xrel(2)

    xrel = (x-x0)/sigma
    lpl_gaussian_2d = 4/sigma**2 * (sum(xrel**2) - 1) * &
         gaussian_2d(x, x0, sigma)
  end function lpl_gaussian_2d

  ! Fourth derivative of a Gaussian
  real(dp) function d4_gaussian_2d(x, x0, sigma)
    real(dp), intent(in) :: x(2), x0(2), sigma
    real(dp)             :: xrel(2)

    xrel = (x-x0)/sigma
    d4_gaussian_2d = gaussian_2d(x, x0, sigma) / sigma**4 * ( &
         16 * sum(xrel**4) - 48 * sum(xrel**2) + 24)
  end function d4_gaussian_2d

  ! Approximate discretization error of a Gaussian
  real(dp) function analytic_fourth(x)
    real(dp), intent(in) :: x(2)
    integer              :: n

    analytic_fourth = 0
    do n = 1, n_gaussians
       analytic_fourth = analytic_fourth + g_params(1, n) * &
            d4_gaussian_2d(x, g_params(2:3, n), g_params(4, n))
    end do
  end function analytic_fourth

end program
