!> \example test_mg_cyl_diel.f90
!>
!> Example showing how to use m_a2_mg in cylindrical
!> coordinates with a change in eps, and compare with an
!> analytic solution.
program test_mg_cyl_diel
  use m_a2_t
  use m_a2_core
  use m_a2_mg
  use m_a2_utils
  use m_a2_io

  implicit none

  integer, parameter :: box_size     = 8
  integer, parameter :: n_boxes_base = 1
  integer, parameter :: i_phi = 1, i_tmp = 2
  integer, parameter :: i_rhs = 3, i_err = 4
  integer, parameter :: i_eps = 5, i_sol = 6

  ! The manufactured solution exists of two Gaussians here.
  ! For each Gaussian, 4 constants are used: pre-factor, x0, y0, sigma.
  integer, parameter :: n_gaussians = 2
  real(dp), parameter :: g_params(4, n_gaussians) = reshape(&
       [1.0_dp, 0.25_dp, 0.25_dp, 0.04_dp, &
       1.0_dp, 0.75_dp, 0.75_dp, 0.04_dp], [4,2])

  type(a2_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: i, id
  integer            :: ix_list(2, n_boxes_base)
  integer            :: nb_list(4, n_boxes_base)
  real(dp)           :: dr, max_res, min_res
  character(len=40)  :: fname, var_names(6)
  type(mg2_t)        :: mg

  var_names(i_phi) = "phi"
  var_names(i_tmp) = "tmp"
  var_names(i_rhs) = "rhs"
  var_names(i_err) = "err"
  var_names(i_eps) = "eps"
  var_names(i_sol) = "sol"

  dr = 1.0_dp / box_size

  ! Initialize tree
  call a2_init(tree, box_size, n_var_cell=6, n_var_face=0, &
       dr = dr, coarsen_to = 2, coord=a5_cyl, n_boxes=10*1000)

  id = 1
  ix_list(:, id) = [1,1]         ! Set index of boxnn
  nb_list(:, id) = -1            ! Dirichlet zero -> -1

  call a2_set_base(tree, ix_list, nb_list)

  do
     call a2_loop_box(tree, set_init_cond)
     call a2_adjust_refinement(tree, set_ref_flags, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  ! Set rhs and initial guess for phi
  call a2_loop_box(tree, set_init_cond)

  ! Set the multigrid options
  mg%i_phi        = i_phi
  mg%i_tmp        = i_tmp
  mg%i_rhs        = i_rhs
  mg%i_eps        = i_eps
  mg%n_cycle_down = 2
  mg%n_cycle_up   = 2
  mg%n_cycle_base = 2
  mg%sides_bc     => sides_bc
  mg%box_op       => mg2_auto_op
  mg%box_gsrb     => mg2_auto_gsrb
  mg%box_corr     => mg2_auto_corr

  call mg2_init_mg(mg)

  do i = 1, 12
     ! call mg2_fas_vcycle(tree, mg, tree%n_lvls)
     call mg2_fas_fmg(tree, mg, .true., i == 1)
     call a2_loop_box(tree, set_err)
     call a2_tree_min_cc(tree, i_tmp, min_res)
     call a2_tree_max_cc(tree, i_tmp, max_res)
     print *, i, max(abs(min_res), abs(max_res))
     write(fname, "(A,I0)") "test_mg_cyl_diel_", i
     ! call a2_write_silo(tree, trim(fname), var_names, i, 0.0_dp)
  end do

  print *, "max_id", tree%max_id
  print *, "n_cells", tree%max_id * tree%n_cell**2

  call a2_destroy(tree)

contains

  subroutine set_ref_flags(boxes, id, ref_flags)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flags(:)
    integer                  :: nc
    real(dp)                 :: max_crv

    nc = boxes(id)%n_cell
    max_crv = boxes(id)%dr**2 * maxval(abs(boxes(id)%cc(1:nc, 1:nc, i_rhs) / &
         boxes(id)%cc(1:nc, 1:nc, i_eps)))

    if (max_crv > 5.0e-4_dp) then
       ref_flags(id) = a5_do_ref
    end if
  end subroutine set_ref_flags

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2), grad(2), dr, qbnd, tmp

    nc                  = box%n_cell
    box%cc(:, :, i_phi) = 0
    dr                  = box%dr

    do j = 0, nc+1
       do i = 0, nc+1
          xy = a2_r_cc(box, [i,j])

          if (xy(1) < 0.5_dp .and. xy(2) < 0.5_dp) then
             box%cc(i, j, i_eps) = 100.0_dp
          else
             box%cc(i, j, i_eps) = 1.0_dp
          end if
          box%cc(i, j, i_rhs) = rhs(xy) * box%cc(i, j, i_eps)
          box%cc(i, j, i_sol) = phi_sol(xy)
       end do
    end do

    ! Set surface charge in r-direction
    do j = 1, nc
       do i = 0, nc
          xy = a2_rr_cc(box, [i + 0.5_dp, real(j, dp)])
          grad = phi_grad(xy)
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
          xy = a2_rr_cc(box, [real(i, dp), j + 0.5_dp])
          grad = phi_grad(xy)
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

  subroutine set_err(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, i_err) = box%cc(i, j, i_phi) - box%cc(i, j, i_sol)
       end do
    end do
  end subroutine set_err

  subroutine sides_bc(boxes, id, nb, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    real(dp)                    :: xy(2)
    integer                     :: n, nc

    nc = boxes(id)%n_cell

    select case (nb)
    case (a2_nb_lx)
       boxes(id)%cc(0, 1:nc, iv) = boxes(id)%cc(1, 1:nc, iv)
    case (a2_nb_hx)
       do n = 1, nc
          xy = a2_rr_cc(boxes(id), [nc+0.5_dp, real(n, dp)])
          boxes(id)%cc(nc+1, n, iv) = 2 * phi_sol(xy) - boxes(id)%cc(nc, n, iv)
       end do
    case (a2_nb_ly)
       do n = 1, nc
          xy = a2_rr_cc(boxes(id), [real(n, dp), 0.5_dp])
          boxes(id)%cc(n, 0, iv) = 2 * phi_sol(xy) - boxes(id)%cc(n, 1, iv)
       end do
    case (a2_nb_hy)
       do n = 1, nc
          xy = a2_rr_cc(boxes(id), [real(n, dp), nc+0.5_dp])
          boxes(id)%cc(n, nc+1, iv) = 2 * phi_sol(xy) - boxes(id)%cc(n, nc, iv)
       end do
    end select
  end subroutine sides_bc

  real(dp) function phi_sol(x)
    real(dp), intent(in) :: x(2)
    integer :: n

    phi_sol = 0
    do n = 1, n_gaussians
       phi_sol = phi_sol + g_params(1, n) * &
            gaussian_2d(x, g_params(2:3, n), g_params(4, n))
    end do
  end function phi_sol

  function phi_grad(x) result(grad)
    real(dp), intent(in) :: x(2)
    real(dp) :: grad(2)
    integer :: n

    grad = 0
    do n = 1, n_gaussians
       grad = grad + g_params(1, n) * &
            grad_gaussian_2d(x, g_params(2:3, n), g_params(4, n))
    end do
  end function phi_grad

  real(dp) function rhs(x)
    real(dp), intent(in) :: x(2)
    integer :: n

    rhs = 0
    do n = 1, n_gaussians
       rhs = rhs + g_params(1, n) * &
            lpl_gaussian_2d(x, g_params(2:3, n), g_params(4, n))
    end do
  end function rhs

  real(dp) function gaussian_2d(x, x0, sigma)
    real(dp), intent(in) :: x(2), x0(2), sigma
    real(dp) :: xrel(2)
    xrel = (x-x0)/sigma
    gaussian_2d = exp(-sum(xrel**2))
  end function gaussian_2d

  function grad_gaussian_2d(x, x0, sigma) result(grad)
    real(dp), intent(in) :: x(2), x0(2), sigma
    real(dp) :: xrel(2), grad(2)
    xrel = (x-x0)/sigma
    grad = -2 * xrel/sigma * gaussian_2d(x, x0, sigma)
  end function grad_gaussian_2d

  real(dp) function lpl_gaussian_2d(x, x0, sigma)
    real(dp), intent(in) :: x(2), x0(2), sigma
    real(dp) :: xrel(2)
    xrel = (x-x0)/sigma

    lpl_gaussian_2d = 4/sigma**2 * (sum(xrel**2) - 1.0_dp - &
         0.5_dp * (x(1)-x0(1))/x(1)) * gaussian_2d(x, x0, sigma)
  end function lpl_gaussian_2d

end program test_mg_cyl_diel
