program test_mg
  use m_afivo_2d
  use m_mg2d

  implicit none

  integer, parameter :: dp           = kind(0.0d0)
  integer, parameter :: box_size     = 8
  integer, parameter :: n_boxes_base = 1
  integer, parameter :: i_phi = 1, i_tmp = 2
  integer, parameter :: i_rhs = 3, i_res = 4
  integer, parameter :: i_err = 5

  ! The manufactured solution exists of two Gaussians here.
  ! For each Gaussian, 4 constants are used: pre-factor, x0, y0, sigma.
  integer, parameter :: n_gaussians = 2
  real(dp), parameter :: g_params(4, n_gaussians) = reshape(&
       [1.0_dp, 0.25_dp, 0.25_dp, 0.05_dp, &
       1.0_dp, 0.75_dp, 0.75_dp, 0.05_dp], [4,2])

  type(a2_t)         :: tree
  integer            :: i, id, n_changes
  integer            :: ix_list(2, n_boxes_base)
  integer            :: nb_list(4, n_boxes_base)
  integer            :: n_boxes_max  = 20*1000
  integer            :: n_lvls_max   = 20
  real(dp)           :: dr
  character(len=40)  :: fname, var_names(5)
  type(mg2_t)        :: mg

  var_names(i_phi) = "phi"
  var_names(i_tmp) = "tmp"
  var_names(i_rhs) = "rhs"
  var_names(i_res) = "res"
  var_names(i_err) = "err"

  dr = 1.0_dp / box_size

  ! Initialize tree
  call a2_init(tree, n_lvls_max, n_boxes_max, box_size, n_var_cell=5, &
       n_var_face=0, dr = dr, r_min = [0.0_dp, 0.0_dp])

  id = 1
  ix_list(:, id) = [1,1]         ! Set index of boxnn
  nb_list(:, id) = -1            ! Dirichlet zero -> -1

  call a2_set_base(tree, ix_list, nb_list)

  do i = 1, 20
     call a2_adjust_refinement(tree, ref_func_init, n_changes)
     if (n_changes == 0) exit
  end do

  ! Set rhs and initial guess for phi
  call a2_loop_box(tree, set_init_cond)

  ! Set the multigrid options
  call mg2d_set(mg, i_phi, i_tmp, i_rhs, i_res, 2, 2, 2, &
       sides_bc, a2_corners_extrap, mg2d_lpl_box, mg2d_gsrb_lpl_box)

  ! Create a "subtree" with coarser levels than tree
  call mg2d_create_subtree(tree)

  ! Restrict from children recursively
  call mg2d_restrict_trees(tree, i_rhs, mg)
  call mg2d_restrict_trees(tree, i_phi, mg)

  !$omp parallel
  do i = 1, 10
     ! call mg2d_fas_vcycle(tree, mg, tree%n_lvls)
     call mg2d_fas_fmg(tree, mg)
     call a2_loop_box(tree, set_err)
     write(fname, "(A,I0,A)") "test_mg2_", i, ".vtu"
     call a2_write_tree(tree, trim(fname), var_names, i, 0.0_dp)
  end do
  !$omp end parallel

  print *, "max_id", tree%max_id
  print *, "n_cells", tree%max_id * tree%n_cell**2

  call a2_destroy(tree)

contains

  integer function ref_func_init(box)
    type(box2_t), intent(in) :: box
    integer                  :: n

    ref_func_init = a5_rm_ref

    if (box%lvl < 4) then
       ref_func_init = a5_do_ref
    else if (box%lvl < 8) then
       do n = 1, n_gaussians
          if (norm2(a2_r_center(box) - g_params(2:3, n)) < 0.1_dp) then
             ref_func_init = a5_do_ref
             exit
          end if
       end do
    end if
  end function ref_func_init

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    box%cc(:, :, i_phi) = 0

    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, i_rhs) = rhs(xy)
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
          box%cc(i, j, i_err) = box%cc(i, j, i_phi) - phi_sol(xy)
       end do
    end do
  end subroutine set_err

  subroutine sides_bc(boxes, id, nb, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    real(dp)                    :: xy(2)
    integer                     :: n, nc

    nc = boxes(id)%n_cell

    if (boxes(id)%neighbors(nb) == -1) then
       select case (nb)
       case (a2_nb_lx)
          do n = 1, nc
             xy = a2_rr_cc(boxes(id), [0.5_dp, real(n, dp)])
             boxes(id)%cc(0, n, iv) = 2 * phi_sol(xy) - boxes(id)%cc(1, n, iv)
          end do
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
    end if
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
    real(dp), parameter  :: fac = sqrt(acos(-1.0_dp))
    gaussian_2d = exp(-0.5_dp/sigma**2 * sum((x-x0)**2)) / (sigma*fac)
  end function gaussian_2d

  real(dp) function lpl_gaussian_2d(x, x0, sigma)
    real(dp), intent(in) :: x(2), x0(2), sigma
    lpl_gaussian_2d = 4 * (-0.5_dp/sigma**2 + 0.25_dp/sigma**4 * sum((x-x0)**2)) * &
         gaussian_2d(x, x0, sigma)
  end function lpl_gaussian_2d

end program test_mg
