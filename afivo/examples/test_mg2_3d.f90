!> \example test_mg2_3d.f90
!> Example showing how to use m_mg_3d, and compare with an analytic solution.
program test_mg2_3d
  use m_afivo_3d
  use m_mg_3d

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
  real(dp), parameter :: g_params(5, n_gaussians) = reshape(&
       [1.0_dp, 0.25_dp, 0.25_dp, 0.25_dp, 0.05_dp, &
       1.0_dp, 0.75_dp, 0.75_dp, 0.75_dp, 0.05_dp], [5,2])

  type(a3_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: i, id
  integer            :: ix_list(3, n_boxes_base)
  integer            :: nb_list(6, n_boxes_base)
  real(dp)           :: dr
  character(len=40)  :: fname, var_names(5)
  type(mg3_t)        :: mg

  var_names(i_phi) = "phi"
  var_names(i_tmp) = "tmp"
  var_names(i_rhs) = "rhs"
  var_names(i_res) = "res"
  var_names(i_err) = "err"

  dr = 1.0_dp / box_size

  ! Initialize tree
  call a3_init(tree, box_size, n_var_cell=5, &
       n_var_face=0, dr = dr, coarsen_to=2)

  id = 1
  ix_list(:, id) = [1,1,1]       ! Set index of boxnn
  nb_list(:, id) = -1            ! Dirichlet zero -> -1

  print *, "Setting up base"
  call a3_set_base(tree, ix_list, nb_list)

  print *, "Doing initial refinement"
  do i = 1, 20
     call a3_adjust_refinement(tree, set_ref_flags, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  print *, "Set rhs and initial guess for phi"
  call a3_loop_box(tree, set_init_cond)

  print *, "Set the multigrid options"
  mg%i_phi        = i_phi
  mg%i_tmp        = i_tmp
  mg%i_rhs        = i_rhs
  mg%i_res        = i_res
  mg%n_cycle_down = 2
  mg%n_cycle_up   = 2
  mg%n_cycle_base = 2
  mg%sides_bc     => sides_bc

  call mg3_init_mg(mg)

  print *, "Restrict from children recursively"
  call a3_restrict_tree(tree, i_rhs)
  call a3_restrict_tree(tree, i_phi)

  print *, "Do multigrid"
  do i = 1, 10
     ! call mg3_fas_vcycle(tree, mg, tree%n_lvls)
     call mg3_fas_fmg(tree, mg)
     call a3_loop_box(tree, set_err)
     write(fname, "(A,I0,A)") "test_mg2_3d_", i, ".vtu"
     call a3_write_vtk(tree, trim(fname), var_names, i, 0.0_dp)
  end do

  print *, "max_id", tree%max_id
  print *, "n_cells", tree%max_id * tree%n_cell**3

  call a3_destroy(tree)

contains

  subroutine set_ref_flags(boxes, id, ref_flags)
    type(box3_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flags(:)
    integer                  :: n

    ref_flags(id) = a5_rm_ref

    if (boxes(id)%lvl < 3) then
       ref_flags(id) = a5_do_ref
    else if (boxes(id)%lvl < 5) then
       do n = 1, n_gaussians
          if (norm2(a3_r_center(boxes(id)) - g_params(2:4, n)) < 0.25_dp) then
             ref_flags(id) = a5_do_ref
             exit
          end if
       end do
    end if
  end subroutine set_ref_flags

  subroutine set_init_cond(box)
    type(box3_t), intent(inout) :: box
    integer                     :: i, j, k, nc
    real(dp)                    :: xyz(3)

    nc = box%n_cell
    box%cc(:, :, :, i_phi) = 0

    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             xyz = a3_r_cc(box, [i,j,k])
             box%cc(i, j, k, i_rhs) = rhs(xyz)
          end do
       end do
    end do
  end subroutine set_init_cond

  subroutine set_err(box)
    type(box3_t), intent(inout) :: box
    integer                     :: i, j, k, nc
    real(dp)                    :: xyz(3)

    nc = box%n_cell
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             xyz = a3_r_cc(box, [i,j,k])
             box%cc(i, j, k, i_err) = box%cc(i, j, k, i_phi) - phi_sol(xyz)
          end do
       end do
    end do
  end subroutine set_err

  subroutine sides_bc(boxes, id, nb, iv)
    type(box3_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    real(dp)                    :: xyz(3), loc
    integer                     :: ix, dix, i, j, k, nc

    nc = boxes(id)%n_cell

    if (a3_nb_low(nb)) then
       ix = 0
       dix = 1
    else
       ix = nc+1
       dix = -1
    end if

    loc = ix + 0.5_dp * dix

    select case (a3_nb_dim(nb))
    case (1)
       do k = 1, nc
          do j = 1, nc
             xyz = a3_rr_cc(boxes(id), [loc, real(j, dp), real(k, dp)])
             boxes(id)%cc(ix, j, k, iv) = 2 * phi_sol(xyz) - boxes(id)%cc(ix+dix, j, k, iv)
          end do
       end do
    case (2)
       do k = 1, nc
          do i = 1, nc
             xyz = a3_rr_cc(boxes(id), [real(i, dp), loc, real(k, dp)])
             boxes(id)%cc(i, ix, k, iv) = 2 * phi_sol(xyz) - boxes(id)%cc(i, ix+dix, k, iv)
          end do
       end do
    case (3)
       do j = 1, nc
          do i = 1, nc
             xyz = a3_rr_cc(boxes(id), [real(i, dp), real(j, dp), loc])
             boxes(id)%cc(i, j, ix, iv) = 2 * phi_sol(xyz) - boxes(id)%cc(i, j, ix+dix, iv)
          end do
       end do
    end select
  end subroutine sides_bc

  real(dp) function phi_sol(x)
    real(dp), intent(in) :: x(3)
    integer :: n

    phi_sol = 0
    do n = 1, n_gaussians
       phi_sol = phi_sol + g_params(1, n) * &
            gaussian_3d(x, g_params(2:4, n), g_params(5, n))
    end do
  end function phi_sol

  real(dp) function rhs(x)
    real(dp), intent(in) :: x(3)
    integer :: n

    rhs = 0
    do n = 1, n_gaussians
       rhs = rhs + g_params(1, n) * &
            lpl_gaussian_3d(x, g_params(2:4, n), g_params(5, n))
    end do
  end function rhs

  real(dp) function gaussian_3d(x, x0, sigma)
    real(dp), intent(in) :: x(3), x0(3), sigma
    real(dp) :: xrel(3)
    xrel = (x-x0)/sigma
    gaussian_3d = exp(-sum(xrel**2))
  end function gaussian_3d

  real(dp) function lpl_gaussian_3d(x, x0, sigma)
    real(dp), intent(in) :: x(3), x0(3), sigma
    real(dp) :: xrel(3)
    xrel = (x-x0)/sigma
    lpl_gaussian_3d = 4/sigma**2 * (sum(xrel**2) - 1.5_dp) * &
         gaussian_3d(x, x0, sigma)
  end function lpl_gaussian_3d

end program
