#include "../src/cpp_macros.h"
!> \example electrode_dielectric.f90
!>
!> Example showing how to include a dielectric surface
program electrode_dielectric
  use m_af_all

  implicit none

  integer, parameter :: box_size     = 8
  integer, parameter :: n_iterations = 10
  integer            :: i_phi
  integer            :: i_rhs
  integer            :: i_tmp
  integer            :: i_lsf
  integer            :: i_eps
  integer            :: i_field
  integer            :: i_field_norm

  type(af_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: n, mg_iter, coord, n_args
  real(dp)           :: residu
  character(len=100) :: fname
  type(mg_t)         :: mg

  if (NDIM /= 2) error stop "Example only set up for 2D"

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)
  call af_add_cc_variable(tree, "lsf", ix=i_lsf)
  call af_add_cc_variable(tree, "eps", ix=i_eps)
  call af_add_cc_variable(tree, "field_norm", ix=i_field_norm)
  call af_add_fc_variable(tree, "field", ix=i_field)

  call af_set_cc_methods(tree, i_lsf, funcval=set_lsf)
  call af_set_cc_methods(tree, i_rhs, af_bc_neumann_zero)

  ! If an argument is given, switch to cylindrical coordinates in 2D
  n_args = command_argument_count()
  if (NDIM == 2 .and. n_args == 1) then
     coord = af_cyl
  else
     coord = af_xyz
  end if

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [DTIMES(1.0_dp)], &
       4 * [DTIMES(box_size)], &
       coord=coord)

  do n = 1, 4
     call af_adjust_refinement(tree, ref_routine, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  tree%mg_i_lsf = i_lsf
  tree%mg_i_eps = i_eps
  mg%i_phi = i_phi
  mg%i_rhs = i_rhs
  mg%i_tmp = i_tmp
  mg%sides_bc => af_bc_dirichlet_zero
  mg%lsf => get_lsf
  mg%lsf_boundary_value = 1.0_dp

  call mg_init(tree, mg)

  do mg_iter = 1, n_iterations
     call mg_fas_fmg(tree, mg, .true., mg_iter>1)
     call mg_compute_phi_gradient(tree, mg, i_field, 1.0_dp, i_field_norm)

     ! Determine the minimum and maximum residual and error
     call af_tree_maxabs_cc(tree, i_tmp, residu)
     write(*, "(I8,Es14.5)") mg_iter, residu

     write(fname, "(A,I0)") "output/electrode_dielectric_" // &
          DIMNAME // "_", mg_iter
     call af_write_silo(tree, trim(fname))
  end do

contains

  ! Set refinement flags for box
  subroutine ref_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))

    if (box%lvl < 9 - 2 * NDIM .and. box%r_min(2) < 0.5_dp) then
       cell_flags(DTIMES(:)) = af_do_ref
    else
       cell_flags(DTIMES(:)) = af_keep_ref
    end if
  end subroutine ref_routine

  ! This routine sets the level set function
  subroutine set_lsf(box, iv)
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: iv
    real(dp)                   :: rr(NDIM)
    integer                    :: IJK, nc

    nc = box%n_cell

    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, iv) = get_lsf(rr)
       if (rr(1) < 0.5_dp) then
          box%cc(IJK, i_eps) = 20.0_dp
       else
          box%cc(IJK, i_eps) = 1.0_dp
       end if
    end do; CLOSE_DO
  end subroutine set_lsf

  real(dp) function get_lsf(r)
    real(dp), intent(in) :: r(NDIM)
    real(dp)             :: r0(NDIM), r1(NDIM)
    real(dp)             :: radius

    ! Start and end point of line
    r0(:)   = 0.3_dp
    r1(:)   = 0.5_dp
    radius  = 0.02_dp
    get_lsf = dist_vec_line(r, r0, r1, NDIM) - radius
  end function get_lsf

  !> Compute distance vector between point and its projection onto a line
  !> between r0 and r1
  function dist_vec_line(r, r0, r1, n_dim) result(dist)
    integer, intent(in)  :: n_dim
    real(dp), intent(in) :: r(n_dim), r0(n_dim), r1(n_dim)
    real(dp)             :: dist_vec(n_dim), frac
    real(dp)             :: dist, line_len2

    line_len2 = sum((r1 - r0)**2)
    frac = sum((r - r0) * (r1 - r0))

    if (frac <= 0.0_dp) then
       dist_vec = r - r0
    else if (frac >= line_len2) then
       dist_vec = r - r1
    else
       dist_vec = r - (r0 + frac/line_len2 * (r1 - r0))
    end if

    dist = norm2(dist_vec)
  end function dist_vec_line

end program electrode_dielectric
