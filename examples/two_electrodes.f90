#include "../src/cpp_macros.h"
!> \example two_electrodes.f90
!>
!> Test with two electrodes
program two_electrodes
  use m_af_all
  use m_config
  use omp_lib

  implicit none

  integer            :: box_size         = 8
  integer            :: n_iterations     = 10
  integer            :: max_refine_level = 7
  integer            :: min_refine_level = 5
  integer            :: i_phi
  integer            :: i_rhs
  integer            :: i_tmp
  integer            :: i_lsf
  integer            :: i_field
  integer            :: i_field_norm

  type(af_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: n, mg_iter
  real(dp)           :: residu
  character(len=150) :: fname
  type(mg_t)         :: mg
  real(dp)           :: mem_limit_gb = 8.0_dp
  real(dp)           :: t0, t_sum

  real(dp) :: rod_r0(NDIM), rod_r1(NDIM), rod_radius
  real(dp) :: sphere_r0(NDIM), sphere_radius

  rod_r0 = 0.5_dp
  rod_r0(NDIM) = 0.7_dp
  rod_r1 = 0.5_dp
  rod_r1(NDIM) = 1.0_dp
  rod_radius = 0.05_dp

  sphere_r0 = 0.5_dp
  sphere_r0(NDIM) = 0.0_dp
  sphere_radius = 0.25_dp

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)
  call af_add_cc_variable(tree, "lsf", ix=i_lsf)
  call af_add_cc_variable(tree, "field_norm", ix=i_field_norm)
  call af_add_fc_variable(tree, "field", ix=i_field)

  call af_set_cc_methods(tree, i_lsf, funcval=set_lsf)

  tree%mg_i_lsf = i_lsf
  mg%i_phi = i_phi
  mg%i_rhs = i_rhs
  mg%i_tmp = i_tmp
  mg%sides_bc => af_bc_dirichlet_zero
  mg%lsf_boundary_value = 1.0_dp
  mg%lsf => get_lsf
  mg%lsf_length_scale = 1e-2_dp
  mg%lsf_boundary_function => get_boundary_value
  mg%sides_bc => bc_sides

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [DTIMES(1.0_dp)], &
       [DTIMES(box_size)], &
       mem_limit_gb=mem_limit_gb)

  call af_refine_up_to_lvl(tree, min_refine_level)
  call mg_init(tree, mg)

  do n = 1, 100
     call af_adjust_refinement(tree, ref_routine, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  call af_print_info(tree)

  t_sum = 0
  write(*, "(A,A4,4A14)") "# ", "iter", "residu"

  do mg_iter = 1, n_iterations
     t0 = omp_get_wtime()
     call mg_fas_fmg(tree, mg, .true., mg_iter>1)
     t_sum = t_sum + omp_get_wtime() - t0

     call mg_compute_phi_gradient(tree, mg, i_field, 1.0_dp, i_field_norm)
     call af_loop_box(tree, set_field_to_zero_inside)

     ! Determine the minimum and maximum residual and error
     call af_tree_maxabs_cc(tree, i_tmp, residu)
     write(*, "(A,i4,4E14.5)") "# ", mg_iter, residu

     write(fname, "(A,I0)") "output/two_electrodes_" // DIMNAME // "_", mg_iter
     call af_write_silo(tree, trim(fname))
  end do

  write(*, "(A,E14.5)") " seconds_per_FMG_cycle", t_sum/n_iterations

  call af_stencil_print_info(tree)

contains

  ! Set refinement flags for box
  subroutine ref_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))
    integer                 :: IJK, nc
    real(dp)                :: r(NDIM), r_ref(NDIM)

    nc = box%n_cell
    cell_flags(DTIMES(:)) = af_keep_ref

    ! Refine around this point
    r_ref = rod_r0
    r_ref(NDIM) = r_ref(NDIM) - rod_radius

    if (box%lvl < max_refine_level) then
       do KJI_DO(1, nc)
          r = af_r_cc(box, [IJK])
          if (norm2(r - r_ref) < 2 * rod_radius) then
             cell_flags(DTIMES(:)) = af_do_ref
          end if
       end do; CLOSE_DO
    end if
  end subroutine ref_routine

  ! This routine sets the level set function
  subroutine set_lsf(box, iv)
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: iv
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM)

    nc = box%n_cell

    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, iv) = get_lsf(rr)
    end do; CLOSE_DO
  end subroutine set_lsf

  ! The value at the lsf boundaries
  real(dp) function get_boundary_value(r)
    real(dp), intent(in) :: r(NDIM)

    if (r(NDIM) > 0.5_dp) then
       get_boundary_value = 1.0_dp
    else
       get_boundary_value = 0.0_dp
    end if
  end function get_boundary_value

  real(dp) function get_lsf(r)
    real(dp), intent(in) :: r(NDIM)

    if (r(NDIM) > 0.5_dp) then
       get_lsf = GM_dist_line(r, rod_r0, rod_r1, NDIM) - rod_radius
    else
       get_lsf = norm2(r - sphere_r0) - sphere_radius
    end if
  end function get_lsf

  !> To set boundary conditions at the sides of the domain
  subroutine bc_sides(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    bc_type = af_bc_dirichlet

    if (af_neighb_dim(nb) == NDIM) then
       if (af_neighb_low(nb)) then
          bc_type = af_bc_dirichlet
          bc_val = 0.0_dp
       else
          bc_type = af_bc_dirichlet
          bc_val  = 1.0_dp
       end if
    else
       bc_type = af_bc_neumann
       bc_val = 0.0_dp
    end if
  end subroutine bc_sides

  !> To set the electric field to zero inside the object
  subroutine set_field_to_zero_inside(box)
    type(box_t), intent(inout) :: box

    where (box%cc(DTIMES(:), i_lsf) < 0)
       box%cc(DTIMES(:), i_field_norm) = 0
    end where
  end subroutine set_field_to_zero_inside

  !> Compute distance vector between point and its projection onto a line
  !> between r0 and r1
  subroutine GM_dist_vec_line(r, r0, r1, n_dim, dist_vec, frac)
    integer, intent(in)   :: n_dim
    real(dp), intent(in)  :: r(n_dim), r0(n_dim), r1(n_dim)
    real(dp), intent(out) :: dist_vec(n_dim)
    real(dp), intent(out) :: frac !< Fraction [0,1] along line
    real(dp)              :: line_len2

    line_len2 = sum((r1 - r0)**2)
    frac = sum((r - r0) * (r1 - r0))

    if (frac <= 0.0_dp) then
       frac = 0.0_dp
       dist_vec = r - r0
    else if (frac >= line_len2) then
       frac = 1.0_dp
       dist_vec = r - r1
    else
       dist_vec = r - (r0 + frac/line_len2 * (r1 - r0))
       frac = sqrt(frac / line_len2)
    end if
  end subroutine GM_dist_vec_line

  function GM_dist_line(r, r0, r1, n_dim) result(dist)
    integer, intent(in)  :: n_dim
    real(dp), intent(in) :: r(n_dim), r0(n_dim), r1(n_dim)
    real(dp)             :: dist, dist_vec(n_dim), frac
    call GM_dist_vec_line(r, r0, r1, n_dim, dist_vec, frac)
    dist = norm2(dist_vec)
  end function GM_dist_line

end program
