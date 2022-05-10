#include "../src/cpp_macros.h"
!> \example poisson_lsf_test.f90
!>
!> Test of the level-set functionality of the Poisson solver
program poisson_lsf_test
  use m_af_all
  use m_config

  implicit none

  integer            :: box_size         = 8
  integer, parameter :: n_iterations     = 10
  integer            :: max_refine_level = 4
  integer            :: min_refine_level = 2
  integer            :: i_phi
  integer            :: i_rhs
  integer            :: i_tmp
  integer            :: i_lsf
  integer            :: i_lsf_mask
  integer            :: i_error
  integer            :: i_field
  integer            :: i_field_norm

  logical :: cylindrical = .false.
  logical :: write_numpy = .false.
  integer :: refinement_type = 1

  ! Which shape to use, 1 = circle, 2 = heart, 3 = rhombus, 4-5 = triangle
  integer             :: shape             = 1

  real(dp) :: sharpness_t             = 4.0_dp
  real(dp) :: boundary_value          = 1.0_dp
  real(dp) :: solution_coeff          = 1.0_dp
  real(dp) :: solution_radius         = 0.25_dp
  real(dp) :: solution_r0(NDIM)       = 0.5_dp
  real(dp) :: solution_smallest_width = 1e-3_dp

  type(af_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: n, mg_iter, coord
  real(dp)           :: residu, max_error, max_field
  character(len=150) :: fname
  character(len=20)  :: output_suffix = ""
  type(mg_t)         :: mg
  type(CFG_t)        :: cfg
  real(dp)           :: error_squared, rmse
  logical            :: write_output = .true.

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)
  call af_add_cc_variable(tree, "lsf", ix=i_lsf)
  call af_add_cc_variable(tree, "lsf_mask", ix=i_lsf_mask)
  call af_add_cc_variable(tree, "error", ix=i_error)
  call af_add_cc_variable(tree, "field_norm", ix=i_field_norm)
  call af_add_fc_variable(tree, "field", ix=i_field)

  call af_set_cc_methods(tree, i_lsf, funcval=set_lsf)
  call af_set_cc_methods(tree, i_lsf_mask, af_bc_neumann_zero)
  call af_set_cc_methods(tree, i_field_norm, af_bc_neumann_zero)

  call CFG_update_from_arguments(cfg)

  call CFG_add_get(cfg, "cylindrical", cylindrical, &
       "Whether to use cylindrical coordinates")
  call CFG_add_get(cfg, "write_output", write_output, &
       "Whether to write Silo output")
  call CFG_add_get(cfg, "write_numpy", write_numpy, &
       "Whether to write Numpy output")
  call CFG_add_get(cfg, "output_suffix", output_suffix, &
       "Suffix to add to output")
  call CFG_add_get(cfg, "max_refine_level", max_refine_level, &
       "Maximum refinement level")
  call CFG_add_get(cfg, "min_refine_level", min_refine_level, &
       "Minimum refinement level")
  call CFG_add_get(cfg, "refinement_type", refinement_type, &
       "Type of refinement (1: uniform, 2: only boundary)")
  call CFG_add_get(cfg, "shape", shape, &
       "Index of the shape used")
  call CFG_add_get(cfg, "box_size", box_size, &
       "Size of grid boxes")
  call CFG_add_get(cfg, "solution%sharpness", sharpness_t, &
       "Sharpness of solution (for some cases)")
  call CFG_add_get(cfg, "solution%r0", solution_r0, &
       "Origin of solution")
  call CFG_add_get(cfg, "solution%radius", solution_radius, &
       "Solution radius")
  call CFG_add_get(cfg, "solution%coeff", solution_coeff, &
       "Linear coefficient for solution")
  call CFG_add_get(cfg, "solution%smallest_width", solution_smallest_width, &
       "Smallest width to resolve in the solution")
  call CFG_add_get(cfg, "boundary_value", boundary_value, &
       "Value for Dirichlet boundary condition")

  call CFG_check(cfg)
  ! call CFG_write(cfg, "stdout")

  if (cylindrical) then
     coord = af_cyl
     ! Place solution on axis
     solution_r0(1) = 0.0_dp
  else
     coord = af_xyz
  end if

  tree%mg_i_lsf = i_lsf
  mg%i_phi = i_phi
  mg%i_rhs = i_rhs
  mg%i_tmp = i_tmp
  mg%sides_bc => bc_solution
  mg%lsf_boundary_value = boundary_value
  mg%lsf_dist => mg_lsf_dist_gss
  mg%lsf => get_lsf
  mg%lsf_length_scale = solution_smallest_width

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [DTIMES(1.0_dp)], &
       [DTIMES(box_size)], &
       coord=coord)

  call af_refine_up_to_lvl(tree, min_refine_level)
  call mg_init(tree, mg)

  do n = 1, 100
     call mg_set_operators_tree(tree, mg)
     call af_adjust_refinement(tree, ref_routine, ref_info, ref_buffer=1)
     if (ref_info%n_add == 0) exit
  end do

  call af_print_info(tree)

  write(*, "(A,A4,4A14)") "# ", "iter", "residu", "max_error", "rmse", "max_field"

  do mg_iter = 1, n_iterations
     call mg_fas_fmg(tree, mg, .true., mg_iter>1)
     call mg_compute_phi_gradient(tree, mg, i_field, 1.0_dp, i_field_norm)
     call af_gc_tree(tree, [i_field_norm])
     call af_loop_box(tree, set_error)

     ! Set mask where lsf changes sign
     call af_loop_box(tree, set_lsf_mask)
     call af_gc_tree(tree, [i_lsf_mask])

     ! Determine the minimum and maximum residual and error
     call af_tree_maxabs_cc(tree, i_tmp, residu)
     call af_tree_maxabs_cc(tree, i_error, max_error)
     call af_tree_maxabs_cc(tree, i_field_norm, max_field)
     call af_tree_sum_cc(tree, i_error, error_squared, 2)
     rmse = sqrt(error_squared/af_total_volume(tree))
     write(*, "(A,i4,4E14.5)") "# ", mg_iter, residu, max_error, rmse, max_field

     write(fname, "(A,I0,A)") "output/poisson_lsf_test_" // DIMNAME // &
          "_", mg_iter, trim(output_suffix)
     if (write_output) then
        call af_write_silo(tree, trim(fname))
     end if

     if (write_numpy) then
        call af_write_numpy(tree, trim(fname)//".npz", &
             ixs_cc=[i_phi, i_lsf_mask])
     end if
  end do

  call af_stencil_print_info(tree)

contains

  ! Set refinement flags for box
  subroutine ref_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))
    real(dp)                :: tmp(DTIMES(box%n_cell)), dr_min
    integer                 :: nc, IJK, ix_mask, n, cix(NDIM)

    nc = box%n_cell
    cell_flags(DTIMES(:)) = af_keep_ref

    select case (refinement_type)
    case (1)
       ! Uniform refinement
       if (box%lvl < max_refine_level) then
          cell_flags(DTIMES(:)) = af_do_ref
       end if

    case (2)
       ! Only refine at boundary
       ix_mask = af_stencil_index(box, mg_lsf_mask_key)
       if (ix_mask /= af_stencil_none .and. box%lvl < max_refine_level) then
          do n = 1, size(box%stencils(ix_mask)%sparse_ix, 2)
             cix = box%stencils(ix_mask)%sparse_ix(:, n)
             cell_flags(DINDEX(cix)) = af_do_ref
          end do
       end if

    case (3)
       ! 'Bad' refinement to test the method
       if (norm2(box%r_min - solution_r0) < solution_radius .and. &
            box%lvl < max_refine_level) then
          cell_flags(DTIMES(:)) = af_do_ref
       end if
    case (4)
       ! Let refinement depend on radius
       dr_min = minval(af_lvl_dr(tree, max_refine_level))

       do KJI_DO(1, nc)
          tmp(IJK) = max(1.0_dp, norm2(af_r_cc(box, [IJK]) - solution_r0) &
               / solution_radius)
       end do; CLOSE_DO

       if (minval(box%dr) > minval(tmp) * dr_min .and. box%lvl < max_refine_level) then
          cell_flags(DTIMES(:)) = af_do_ref
       end if
    case default
       error stop "Invalid refinement_type"
    end select
  end subroutine ref_routine

  real(dp) function solution(r)
    real(dp), intent(in) :: r(NDIM)
    real(dp) :: distance, lsf

    select case (shape)
    case (1)
       ! Relative distance
       distance = norm2(r-solution_r0) / solution_radius

       ! Let values increase for distance -> infinity
       if (distance < 1.0_dp) then
          solution = boundary_value
       else if (NDIM == 1) then
          solution = boundary_value + solution_coeff * (distance - 1.0_dp)
       else if (NDIM == 2 .and. coord == af_xyz) then
          solution = boundary_value + solution_coeff * log(distance)
       else
          solution = boundary_value + solution_coeff * (1 - 1/distance)
       end if
    case (5)
       ! Triangle
       lsf = get_lsf(r)
       if (lsf <= 0.0_dp) then
          solution = boundary_value
       else
          solution = boundary_value * exp(-lsf)
       end if
    case default
       solution = 0.0_dp
    end select
  end function solution

  ! This routine sets the level set function
  subroutine set_lsf(box, iv)
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: iv
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM), norm_dr

    nc = box%n_cell
    norm_dr = norm2(box%dr)

    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, iv) = get_lsf(rr)
    end do; CLOSE_DO
  end subroutine set_lsf

  real(dp) function get_lsf(rr)
    real(dp), intent(in) :: rr(NDIM)
    real(dp)             :: qq(NDIM)
#if NDIM > 1
    real(dp)             :: dist1, dist2
#endif

    qq = (rr - solution_r0)/solution_radius

#if NDIM == 3
    ! Generalize shapes to 3D
    qq(1) = norm2(qq([1, 2]))
    qq(2) = qq(3)
    qq(3) = 0.0_dp
#endif

    select case (shape)
    case (1)
       get_lsf = norm2(qq) - 1.0_dp
#if NDIM > 1
    case (2)
       ! Heart centered on r0
       get_lsf = (qq(1)**2 + qq(2)**2 - 1)**3 - &
            qq(1)**2 * qq(2)**3
    case (3)
       ! Astroid
       get_lsf = abs(qq(1))**(2.0_dp/3) / 0.8_dp + &
            abs(qq(2))**(2.0_dp/3) / 1.5_dp - 0.8_dp
    case (4)
       ! Rhombus
       ! sharpness_t -> for sharpness of the top angle,
       ! larger sharpness_t equals more acute angle
       get_lsf = sharpness_t * abs(qq(1)) + abs(qq(2)) - 1.5_dp
    case (5)
       ! Triangle, uses signed distance from the triangle
       dist1 = GM_dist_line(rr, [solution_r0(1) - solution_r0(2)/sharpness_t, 0.0_dp], &
            solution_r0, 2)
       dist2 = GM_dist_line(rr, [solution_r0(1) + solution_r0(2)/sharpness_t, 0.0_dp], &
            solution_r0, 2)

       ! Determine sign of lsf function
       qq = rr - solution_r0
       get_lsf = sharpness_t * abs(qq(1)) + qq(2)

       ! Use sign in front of minimum distance
       get_lsf = sign(min(dist1, dist2), get_lsf)
    case (6)
       ! Heart v2
       get_lsf = qq(1)**2 + (qq(2) - abs(qq(1))**(2/3.0_dp))**2 - 1
    case (7)
       ! Spheroid
       get_lsf = sqrt(sharpness_t * qq(1)**2 + qq(2)**2) - 1.0_dp
#endif
    case default
       error stop "Invalid shape"
    end select

  end function get_lsf

  subroutine set_error(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM)

    nc = box%n_cell
    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, i_error) = box%cc(IJK, i_phi) - solution(rr)
    end do; CLOSE_DO
  end subroutine set_error

  subroutine set_lsf_mask(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc, ix

    nc = box%n_cell
    box%cc(DTIMES(:), i_lsf_mask) = 0

    ix = af_stencil_index(box, mg%operator_key)

    if (ix == af_stencil_none) error stop "No stencil stored"

    if (allocated(box%stencils(ix)%f)) then
       where (abs(box%stencils(ix)%f) > 0)
          box%cc(DTIMES(1:nc), i_lsf_mask) = 1
       end where
    end if

    ! This code can get the lsf mask that is used to check for roots
    ! ix_mask = af_stencil_index(box, mg_lsf_mask_key)
    ! if (ix_mask /= af_stencil_none) then
    !    do n = 1, size(box%stencils(ix_mask)%sparse_ix, 2)
    !       ix = box%stencils(ix_mask)%sparse_ix(:, n)
    !       box%cc(DINDEX(ix), i_lsf_mask) = 1
    !    end do
    ! end if
  end subroutine set_lsf_mask

  ! This routine sets boundary conditions for a box
  subroutine bc_solution(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type
    integer                 :: n

    bc_type = af_bc_dirichlet

    do n = 1, box%n_cell**(NDIM-1)
       bc_val(n) = solution(coords(:, n))
    end do
  end subroutine bc_solution

#if NDIM > 1
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
#endif

end program poisson_lsf_test
