#include "../src/cpp_macros.h"
!> \example poisson_lsf_test.f90
!>
!> Test of the level-set functionality of the Poisson solver
program poisson_lsf_test
  use m_af_all
  use m_dielectric

  implicit none

  integer, parameter :: box_size         = 8
  integer, parameter :: n_iterations     = 10
  integer            :: max_refine_level = 3
  integer            :: i_phi
  integer            :: i_rhs
  integer            :: i_tmp
  integer            :: i_lsf
  integer            :: i_error
  integer            :: i_field
  integer            :: i_field_norm

  real(dp), parameter :: boundary_value    = 1.0_dp
  real(dp), parameter :: solution_coeff    = 1.0_dp
  real(dp), parameter :: solution_radius   = 0.25_dp
  real(dp)            :: solution_r0(NDIM) = 0.5_dp

  type(af_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: n, mg_iter, coord, n_args
  real(dp)           :: residu, max_error, max_field
  character(len=100) :: fname, argv
  type(mg_t)         :: mg

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)
  call af_add_cc_variable(tree, "lsf", ix=i_lsf)
  call af_add_cc_variable(tree, "error", ix=i_error)
  call af_add_cc_variable(tree, "field_norm", ix=i_field_norm)
  call af_add_fc_variable(tree, "field", ix=i_field)

  call af_set_cc_methods(tree, i_lsf, funcval=set_lsf)
  call af_set_cc_methods(tree, i_field_norm, af_bc_neumann_zero)

  ! If an argument is given, switch to cylindrical coordinates in 2D
  n_args = command_argument_count()
  if (n_args > 2) &
       stop "Usage: ./poisson_lsf_test [cyl] [max_refinement_level]"

  coord = af_xyz
  do n = 1, n_args
     call get_command_argument(n, argv)
     if (argv == 'cyl') then
        coord = af_cyl
        ! Place solution on axis
        solution_r0(1) = 0.0_dp
     else
        read(argv, *) max_refine_level
     end if
  end do

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [DTIMES(1.0_dp)], &
       4 * [DTIMES(box_size)], &
       coord=coord)

  do n = 1, 100
     call af_adjust_refinement(tree, ref_routine, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  call af_print_info(tree)

  mg%i_phi    = i_phi
  mg%i_rhs    = i_rhs
  mg%i_tmp    = i_tmp
  mg%i_lsf    = i_lsf
  mg%sides_bc => bc_solution
  mg%lsf_boundary_value = boundary_value

  call mg_init(tree, mg)

  do mg_iter = 1, n_iterations
     call mg_fas_fmg(tree, mg, .true., mg_iter>1)
     call mg_compute_phi_gradient(tree, mg, i_field, 1.0_dp, i_field_norm)
     call af_gc_tree(tree, [i_field_norm])
     call af_loop_box(tree, set_error)

     ! Determine the minimum and maximum residual and error
     call af_tree_maxabs_cc(tree, i_tmp, residu)
     call af_tree_maxabs_cc(tree, i_error, max_error)
     call af_tree_maxabs_cc(tree, i_field_norm, max_field)
     write(*, "(I8,3E14.5)") mg_iter, residu, max_error, max_field

     write(fname, "(A,I0)") "poisson_lsf_test_" // DIMNAME // "_", mg_iter
     call af_write_silo(tree, trim(fname), dir="output")
  end do

contains

  ! Set refinement flags for box
  subroutine ref_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))
    integer                 :: nc
    integer, parameter      :: refinement_type = 1

    nc = box%n_cell

    select case (refinement_type)
    case (1)
       ! Uniform refinement
       if (box%lvl < max_refine_level) then
          cell_flags(DTIMES(:)) = af_do_ref
       else
          cell_flags(DTIMES(:)) = af_keep_ref
       end if

    case (2)
       ! Only refine at boundary
       if (box%lvl < max_refine_level .and. &
            minval(box%cc(DTIMES(:), i_lsf)) * &
            maxval(box%cc(DTIMES(:), i_lsf)) <= 0) then
          cell_flags(DTIMES(:)) = af_do_ref
       else
          cell_flags(DTIMES(:)) = af_keep_ref
       end if

    case (3)
       ! 'Bad' refinement to test the method
       if (norm2(box%r_min - solution_r0) < solution_radius .and. &
            box%lvl < max_refine_level) then
          cell_flags(DTIMES(:)) = af_do_ref
       else
          cell_flags(DTIMES(:)) = af_keep_ref
       end if
    end select
  end subroutine ref_routine

  real(dp) function solution(r)
    real(dp), intent(in) :: r(NDIM)
    real(dp) :: distance

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
  end function solution

  ! This routine sets the level set function
  subroutine set_lsf(box, iv)
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: iv
    integer                    :: IJK, nc
    real(dp)                   :: distance, rr(NDIM)

    nc = box%n_cell

    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       distance = norm2(rr-solution_r0) / solution_radius
       box%cc(IJK, iv) = distance - 1.0_dp
    end do; CLOSE_DO
  end subroutine set_lsf

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

end program
