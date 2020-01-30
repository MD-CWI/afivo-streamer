#include "../src/cpp_macros.h"
!> \example dielectric_surface.f90
!>
!> Example showing how to include a dielectric surface
program dielectric_surface
  use m_af_all
  use m_dielectric

  implicit none

  integer, parameter :: box_size     = 8
  integer, parameter :: n_iterations = 10
  integer            :: i_phi
  integer            :: i_rhs
  integer            :: i_tmp
  integer            :: i_eps
  integer            :: i_fld_norm_fc
  integer            :: i_fld_norm_cc
  integer            :: i_fld_cc(NDIM)
  integer            :: i_fld_fc

  ! The dielectric constant used in this example
  double precision, parameter :: epsilon_high = 1000.0_dp

  ! Where the interface is located
  real(dp), parameter :: interface_location = 0.25_dp
  ! Along which dimension the interface occurs
  integer, parameter :: interface_dimension = 1

  type(af_t)         :: tree
  type(ref_info_t)   :: ref_info
  type(dielectric_t) :: dielectric
  integer            :: n, mg_iter
  real(dp)           :: residu
  character(len=100) :: fname
  type(mg_t)         :: mg
  integer, allocatable :: ref_links(:, :)

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)
  call af_add_cc_variable(tree, "eps", ix=i_eps)

  ! Include both cell-centered and face-centered electric field for testing
  call af_add_cc_variable(tree, "fld_cc_x", ix=i_fld_cc(1))
  call af_add_cc_variable(tree, "fld_cc_y", ix=i_fld_cc(2))
#if NDIM == 3
  call af_add_cc_variable(tree, "fld_cc_z", ix=i_fld_cc(3))
#endif
  call af_add_cc_variable(tree, "fld_norm_cc", ix=i_fld_norm_cc)
  call af_add_cc_variable(tree, "fld_norm_fc", ix=i_fld_norm_fc)
  call af_add_fc_variable(tree, "fld_fc", i_fld_fc)

  call af_set_cc_methods(tree, i_eps, af_bc_neumann_zero, &
       af_gc_prolong_copy, prolong=af_prolong_zeroth)

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [DTIMES(1.0_dp)], &
       [DTIMES(box_size)])

  call af_loop_box(tree, set_init_cond)

  do n = 1, 2
     call af_adjust_refinement(tree, ref_routine, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  dielectric%i_eps = i_eps
  call dielectric_initialize(tree, i_eps, dielectric, 1)

  call dielectric_set_values(tree, dielectric, 1, sigma_function)

  do n = 1, 3
     call dielectric_get_refinement_links(dielectric, ref_links)
     call af_adjust_refinement(tree, ref_random, ref_info, ref_links=ref_links)
     call dielectric_update_after_refinement(tree, dielectric, ref_info)
  end do

  call dielectric_surface_charge_to_rhs(tree, dielectric, 1, i_rhs, 1.0_dp)

  mg%i_phi        = i_phi
  mg%i_rhs        = i_rhs
  mg%i_tmp        = i_tmp
  mg%i_eps        = i_eps
  mg%sides_bc     => bc_phi

  call mg_init(tree, mg)

  do mg_iter = 1, n_iterations
     call mg_fas_fmg(tree, mg, .true., mg_iter>1)
     call af_loop_box(tree, compute_fields)
     call dielectric_correct_field_fc(tree, dielectric, 1, i_fld_fc, i_phi, 1.0_dp)
     call dielectric_correct_field_cc(tree, dielectric, 1, i_fld_cc, i_phi, 1.0_dp)
     call af_loop_box(tree, compute_field_norms)

     ! Determine the minimum and maximum residual and error
     call af_tree_maxabs_cc(tree, i_tmp, residu)
     write(*, "(I8,Es14.5)") mg_iter, residu

     write(fname, "(A,I0)") "dielectric_surface_" // DIMNAME // "_", mg_iter
     call af_write_silo(tree, trim(fname), dir="output")
  end do

contains

  ! Set refinement flags for box
  subroutine ref_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))

    if (box%lvl < 3) then
       cell_flags(DTIMES(:)) = af_do_ref
    else
       cell_flags(DTIMES(:)) = af_keep_ref
    end if
  end subroutine ref_routine

  subroutine ref_random(box, cell_flags)
    type(box_t), intent(in) :: box ! A list of all boxes in the tree
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))
    real(dp)                :: rr

    ! Draw a [0, 1) random number
    call random_number(rr)

    if (rr < 0.5_dp**NDIM .and. box%lvl < 6) then
       cell_flags = af_do_ref
    else if (box%lvl > 3) then
       cell_flags = af_rm_ref
    else
       cell_flags = af_keep_ref
    end if
  end subroutine ref_random

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM)

    nc = box%n_cell

    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])

       ! Change epsilon in part of the domain
       if (rr(interface_dimension) < interface_location) then
          box%cc(IJK, i_eps) = epsilon_high
       else
          box%cc(IJK, i_eps) = 1.0_dp
       end if
    end do; CLOSE_DO

  end subroutine set_init_cond

  ! This routine sets boundary conditions for a box
  subroutine bc_phi(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type
    integer                 :: nc

    nc = box%n_cell

    ! Below the solution is specified in the approriate ghost cells
    select case (nb)
    case (2 * interface_dimension - 1)
       bc_type = af_bc_dirichlet
       bc_val = 0.0_dp
    case (2 * interface_dimension)
       bc_type = af_bc_dirichlet
       bc_val = 1.0_dp
    case default
       bc_type = af_bc_neumann
       bc_val = 0.0_dp
    end select
  end subroutine bc_phi

  real(dp) function sigma_function(r)
    real(dp), intent(in) :: r(NDIM)

    ! Screen electric field outside dielectric
    ! sigma_function = -epsilon_high/interface_location

    ! Screen electric field inside dielectric
    ! sigma_function = 1/(1 - interface_location)

    ! Equal electric field on left and right
    ! sigma_function = 1 - epsilon_high

    ! No surface charge
    sigma_function = 0.0_dp
  end function sigma_function

  subroutine compute_fields(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc
    real(dp)                   :: inv_dr(NDIM)

    nc     = box%n_cell
    inv_dr = 1 / box%dr

#if NDIM == 2
    box%fc(1:nc+1, 1:nc, 1, i_fld_fc) = inv_dr(1) * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi))
    box%fc(1:nc, 1:nc+1, 2, i_fld_fc) = inv_dr(2) * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi))
    box%cc(1:nc, 1:nc, i_fld_cc(1)) = 0.5_dp * inv_dr(1) * &
         (box%cc(0:nc-1, 1:nc, i_phi) - box%cc(2:nc+1, 1:nc, i_phi))
    box%cc(1:nc, 1:nc, i_fld_cc(2)) = 0.5_dp * inv_dr(2) * &
         (box%cc(1:nc, 0:nc-1, i_phi) - box%cc(1:nc, 2:nc+1, i_phi))
#elif NDIM == 3
    error stop
#endif

  end subroutine compute_fields

  subroutine compute_field_norms(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc

    nc = box%n_cell

#if NDIM == 2
    box%cc(1:nc, 1:nc, i_fld_norm_fc) = 0.5_dp * sqrt(&
         (box%fc(1:nc, 1:nc, 1, i_fld_fc) + &
         box%fc(2:nc+1, 1:nc, 1, i_fld_fc))**2 + &
         (box%fc(1:nc, 1:nc, 2, i_fld_fc) + &
         box%fc(1:nc, 2:nc+1, 2, i_fld_fc))**2)
    box%cc(1:nc, 1:nc, i_fld_norm_cc) = sqrt(&
         box%cc(1:nc, 1:nc, i_fld_cc(1))**2 + &
         box%cc(1:nc, 1:nc, i_fld_cc(2))**2)
#elif NDIM == 3
    error stop
#endif
  end subroutine compute_field_norms

end program dielectric_surface
