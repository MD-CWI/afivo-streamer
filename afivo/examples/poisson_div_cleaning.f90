#include "../src/cpp_macros.h"
!> \example poisson_div_cleaning.f90
!>
!> Example showing how to use multigrid for ensuring a vector field is
!> divergence free
program poisson_div_cleaning
  use m_af_all
  use m_gaussians

  implicit none

  integer, parameter :: box_size     = 8
  integer, parameter :: n_iterations = 10
  integer            :: domain_size(NDIM)
  real(dp)           :: domain_length(NDIM)
  real(dp)           :: domain_origin(NDIM)
  integer            :: i_phi
  integer            :: i_div_old
  integer            :: i_div_new
  integer            :: i_tmp
  integer            :: if_B_old
  integer            :: if_B_new
  integer            :: i_B(3)
  type(af_t)         :: tree
  type(ref_info_t)   :: refine_info
  integer            :: mg_iter
  character(len=100) :: fname
  type(mg_t)         :: mg

  real(dp), parameter :: alphar_ = 0.64d0
  real(dp), parameter :: ar_     = 0.25d0
  real(dp), parameter :: B0      = 1.0_dp

  print *, "Running poisson_div_cleaning_" // DIMNAME
  print *, "Number of threads", af_get_max_threads()

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "div_old", ix=i_div_old)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)
  call af_add_cc_variable(tree, "div_new", ix=i_div_new)
  call af_add_cc_variable(tree, "Bx", ix=i_B(1))
  call af_add_cc_variable(tree, "By", ix=i_B(2))
  call af_add_cc_variable(tree, "Bz", ix=i_B(3))
  call af_add_fc_variable(tree, "B_old", ix=if_B_old)
  call af_add_fc_variable(tree, "B_new", ix=if_B_new)

  domain_length(:)  = 1.0_dp
  domain_origin(:) = -0.5_dp * domain_length(:)
  domain_size(:) = box_size

  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       domain_origin+domain_length, &
       domain_size, &
       r_min=domain_origin)

  do
     call af_loop_box(tree, set_initial_condition)
     call af_adjust_refinement(tree, refine_routine, refine_info)
     if (refine_info%n_add == 0) exit
  end do

  call af_print_info(tree)

  mg%i_phi    = i_phi ! Solution variable
  mg%i_rhs    = i_div_old ! Right-hand side variable
  mg%i_tmp    = i_tmp ! Variable for temporary space
  mg%sides_bc => af_bc_dirichlet_zero

  call mg_init(tree, mg)

  do mg_iter = 1, n_iterations
     call mg_fas_fmg(tree, mg, set_residual=.true., have_guess=(mg_iter>1))
     call af_loop_box(tree, set_error)
     write(fname, "(A,I0)") "poisson_div_cleaning_" // DIMNAME // "_", mg_iter
     call af_write_silo(tree, trim(fname), dir="output")
  end do

contains

  ! Return the refinement flags for box
  subroutine refine_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))

    if (box%lvl < 3) then
       cell_flags = af_do_ref
    else
       cell_flags = af_keep_ref
    end if
  end subroutine refine_routine

  real(dp) function Bx(xyz)
    real(dp), intent(in) :: xyz(NDIM)
    real(dp)             :: r, t, phi, Bphi
    r    = norm2(xyz(1:2))
    t    = r / ar_
    Bphi = B0 / (t * (1 + t**2)**alphar_) * &
         sqrt(((1 + t**2)**(2 * alphar_) - 1 - 2 * alphar_ * t**2) / &
         (2 * alphar_ - 1))
    phi  = atan2(xyz(2), xyz(1))
    Bx   = -sin(phi) * Bphi
  end function Bx

  real(dp) function By(xyz)
    real(dp), intent(in) :: xyz(NDIM)
    real(dp)             :: r, t, phi, Bphi
    r    = norm2(xyz(1:2))
    t    = r / ar_
    Bphi = B0 / (t * (1 + t**2)**alphar_) * &
         sqrt(((1 + t**2)**(2 * alphar_) - 1 - 2 * alphar_ * t**2) / &
         (2 * alphar_ - 1))
    phi  = atan2(xyz(2), xyz(1))
    By   = cos(phi) * Bphi
  end function By

  real(dp) function Bz(xyz)
    real(dp), intent(in) :: xyz(NDIM)
    real(dp)             :: r, t
    r  = norm2(xyz(1:2))
    t  = r / ar_
    Bz = B0 / (1 + t**2)**alphar_
  end function Bz

  ! This routine sets the initial conditions for each box
  subroutine set_initial_condition(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK, nc
    real(dp)                   :: xyz(NDIM)

    nc = box%n_cell

    ! x-component
    do k = 1, nc; do j = 1, nc; do i = 1, nc+1
       box%fc(IJK, 1, if_B_old) = Bx(af_r_fc(box, 1, [IJK]))
    end do; end do; end do

    ! y-component
    do k = 1, nc; do j = 1, nc+1; do i = 1, nc
       box%fc(IJK, 2, if_B_old) = By(af_r_fc(box, 2, [IJK]))
    end do; end do; end do

    ! z-component
    do k = 1, nc+1; do j = 1, nc; do i = 1, nc
       box%fc(IJK, 3, if_B_old) = Bz(af_r_fc(box, 3, [IJK]))
    end do; end do; end do

    do KJI_DO(1, nc)
       xyz = af_r_cc(box, [IJK])

       ! Cell-centered field for visualization
       box%cc(IJK, i_B(:)) = [Bx(xyz), By(xyz), Bz(xyz)]

       ! Initial divergence of field
       box%cc(IJK, i_div_old) = &
            (box%fc(i+1, j, k, 1, if_B_old) - &
            box%fc(i, j, k, 1, if_B_old)) / box%dr(1) + &
            (box%fc(i, j+1, k, 2, if_B_old) - &
            box%fc(i, j, k, 2, if_B_old)) / box%dr(2) + &
            (box%fc(i, j, k+1, 3, if_B_old) - &
            box%fc(i, j, k, 3, if_B_old)) / box%dr(3)
    end do; CLOSE_DO
  end subroutine set_initial_condition

  subroutine set_error(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK, nc
    real(dp)                   :: inv_dr(NDIM)

    nc     = box%n_cell
    inv_dr = 1 / box%dr

    do k = 1, nc; do j = 1, nc; do i = 1, nc+1
       ! x-component
       box%fc(IJK, 1, if_B_new) = box%fc(IJK, 1, if_B_old) - inv_dr(1) * &
            (box%cc(i, j, k, i_phi) - box%cc(i-1, j, k, i_phi))
    end do; end do; end do

    do k = 1, nc; do j = 1, nc+1; do i = 1, nc
       ! y-component
       box%fc(IJK, 2, if_B_new) = box%fc(IJK, 2, if_B_old) - inv_dr(2) * &
            (box%cc(i, j, k, i_phi) - box%cc(i, j-1, k, i_phi))
    end do; end do; end do

    do k = 1, nc+1; do j = 1, nc; do i = 1, nc
       ! z-component
       box%fc(IJK, 3, if_B_new) = box%fc(IJK, 3, if_B_old) - inv_dr(3) * &
            (box%cc(i, j, k, i_phi) - box%cc(i, j, k-1, i_phi))
    end do; end do; end do

    ! Divergence after cleaning
    do KJI_DO(1, nc)
       box%cc(IJK, i_div_new) = &
            (box%fc(i+1, j, k, 1, if_B_new) - &
            box%fc(i, j, k, 1, if_B_new)) / box%dr(1) + &
            (box%fc(i, j+1, k, 2, if_B_new) - &
            box%fc(i, j, k, 2, if_B_new)) / box%dr(2) + &
            (box%fc(i, j, k+1, 3, if_B_new) - &
            box%fc(i, j, k, 3, if_B_new)) / box%dr(3)
    end do; CLOSE_DO
  end subroutine set_error

end program poisson_div_cleaning
