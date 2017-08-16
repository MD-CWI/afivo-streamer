#include "../afivo/src/cpp_macros_$Dd.h"
!> Module for photoionization with the Helmholtz approximation
!>
!> The equations solved are nabla^2 phi_i - lambda_i^2 * phi_i = f, where f is
!> the photon production term, and there can be modes i = 1, ..., N. The photon
!> absorption term is then given by sum(c_i * phi_i).
!>
!> For the case N=2, the following coefficients can be used:
!> lambdas = 4425.36_dp, 750.06_dp 1/(m bar)
!> coeffs = 337555.5_dp, 19972.0_dp 1/(m bar)^2
!>
!> TODO: look up values for the case N=3
module m_photoi_helmh_$Dd
  use m_a$D_all
  use m_streamer

  implicit none
  private

  type(mg$D_t) :: mg_helm       ! Multigrid option struct

  ! Lambda squared, used internally by the routines
  real(dp) :: lambda2

  !> How many V-cycles to perform to update each mode
  integer :: num_vcycles = 2

  !> Number of Helmholtz modes to include
  integer               :: n_modes = 2
  !> Lambdas to use for lpl(phi) - lambda*phi = f; unit 1/(m bar)
  real(dp), allocatable :: lambdas(:)
  !> Weights corresponding to the lambdas; unit 1/(m bar)^2
  real(dp), allocatable :: coeffs(:)

  integer, allocatable :: i_modes(:)

  logical :: have_guess = .false.

  public :: photoi_helmh_initialize
  public :: photoi_helmh_compute
  public :: photoi_helmh_bc

contains

  !> Initialize options for Helmholtz photoionization
  subroutine photoi_helmh_initialize(cfg, is_used)
    use m_config
    use m_units_constants
    type(CFG_t), intent(inout) :: cfg
    logical, intent(in)        :: is_used
    integer                    :: n
    character(len=12)          :: name

    call CFG_add(cfg, "photoi_helmh%lambdas", [4425.36_dp, 750.06_dp], &
         "Lambdas to use for lpl(phi) - lambda*phi = f; unit 1/(m bar)", &
         .true.)
    call CFG_add(cfg, "photoi_helmh%coeffs", [337555.5_dp, 19972.0_dp], &
         "Weights corresponding to the lambdas; unit 1/(m bar)^2", &
         .true.)

    call CFG_add(cfg, "photoi_helmh%num_vcycles", num_vcycles, &
         "How many V-cycles to perform to update each mode")

    if (is_used) then
       call CFG_get_size(cfg, "photoi_helmh%lambdas", n_modes)
       allocate(lambdas(n_modes))
       allocate(coeffs(n_modes))
       call CFG_get(cfg, "photoi_helmh%lambdas", lambdas)
       call CFG_get(cfg, "photoi_helmh%coeffs", coeffs)
       call CFG_get(cfg, "photoi_helmh%num_vcycles", num_vcycles)

       ! Convert to correct units by multiplying with pressure in bar
       lambdas = lambdas * ST_gas_pressure  ! 1/m
       coeffs  = coeffs * ST_gas_pressure**2 ! 1/m^2

       ! Add required variables
       allocate(i_modes(n_modes))
       do n = 1, n_modes
          write(name, "(A,I0)") "helmh_", n
          i_modes(n) = ST_add_cc_variable(trim(name), .false.)
       end do

       ! Now set the multigrid options
       mg_helm%i_phi = i_modes(1) ! Will updated later on
       mg_helm%i_rhs = i_electron_old
       mg_helm%i_tmp = i_pos_ion_old

       ! Todo: check what good b.c. are
       mg_helm%sides_bc => a$D_bc_dirichlet_zero

#if $D == 2
       if (ST_cylindrical) then
          mg_helm%box_op => helmholtz_cyl_operator
          mg_helm%box_gsrb => helmholtz_cyl_gsrb
       else
          mg_helm%box_op => helmholtz_operator
          mg_helm%box_gsrb => helmholtz_gsrb
       end if
#elif $D == 3
       mg_helm%box_op => helmholtz_operator
       mg_helm%box_gsrb => helmholtz_gsrb
#endif

       call mg$D_init_mg(mg_helm)
    end if
  end subroutine photoi_helmh_initialize

  subroutine photoi_helmh_compute(tree)
    type(a$D_t), intent(inout) :: tree

    integer :: n, lvl, i, id

    call a$D_tree_clear_cc(tree, i_photo)

    do n = 1, n_modes
       lambda2 = lambdas(n)**2

       mg_helm%i_phi = i_modes(n)

       if (.not. have_guess) then
          call mg$D_fas_fmg(tree, mg_helm, .false., have_guess)
       else
          have_guess = .true.
          do i = 1, num_vcycles
             call mg$D_fas_vcycle(tree, mg_helm, .false.)
          end do
       end if

       !$omp parallel private(lvl, i, id)
       do lvl = 1, tree%highest_lvl
          !$omp do
          do i = 1, size(tree%lvls(lvl)%leaves)
             id = tree%lvls(lvl)%leaves(i)
             tree%boxes(id)%cc(DTIMES(:), i_photo) = &
                  tree%boxes(id)%cc(DTIMES(:), i_photo) - &
                  coeffs(n) * tree%boxes(id)%cc(DTIMES(:), i_modes(n))
          end do
          !$omp end do
       end do
       !$omp end parallel
    end do
  end subroutine photoi_helmh_compute

  subroutine photoi_helmh_bc(box, nb, iv, bc_type)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: nb ! Direction for the boundary condition
    integer, intent(in)         :: iv ! Index of variable
    integer, intent(out)        :: bc_type ! Type of boundary condition
    integer                     :: nc

    nc = box%n_cell

    select case (nb)
    case (a$D_neighb_lowx, a$D_neighb_highx)
       call a$D_bc_neumann_zero(box, nb, iv, bc_type)
#if $D == 2
    case (a$D_neighb_lowy, a$D_neighb_highy)
       call a$D_bc_dirichlet_zero(box, nb, iv, bc_type)
#elif $D == 3
    case (a$D_neighb_lowy, a$D_neighb_highy)
       call a$D_bc_neumann_zero(box, nb, iv, bc_type)
    case (a$D_neighb_lowz, a$D_neighb_highz)
       call a$D_bc_dirichlet_zero(box, nb, iv, bc_type)
#endif
    end select
  end subroutine photoi_helmh_bc

#if $D == 2
  !> Perform cylindrical Laplacian operator on a box
  subroutine helmholtz_cyl_operator(box, i_out, mg)
    type(box2_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    type(mg2_t), intent(in)     :: mg !< Multigrid options
    integer                     :: i, j, nc, i_phi, ioff
    real(dp)                    :: inv_dr_sq, rfac(2)

    nc        = box%n_cell
    inv_dr_sq = 1 / box%dr**2
    i_phi     = mg%i_phi
    ioff      = (box%ix(1)-1) * nc

    do j = 1, nc
       do i = 1, nc
          rfac = [i+ioff-1, i+ioff] / (i+ioff-0.5_dp)
          box%cc(i, j, i_out) = ( &
               rfac(1) * box%cc(i-1, j, i_phi) + &
               rfac(2) * box%cc(i+1, j, i_phi) + &
               box%cc(i, j-1, i_phi) + box%cc(i, j+1, i_phi) - &
               4 * box%cc(i, j, i_phi)) * inv_dr_sq - &
               lambda2 * box%cc(i, j, i_phi)
       end do
    end do
  end subroutine helmholtz_cyl_operator

  !> Perform Gauss-Seidel relaxation on box for a cylindrical Helmholtz operator
  subroutine helmholtz_cyl_gsrb(box, redblack_cntr, mg)
    type(box2_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    type(mg2_t), intent(in)     :: mg !< Multigrid options
    integer                     :: i, i0, j, nc, i_phi, i_rhs, ioff
    real(dp)                    :: dx2, rfac(2)

    dx2   = box%dr**2
    nc    = box%n_cell
    i_phi = mg%i_phi
    i_rhs = mg%i_rhs
    ioff  = (box%ix(1)-1) * nc

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          rfac = [i+ioff-1, i+ioff] / (i+ioff-0.5_dp)
          box%cc(i, j, i_phi) = 1/(4 + lambda2 * dx2) * ( &
               rfac(1) * box%cc(i-1, j, i_phi) + &
               rfac(2) * box%cc(i+1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi) - &
               dx2 * box%cc(i, j, i_rhs))
       end do
    end do
  end subroutine helmholtz_cyl_gsrb
#endif

  !> Perform Helmholtz operator on a box
  subroutine helmholtz_operator(box, i_out, mg)
    type(box$D_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    type(mg$D_t), intent(in)     :: mg !< Multigrid options
    integer                     :: i, j, nc, i_phi
    real(dp)                    :: inv_dr_sq
#if $D == 3
    integer                     :: k
#endif

    nc = box%n_cell
    inv_dr_sq = 1 / box%dr**2
    i_phi = mg%i_phi

#if $D == 2
    do j = 1, nc
       do i = 1, nc
          box%cc(i, j, i_out) = inv_dr_sq * (box%cc(i-1, j, i_phi) + &
               box%cc(i+1, j, i_phi) + box%cc(i, j-1, i_phi) + &
               box%cc(i, j+1, i_phi) - 4 * box%cc(i, j, i_phi)) - &
               lambda2 * box%cc(i, j, i_phi)
       end do
    end do
#elif $D == 3
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             box%cc(i, j, k, i_out) = inv_dr_sq * (box%cc(i-1, j, k, i_phi) + &
                  box%cc(i+1, j, k, i_phi) + box%cc(i, j-1, k, i_phi) + &
                  box%cc(i, j+1, k, i_phi) + box%cc(i, j, k-1, i_phi) + &
                  box%cc(i, j, k+1, i_phi) - 6 * box%cc(i, j, k, i_phi)) - &
                  lambda2 * box%cc(i, j, k, i_phi)
          end do
       end do
    end do
#endif
  end subroutine helmholtz_operator

  !> Perform Gauss-Seidel relaxation on box for a Helmholtz operator
  subroutine helmholtz_gsrb(box, redblack_cntr, mg)
    type(box$D_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    type(mg$D_t), intent(in)     :: mg !< Multigrid options
    integer                     :: i, i0, j, nc, i_phi, i_rhs
    real(dp)                    :: dx2
#if $D == 3
    integer                     :: k
#endif

    dx2   = box%dr**2
    nc    = box%n_cell
    i_phi = mg%i_phi
    i_rhs = mg%i_rhs

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
#if $D == 2
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          box%cc(i, j, i_phi) = 1/(4 + lambda2 * dx2) * ( &
               box%cc(i+1, j, i_phi) + box%cc(i-1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi) - &
               dx2 * box%cc(i, j, i_rhs))
       end do
    end do
#elif $D == 3
    do k = 1, nc
       do j = 1, nc
          i0 = 2 - iand(ieor(redblack_cntr, k+j), 1)
          do i = i0, nc, 2
             box%cc(i, j, k, i_phi) = 1/(6 + lambda2 * dx2) * ( &
                  box%cc(i+1, j, k, i_phi) + box%cc(i-1, j, k, i_phi) + &
                  box%cc(i, j+1, k, i_phi) + box%cc(i, j-1, k, i_phi) + &
                  box%cc(i, j, k+1, i_phi) + box%cc(i, j, k-1, i_phi) - &
                  dx2 * box%cc(i, j, k, i_rhs))
          end do
       end do
    end do
#endif
  end subroutine helmholtz_gsrb

end module m_photoi_helmh_$Dd


