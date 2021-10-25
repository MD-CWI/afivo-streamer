#include "../afivo/src/cpp_macros.h"
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
module m_photoi_helmh
  use m_af_all
  use m_streamer

  implicit none
  private

  type(mg_t), allocatable :: mg_helm(:) ! Multigrid option struct

  !> Maximum number of FMG cycles to perform to update each mode
  integer, parameter :: max_fmg_cycles = 10

  !> Maximum residual relative to max(|rhs|)
  real(dp) :: max_rel_residual = 1.0d-2

  !> Number of Helmholtz modes to include
  integer               :: n_modes = 2
  !> Lambdas to use for lpl(phi) - lambda*phi = f; unit 1/(m bar)
  real(dp), allocatable :: lambdas(:)
  !> Weights corresponding to the lambdas; unit 1/(m bar)^2
  real(dp), allocatable :: coeffs(:)

  integer, allocatable :: i_modes(:)

  public :: photoi_helmh_initialize
  public :: photoi_helmh_compute
  public :: photoi_helmh_bc

contains

  !> Initialize options for Helmholtz photoionization
  subroutine photoi_helmh_initialize(tree, cfg, is_used, eta)
    use m_config
    use m_units_constants
    use m_gas
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg
    logical, intent(in)        :: is_used !< Whether the module is used
    real(dp), intent(in)       :: eta !< Efficiency
    integer                    :: n, ix
    character(len=12)          :: name
    character(len=12)          :: author = "Bourdon-3"
    real(dp)                   :: frac_O2, dummy(0)

    call CFG_add_get(cfg, "photoi_helmh%author", author, &
         "Can be Luque (default), Bourdon-2, Bourdon-3 or custom")
    call CFG_add(cfg, "photoi_helmh%lambdas", dummy, &
         "Lambdas to use in Helmholtz eq; unit 1/(m bar)", &
         .true.)
    call CFG_add(cfg, "photoi_helmh%coeffs", dummy, &
         "Weights corresponding to the lambdas; unit 1/(m bar)^2", &
         .true.)
    call CFG_add_get(cfg, "photoi_helmh%max_rel_residual", max_rel_residual, &
         "Maximum residual relative to max(|rhs|)")

    ! Exit here if the module is not used
    if (.not. is_used) return

    ! Whether oxygen is required is tested below, the idea being that we will in
    ! the future add photoionization cases for other gases
    ix = gas_index("O2")
    if (ix /= -1) then
       frac_O2 = gas_fractions(ix)
    else
       frac_O2 = 0.0_dp         ! No oxygen
    end if

    select case (author)
    case ("Luque")
       if (frac_O2 <= 0.0_dp) error stop "Photoionization: no oxygen present"
       n_modes = 2
       lambdas = [4425.38_dp, 750.06_dp]
       coeffs  = [337557.38_dp, 19972.14_dp]

       ! Convert to correct units by multiplying with pressure in bar
       lambdas = lambdas * (frac_O2/0.2_dp) * gas_pressure  ! 1/m
       coeffs  = coeffs * ((frac_O2/0.2_dp) * gas_pressure)**2 ! 1/m^2

       ! Note that for Luque photoi_eta must be 1.0 (since we must not multiply
       ! coeffs of Luque by photoi_eta)
       if (abs(eta - 1.0_dp) > 0) &
            error stop "With Luque photoionization, photoi%eta should be 1.0"
    case ("Bourdon-2")
       if (frac_O2 <= 0.0_dp) error stop "Photoionization: no oxygen present"
       !> Bourdon two terms lambdas in SI units are: [7305.62,44081.25]
       !> Bourdon two terms coeffs in SI units are: [11814508.38,998607256]
       n_modes = 2
       lambdas = [7305.62_dp, 44081.25_dp]
       coeffs  = [11814508.38_dp, 998607256.0_dp]

       ! Convert to correct units by multiplying with pressure in bar
       lambdas = lambdas * frac_O2 * gas_pressure  ! 1/m
       ! Important note: in the case of Bourdon coeffs must be multiplied by
       ! photoi_eta, however we do not multiply it here but it is considered in
       ! set_photoionization_rate calculation as [photoi_eta * quench_fac]
       ! please check m_photoi.f90
       coeffs  = coeffs * (frac_O2 * gas_pressure)**2 ! 1/m^2
    case ("Bourdon-3")
       if (frac_O2 <= 0.0_dp) error stop "Photoionization: no oxygen present"
       !> Bourdon three terms lambdas in SI units are: [4147.85   10950.93  66755.67]
       !> bourdon three terms coeffs in SI units are: [ 1117314.935  28692377.5  2748842283 ]
       n_modes = 3
       lambdas = [4147.85_dp, 10950.93_dp, 66755.67_dp]
       coeffs  = [1117314.935_dp, 28692377.5_dp, 2748842283.0_dp]

       ! Convert to correct units by multiplying with pressure in bar
       lambdas = lambdas * frac_O2 * gas_pressure  ! 1/m
       ! Important note: in the case of Bourdon coeffs must be multiplied by
       ! photoi_eta, however we do not multiply it here but it is considered in
       ! set_photoionization_rate calculation as [photoi_eta * quench_fac]
       ! please check m_photoi.f90
       coeffs  = coeffs * (frac_O2 * gas_pressure)**2 ! 1/m^2
    case ("custom")
       call CFG_get_size(cfg, "photoi_helmh%lambdas", n_modes)
       allocate(lambdas(n_modes))
       allocate(coeffs(n_modes))
       call CFG_get(cfg, "photoi_helmh%lambdas", lambdas)
       call CFG_get(cfg, "photoi_helmh%coeffs", coeffs)
       ! Convert to correct units by multiplying with pressure in bar, but
       ! perform no other normalization.
       lambdas = lambdas * gas_pressure  ! 1/m
       coeffs  = coeffs * gas_pressure**2 ! 1/m^2
    case default
       print *, "Unknown photoi_helmh_author: ", trim(author)
       error stop
    end select

    ! Add required variables
    allocate(i_modes(n_modes))
    do n = 1, n_modes
       write(name, "(A,I0)") "helmh_", n
       call af_add_cc_variable(tree, trim(name), write_out=.false., ix=i_modes(n))
    end do

    ! Now set the multigrid options
    allocate(mg_helm(n_modes))

    do n = 1, n_modes
       mg_helm(n)%i_phi = i_modes(n)
       mg_helm(n)%i_rhs = i_rhs
       mg_helm(n)%i_tmp = i_tmp
       mg_helm(n)%sides_bc => photoi_helmh_bc
       mg_helm(n)%helmholtz_lambda = lambdas(n)**2
       mg_helm(n)%operator_type = mg_normal_operator
       mg_helm(n)%prolongation_type = mg_prolong_linear
    end do
  end subroutine photoi_helmh_initialize

  subroutine photoi_helmh_compute(tree, i_photo)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: i_photo

    integer :: n, lvl, i, id
    real(dp) :: residu, max_rhs

    call af_tree_clear_cc(tree, i_photo)

    call af_tree_maxabs_cc(tree, mg_helm(1)%i_rhs, max_rhs)
    max_rhs = max(max_rhs, sqrt(epsilon(1.0_dp)))

    do n = 1, n_modes
       if (.not. mg_helm(n)%initialized) then
          if (n > 1) then
             ! Re-use prolongation stencil
             mg_helm(n)%prolongation_key = mg_helm(1)%prolongation_key
          end if

          call mg_init(tree, mg_helm(n))
       end if

       do i = 1, max_fmg_cycles
          call mg_fas_fmg(tree, mg_helm(n), .true., .true.)
          call af_tree_maxabs_cc(tree, mg_helm(n)%i_tmp, residu)
          ! print *, n, i, residu/max_rhs
          if (residu/max_rhs < max_rel_residual) exit
       end do

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

  !> @todo Think about good (and efficient) boundary conditions for Helmholtz
  !> equations, see also Bourdon et al. PSST 2017. Setting a Dirichlet zero b.c.
  !> is inaccurate if the streamer gets close to that boundary, otherwise it
  !> should be quite reasonable.
  subroutine photoi_helmh_bc(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type
    integer                 :: nc

    nc = box%n_cell

    if (af_neighb_dim(nb) == NDIM) then
       bc_type = af_bc_dirichlet
       bc_val  = 0.0_dp
    else
       bc_type = af_bc_neumann
       bc_val  = 0.0_dp
    end if
  end subroutine photoi_helmh_bc

end module m_photoi_helmh
