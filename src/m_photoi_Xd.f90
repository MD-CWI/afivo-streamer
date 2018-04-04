#include "../afivo/src/cpp_macros_$Dd.h"
!> Top-module for photoionization, which can make use of different methods
module m_photoi_$Dd
  use m_photoi_mc
  use m_photoi_helmh_$Dd
  use m_a$D_all
  use m_streamer

  implicit none
  private

  ! Whether photoionization is enabled
  logical, protected, public :: photoi_enabled = .false.

  ! Which photoionization method to use (helmholtz, montecarlo)
  character(len=ST_slen) :: photoi_method = 'helmholtz'

  ! Photoionization efficiency factor, typically around 0.05-0.1, not for Helmholtz-Luque should be 1.0
  real(dp) :: photoi_eta = 0.05_dp

  ! Update photoionization every N time step
  integer, protected, public :: photoi_per_steps = 10

  public :: photoi_initialize
  public :: photoi_set_src
  ! Imported from Helmholtz module
  public :: photoi_helmh_bc

contains

  !> Initialize photoionization parameters
  subroutine photoi_initialize(cfg)
    use m_config
    type(CFG_t), intent(inout) :: cfg          !< The configuration for the simulation

    call CFG_add_get(cfg, "photoi%enabled", photoi_enabled, &
         "Whether photoionization is enabled")
    call CFG_add_get(cfg, "photoi%per_steps", photoi_per_steps, &
         "Update photoionization every N time step")
    call CFG_add_get(cfg, "photoi%method", photoi_method, &
         "Which photoionization method to use (helmholtz, montecarlo)")
    call CFG_add_get(cfg, "photoi%eta", photoi_eta, &
         "Photoionization efficiency factor, typically around 0.05-0.1")
    if (photoi_eta <= 0.0_dp) error stop "photoi%eta <= 0.0"
    if (photoi_eta > 1.0_dp) error stop "photoi%eta > 1.0"

    if (photoi_enabled) then
       i_photo = ST_add_cc_variable("photo", .not. ST_small_output)
    end if

    select case (photoi_method)
       case ("helmholtz")
          call photoi_helmh_initialize(cfg, .true.)
          call phmc_initialize(cfg, .false.)
       case ("montecarlo")
          call photoi_helmh_initialize(cfg, .false.)
          call phmc_initialize(cfg, .true.)
       case default
          print *, "Unknown photoi_method: ", trim(photoi_method)
          error stop
    end select
  end subroutine photoi_initialize

  !> Sets the photoionization
  subroutine photoi_set_src(tree, dt)
    use m_units_constants

    type(a$D_t), intent(inout)     :: tree
    real(dp), intent(in), optional :: dt
    real(dp), parameter            :: p_quench = 40.0e-3_dp
    real(dp)                       :: quench_fac

    ! Compute quench factor, because some excited species will be quenched by
    ! collisions, preventing the emission of a UV photon
    quench_fac = p_quench / (ST_gas_pressure + p_quench)

    ! Set photon production rate per cell, which is proportional to the
    ! ionization rate.
    call a$D_loop_box_arg(tree, set_photoionization_rate, &
         [photoi_eta * quench_fac], .true.)

    select case (photoi_method)
    case ("helmholtz")
       ! Use Helmholtz approximation
       call photoi_helmh_compute(tree)
    case ("montecarlo")
       if (phmc_physical_photons) then
#if $D == 2
          call phmc_set_src_$Dd(tree, ST_rng, i_electron_old, &
               i_photo, ST_cylindrical, dt)
#elif $D == 3
          call phmc_set_src_$Dd(tree, ST_rng, i_electron_old, i_photo, dt)
#endif
       else
#if $D == 2
          call phmc_set_src_$Dd(tree, ST_rng, i_electron_old, &
               i_photo, ST_cylindrical)
#elif $D == 3
          call phmc_set_src_$Dd(tree, ST_rng, i_electron_old, i_photo)
#endif
       end if
    end select

  end subroutine photoi_set_src

  !> Sets the photoionization_rate
  subroutine set_photoionization_rate(box, coeff)
    type(box$D_t), intent(inout) :: box
    real(dp), intent(in)        :: coeff(:)
    integer                     :: IJK, nc
    real(dp)                    :: fld, alpha, mobility, tmp
    type(LT_loc_t)              :: loc

    nc = box%n_cell

    do KJI_DO(1,nc)
       fld      = box%cc(IJK, i_electric_fld)
       loc      = LT_get_loc(ST_td_tbl, fld)
       alpha    = LT_get_col_at_loc(ST_td_tbl, i_alpha, loc)
       mobility = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)

       tmp = fld * mobility * alpha * box%cc(IJK, i_electron) * coeff(1)
       if (tmp < 0) tmp = 0
       box%cc(IJK, i_electron_old) = tmp
    end do; CLOSE_DO
  end subroutine set_photoionization_rate

end module m_photoi_$Dd
