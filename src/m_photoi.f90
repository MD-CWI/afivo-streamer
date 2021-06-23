#include "../afivo/src/cpp_macros.h"
!> Top-module for photoionization, which can make use of different methods
module m_photoi
  use m_photoi_mc
  use m_photoi_helmh
  use m_af_all
  use m_streamer
  use m_types

  implicit none
  private

  ! Whether photoionization is enabled
  logical, protected, public :: photoi_enabled = .false.

  ! Whether photoemission is enabled
  logical, protected, public :: photoe_enabled = .false.

  ! Whether photoionization is enabled in gas
  logical, protected, public :: photoi_enabled_ingas = .true.

  ! Which photoionization method to use (helmholtz, montecarlo)
  character(len=20) :: photoi_method = 'helmholtz'

  ! Which photoionization source to use (Zheleznyak, from_species)
  character(len=20) :: photoi_source_type = 'Zheleznyak'

  ! Name of the excited species in case photoi_source_type is 'from_species'
  character(len=string_len) :: photoi_excited_species = 'UNDEFINED'

  ! Index of the excited species
  integer :: i_excited_species = -1

  ! Photoionization efficiency factor, typically around 0.05-0.1, not for Helmholtz-Luque should be 1.0
  real(dp) :: photoi_eta = 0.05_dp

  ! Quenching pressure
  real(dp) :: photoi_quenching_pressure = 40e-3 ! mbar

  ! Decay time in case photoi_source_type is 'from_species'
  real(dp) :: photoi_photoemission_time = 0.0_dp

  ! Update photoionization every N time step
  integer, protected, public :: photoi_per_steps = 5

  ! Update photoemission every N time step
  integer, protected, public :: photoe_per_steps = 10

  ! Optional variable (when using photoionization)
  integer, public, protected :: i_photo = -1 ! Photoionization rate

  public :: photoi_initialize
  public :: photoi_set_src
  ! public :: photoe_set_src

contains

  !> Initialize photoionization parameters
  subroutine photoi_initialize(tree, cfg)
    use m_config
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg !< The configuration for the simulation

    call CFG_add_get(cfg, "photoi%enabled", photoi_enabled, &
         "Whether photoionization is enabled")
    call CFG_add_get(cfg, "photoi%enabled_ingas", photoi_enabled_ingas, &
         "Whether photoionization is enabled in gas")
    call CFG_add_get(cfg, "photoi%per_steps", photoi_per_steps, &
         "Update photoionization every N time step")
    call CFG_add_get(cfg, "photoi%method", photoi_method, &
         "Which photoionization method to use (helmholtz, montecarlo)")
    call CFG_add_get(cfg, "photoi%eta", photoi_eta, &
         "Photoionization efficiency factor, typically around 0.05-0.1")
    call CFG_add_get(cfg, "photoi%quenching_pressure", photoi_quenching_pressure, &
         "Photoionization quenching pressure (bar)")
    call CFG_add_get(cfg, "photoi%source_type", photoi_source_type, &
         "How to compute the photoi. source (Zheleznyak, from_species)")
    call CFG_add_get(cfg, "photoi%excited_species", photoi_excited_species, &
         "Which excited species to use when photoi%source_type = from_species")
    call CFG_add_get(cfg, "photoi%photoemission_time", photoi_photoemission_time, &
         "Photoemission time delay in case photoi_source_type is 'from_species'")

    call CFG_add_get(cfg, "photoe%enabled", photoe_enabled, &
         "Whether photoemission is enabled")
    call CFG_add_get(cfg, "photoe%per_steps", photoe_per_steps, &
         "Update photoemission every N time step")

    if (photoi_enabled) then
       if (photoi_eta <= 0.0_dp) error stop "photoi%eta <= 0.0"
       if (photoi_eta > 1.0_dp) error stop "photoi%eta > 1.0"
       if (photoi_quenching_pressure <= 0.0_dp) &
            error stop "photoi%quenching_pressure <= 0.0"

       call af_add_cc_variable(tree, "photo", ix=i_photo)
       call af_set_cc_methods(tree, i_photo, photoi_helmh_bc, af_gc_interp)

       select case (photoi_source_type)
       case ('Zheleznyak')
          if (photoi_excited_species /= undefined_str) &
             error stop "Cannot use photoi%excited_species with Zheleznyak"
       case ('from_species')
          if (photoi_excited_species == undefined_str) &
             error stop "Please define photoi%excited_species"

          i_excited_species = &
               af_find_cc_variable(tree, trim(photoi_excited_species))
       case default
          error stop "Unknown value for photoi%source_type"
       end select
    end if

    select case (photoi_method)
       case ("helmholtz")
          call photoi_helmh_initialize(tree, cfg, photoi_enabled, photoi_eta)
          call phmc_initialize(cfg, .false.)
       case ("montecarlo")
          call photoi_helmh_initialize(tree, cfg, .false., photoi_eta)
          call phmc_initialize(cfg, photoi_enabled)
       case default
          print *, "Unknown photoi_method: ", trim(photoi_method)
          error stop
    end select

  end subroutine photoi_initialize

  !> Sets the photoionization
  subroutine photoi_set_src(tree, dt)
    use m_units_constants
    use m_gas

    type(af_t), intent(inout)     :: tree
    real(dp), intent(in), optional :: dt
    real(dp)                       :: quench_fac, decay_fraction, eff_decay_time, decay_rate

    ! Compute quench factor, because some excited species will be quenched by
    ! collisions, preventing the emission of a UV photon
    quench_fac = photoi_quenching_pressure / &
         (gas_pressure + photoi_quenching_pressure)

    ! Set photon production rate per cell, which is proportional to the
    ! ionization rate.
    select case (photoi_source_type)
    case ('Zheleznyak')
       call af_loop_box_arg(tree, photoionization_rate_from_alpha, &
            [photoi_eta * quench_fac], .true.)
    case ('from_species')
       ! effective decay time when we have quenching
       eff_decay_time = photoi_photoemission_time * quench_fac
       decay_fraction = 1 - exp(-dt / eff_decay_time)
       if (dt > 1e-6_dp * eff_decay_time) then
          decay_rate = decay_fraction / dt
       else
          decay_rate = 1.0 / eff_decay_time
       end if

       call af_loop_box_arg(tree, photoionization_rate_from_species, &
            [ quench_fac * decay_rate, 1-decay_fraction], .true.)
    end select

    select case (photoi_method)
    case ("helmholtz")
       ! Use Helmholtz approximation
       call photoi_helmh_compute(tree, i_photo)
    case ("montecarlo")
       if (phmc_physical_photons) then
          call phmc_set_src(tree, ST_rng, i_rhs, &
               i_photo, ST_cylindrical, dt)
       else
          call phmc_set_src(tree, ST_rng, i_rhs, &
               i_photo, ST_cylindrical)
       end if
    end select

  end subroutine photoi_set_src

  ! !> Sets the photoemission
  ! subroutine photoe_set_src(tree, dt)
  !   use m_units_constants
  !   use m_gas

  !   type(af_t), intent(inout)     :: tree
  !   real(dp), intent(in), optional :: dt
  !   real(dp), parameter            :: p_quench = 40.0e-3_dp
  !   real(dp)                       :: quench_fac, decay_fraction, decay_rate

  !   ! Compute quench factor, because some excited species will be quenched by
  !   ! collisions, preventing the emission of a UV photon
  !   quench_fac = p_quench / (gas_pressure + p_quench)

  !   ! Set photon production rate per cell, which is proportional to the
  !   ! ionization rate.
  !   select case (photoi_source_type)
  !   case ('Zheleznyak')
  !      call af_loop_box_arg(tree, photoionization_rate_from_alpha, &
  !           [photoi_eta * quench_fac], .true.)
  !   case ('from_species')
  !      decay_fraction = 1 - exp(-dt / eff_decay_time)

  !      if (dt > 1e-6_dp * photoi_decay_time) then
  !         decay_rate = decay_fraction / dt
  !      else
  !         decay_rate = 1 / photoi_decay_time
  !      end if

  !      call af_loop_box_arg(tree, photoionization_rate_from_species, &
  !           [photoi_eta * quench_fac * decay_rate, 1-decay_fraction], .true.)
  !   end select

  !      if (phmc_physical_photons) then
  !         call phe_mc_set_src(tree, ST_rng, i_rhs, &
  !              i_photo, ST_cylindrical, dt)
  !      else
  !         call phe_mc_set_src(tree, ST_rng, i_rhs, &
  !              i_photo, ST_cylindrical)
  !      end if

  ! end subroutine photoe_set_src

  !> Sets the photoionization_rate
  subroutine photoionization_rate_from_alpha(box, coeff)
    use m_streamer
    use m_lookup_table
    use m_transport_data
    use m_gas
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: coeff(:)
    integer                    :: IJK, nc
    real(dp)                   :: fld, Td, alpha, mobility, tmp, gas_dens
    type(LT_loc_t)             :: loc

    nc = box%n_cell

    do KJI_DO(1,nc)
       fld      = box%cc(IJK, i_electric_fld)

       if (gas_constant_density) then
          gas_dens = gas_number_density
       else
          gas_dens = box%cc(IJK, i_gas_dens)
       end if

       Td       = fld * SI_to_Townsend / gas_dens
       loc      = LT_get_loc(td_tbl, Td)
       ! No normalization required: (alpha/N) * (mu/N)
       alpha    = LT_get_col_at_loc(td_tbl, td_alpha, loc)
       mobility = LT_get_col_at_loc(td_tbl, td_mobility, loc)

       tmp = fld * mobility * alpha * box%cc(IJK, i_electron) * coeff(1)
       if (tmp < 0) tmp = 0
       box%cc(IJK, i_rhs) = tmp
    end do; CLOSE_DO
  end subroutine photoionization_rate_from_alpha

  !> Sets the photoionization_rate
  subroutine photoionization_rate_from_species(box, coeff)
    use m_streamer
    use m_lookup_table
    use m_transport_data
    use m_gas
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: coeff(:)
    integer                    :: nc
    real(dp)                   :: source_factor, decay_factor

    nc            = box%n_cell
    source_factor = coeff(1)
    decay_factor  = coeff(2)

    box%cc(DTIMES(1:nc), i_rhs) = source_factor * &
         box%cc(DTIMES(1:nc), i_excited_species)
    box%cc(DTIMES(1:nc), i_excited_species) = decay_factor * &
         box%cc(DTIMES(1:nc), i_excited_species)
  end subroutine photoionization_rate_from_species

end module m_photoi
