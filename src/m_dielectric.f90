#include "../afivo/src/cpp_macros.h"
!> Module with settings and routines to handle dielectrics
!>
!> @todo Make sure species densities are initially zero inside the dielectric
!> @todo Use special prolongation for multigrid when there are surface charges
module m_dielectric
  use m_af_all
  use m_streamer
  use m_types

  implicit none
  private

  ! @todo maybe don't hardcode this?
  integer, parameter, public :: i_photon_flux = 1
  integer, parameter, public :: i_surf_dens = 2

  !> To store dielectric surface
  type(surfaces_t), public :: diel

  ! Maximum travel distance for testing boundary intersection
  real(dp), protected :: photon_step_length = 1.0e-3_dp

  !> Secondary electron emission coefficient for positive ion impact
  real(dp), protected :: gamma_se_ion = 0.1_dp

  !> Secondary electron emission coefficient for high energy photons
  real(dp), protected :: gamma_se_ph_highenergy = 0.1_dp

  !> Secondary electron emission coefficient for low energy photons
  real(dp), protected :: gamma_se_ph_lowenergy = 0.1_dp

  !> Assume photons are not absorbed for photoemission computation
  logical :: photons_no_absorption = .false.

  !> Preset surface charge numbers
  integer, protected :: n_surface_charge

  !> Preset surface charge
  real(dp), allocatable, protected :: preset_charge(:)

  !> Preset surface charge distribution
  real(dp), allocatable, protected :: preset_charge_distribution(:)

  public :: dielectric_initialize
  public :: dielectric_update_surface_charge
  public :: dielectric_photon_emission
  public :: dielectric_photon_absorption
  public :: dielectric_reset_photons

contains

  subroutine dielectric_initialize(tree, cfg)
    use m_config
    type(af_t), intent(in) :: tree
    type(CFG_t), intent(inout) :: cfg

    call CFG_add_get(cfg, "dielectric%photon_step_length", &
         photon_step_length, &
         "Maximum travel distance for testing boundary intersection")
    call CFG_add_get(cfg, "dielectric%gamma_se_ph_highenergy", &
         gamma_se_ph_highenergy, &
         "Secondary electron emission coefficient for high energy photons")
    call CFG_add_get(cfg, "dielectric%gamma_se_ph_lowenergy", &
         gamma_se_ph_lowenergy, &
         "Secondary electron emission coefficient for low energy photons")
    call CFG_add_get(cfg, "dielectric%gamma_se_ion", &
         gamma_se_ion, &
         "Secondary electron emission coefficient for positive ion impact")
    call CFG_add_get(cfg, "dielectric%photons_no_absorption", &
         photons_no_absorption, &
         "Assume photons are not absorbed for photoemission computation")
    call CFG_add(cfg, "dielectric%preset_charge", [0.0_dp], &
         "preset nonuniform surface charge", dynamic_size=.true.)
    call CFG_add(cfg, "dielectric%preset_charge_distribution", [0.0_dp], &
         "The distribution of nonuniform surface charge", dynamic_size=.true.)
    call CFG_get_size(cfg, "dielectric%preset_charge", n_surface_charge)
    allocate(preset_charge(n_surface_charge))
    allocate(preset_charge_distribution(n_surface_charge))
    call CFG_get(cfg, "dielectric%preset_charge", preset_charge)
    call CFG_get(cfg, "dielectric%preset_charge_distribution", preset_charge_distribution)
    preset_charge_distribution = preset_charge_distribution * ST_domain_len(NDIM)

  end subroutine dielectric_initialize

  subroutine dielectric_update_surface_charge(box, surface, dt, s_in, s_out)
    use m_units_constants
    use m_streamer
    type(box_t), intent(inout)          :: box
    type(surface_t), intent(inout) :: surface
    real(dp), intent(in)                :: dt    !< Time step
    integer, intent(in)                 :: s_in  !< Input state
    integer, intent(in)                 :: s_out !< Output state
#if NDIM == 2
    real(dp)                            :: se_flux(box%n_cell)
#elif NDIM == 3
    real(dp)                            :: se_flux(box%n_cell, box%n_cell)
#endif
    integer                             :: nc
    real(dp)                            :: dr

    nc  = box%n_cell
    dr  = box%dr(af_neighb_dim(surface%direction))

    select case (surface%direction)
#if NDIM == 2
    case (af_neighb_lowx)
       surface%sd(:, i_surf_dens+s_out) = surface%sd(:, i_surf_dens+s_in) - &
            dt * matmul(box%fc(1, 1:nc, 1, flux_variables), &
            flux_species_charge)

       if (size(flux_pos_ion) > 0 .and. gamma_se_ion > 0.0_dp) then
          ! Compute secondary emission flux
          se_flux = -gamma_se_ion * sum(box%fc(1, 1:nc, 1, flux_pos_ion), dim=2)
          box%cc(1, 1:nc, i_electron+s_out) = &
               box%cc(1, 1:nc, i_electron+s_out) + dt * se_flux / dr
          surface%sd(:, i_surf_dens+s_out) = surface%sd(:, i_surf_dens+s_out) + &
               dt * se_flux
       end if
    case (af_neighb_highx)
       surface%sd(:, i_surf_dens+s_out) = surface%sd(:, i_surf_dens+s_in) + &
            dt * matmul(box%fc(nc+1, 1:nc, 1, flux_variables), &
            flux_species_charge)

       if (size(flux_pos_ion) > 0 .and. gamma_se_ion > 0.0_dp) then
          ! Compute secondary emission flux
          se_flux = gamma_se_ion * sum(box%fc(nc+1, 1:nc, 1, flux_pos_ion), dim=2)
          box%cc(nc, 1:nc, i_electron+s_out) = &
               box%cc(nc, 1:nc, i_electron+s_out) + dt * se_flux / dr
          surface%sd(:, i_surf_dens+s_out) = surface%sd(:, i_surf_dens+s_out) + &
               dt * se_flux
       end if
    case (af_neighb_lowy)
       surface%sd(:, i_surf_dens+s_out) = surface%sd(:, i_surf_dens+s_in) - &
            dt * matmul(box%fc(1:nc, 1, 2, flux_variables), &
            flux_species_charge)

       if (size(flux_pos_ion) > 0 .and. gamma_se_ion > 0.0_dp) then
          ! Compute secondary emission flux
          se_flux = -gamma_se_ion * sum(box%fc(1:nc, 1, 2, flux_pos_ion), dim=2)
          box%cc(1:nc, 1, i_electron+s_out) = &
               box%cc(1:nc, 1, i_electron+s_out) + dt * se_flux / dr
          surface%sd(:, i_surf_dens+s_out) = surface%sd(:, i_surf_dens+s_out) + &
               dt * se_flux
       end if
    case (af_neighb_highy)
       surface%sd(:, i_surf_dens+s_out) = surface%sd(:, i_surf_dens+s_in) + &
            dt * matmul(box%fc(1:nc, nc+1, 2, flux_variables), &
            flux_species_charge)

       if (size(flux_pos_ion) > 0 .and. gamma_se_ion > 0.0_dp) then
          ! Compute secondary emission flux
          se_flux = gamma_se_ion * sum(box%fc(1:nc, nc+1, 2, flux_pos_ion), dim=2)
          box%cc(1:nc, nc, i_electron+s_out) = &
               box%cc(1:nc, nc, i_electron+s_out) + dt * se_flux / dr
          surface%sd(:, i_surf_dens+s_out) = surface%sd(:, i_surf_dens+s_out) + &
               dt * se_flux
       end if
#elif NDIM == 3
    case default
       error stop
#endif
    end select
  end subroutine dielectric_update_surface_charge

  subroutine dielectric_photon_emission(box, surface, dt, s_in, s_out)
    use m_units_constants
    use m_streamer
    type(box_t), intent(inout)          :: box
    type(surface_t), intent(inout) :: surface
    real(dp), intent(in)                :: dt    !< Time step
    integer, intent(in)                 :: s_in  !< Input state
    integer, intent(in)                 :: s_out !< Output state
    integer                             :: nc
    real(dp)                            :: dr

    nc  = box%n_cell
    dr  = box%dr(af_neighb_dim(surface%direction))

    select case (surface%direction)
#if NDIM == 2
    case (af_neighb_lowx)
       where (box%fc(1, 1:nc, 1, electric_fld) < 0.0_dp)
          box%cc(1, 1:nc, i_electron+s_out) = &
               box%cc(1, 1:nc, i_electron+s_out) + &
               surface%sd(:, i_photon_flux) * dt / dr
          surface%sd(:, i_surf_dens+s_out) = surface%sd(:, i_surf_dens+s_out) + &
               surface%sd(:, i_photon_flux) * dt * UC_elem_charge
       end where
    case (af_neighb_highx)
       where (box%fc(nc, 1:nc, 1, electric_fld) > 0.0_dp)
          box%cc(nc, 1:nc, i_electron+s_out) = &
               box%cc(nc, 1:nc, i_electron+s_out) + &
               surface%sd(:, i_photon_flux) * dt / dr
          surface%sd(:, i_surf_dens+s_out) = surface%sd(:, i_surf_dens+s_out) + &
               surface%sd(:, i_photon_flux) * dt * UC_elem_charge
       end where
    case (af_neighb_lowy)
       where (box%fc(1:nc, 1, 2, electric_fld) < 0.0_dp)
          box%cc(1:nc, 1, i_electron+s_out) = &
               box%cc(1:nc, 1, i_electron+s_out) + &
               surface%sd(:, i_photon_flux) * dt / dr
          surface%sd(:, i_surf_dens+s_out) = surface%sd(:, i_surf_dens+s_out) + &
               surface%sd(:, i_photon_flux) * dt * UC_elem_charge
       end where
    case (af_neighb_highy)
       where (box%fc(1:nc, nc, 2, electric_fld) > 0.0_dp)
          box%cc(1:nc, nc, i_electron+s_out) = &
               box%cc(1:nc, nc, i_electron+s_out) + &
               surface%sd(:, i_photon_flux) * dt / dr
          surface%sd(:, i_surf_dens+s_out) = surface%sd(:, i_surf_dens+s_out) + &
               surface%sd(:, i_photon_flux) * dt * UC_elem_charge
       end where
#elif NDIM == 3
    case default
       error stop
#endif
    end select

  end subroutine dielectric_photon_emission

  !> Determine whether and where photons hit a dielectric, and change their
  !> absorption location to the first cell inside the surface. If
  !> photons_no_absorption is true, assume that photons are not absorbed by the
  !> gas (so extrapolate their path).
  subroutine dielectric_photon_absorption(tree, xyz_start, xyz_end, n_dim, &
       n_photons, photon_weight)
    use m_af_types
    use m_af_interp
    use m_streamer
    type(af_t), intent(in)  :: tree
    integer, intent(in)     :: n_dim
    integer, intent(in)     :: n_photons
    real(dp), intent(in)    :: xyz_start(3, n_photons)
    real(dp), intent(inout) :: xyz_end(3, n_photons)
    real(dp), intent(in)    :: photon_weight
    real(dp)                :: xyz(n_dim), dvec(n_dim)
    real(dp)                :: dvec_small(n_dim), dvec_large(n_dim)
    real(dp)                :: xyz_gas(n_dim), xyz_nogas(n_dim)
    real(dp)                :: xyz_middle(n_dim), eps(1)
    real(dp)                :: travel_distance
    integer                 :: n, n_steps, n_steps_extra, i, k, n_bisect
    logical                 :: success

    ! Determine the number of bisection steps to find the first cell inside the
    ! dielectric, given a photon step length <= photon_step_length
    n_bisect = -ceiling(&
         log(af_min_dr(tree)/photon_step_length) / log(2.0_dp))

    if (photons_no_absorption) then
       n_steps_extra = ceiling(norm2(ST_domain_len) / photon_step_length)
    else
       n_steps_extra = 0
    end if

    do n = 1, n_photons
       xyz             = xyz_start(1:n_dim, n)
       dvec            = xyz_end(1:n_dim, n) - xyz_start(1:n_dim, n)
       travel_distance = norm2(dvec)

       ! Large photon step length
       dvec_large = (dvec/travel_distance) * photon_step_length

       n_steps    = ceiling(travel_distance/photon_step_length)
       ! Normalize direction vector to right length. Possible TODO: near the
       ! boundary of the domain, a photon can fly out crossing on a small
       ! piece of a dielectric.
       dvec_small = dvec / n_steps

       do i = 1, n_steps + n_steps_extra
          if (i <= n_steps) then
             dvec = dvec_small
          else
             dvec = dvec_large
          end if

          xyz = xyz + dvec

          ! If outside, stop
          if (any(xyz < ST_domain_origin .or. &
               xyz > ST_domain_origin + ST_domain_len)) then
             exit
          end if

          ! Get dielectric permittivity
          eps = af_interp0(tree, xyz, [i_eps], success)
          if (.not. success) error stop "photon unexpectedly outside domain"

          ! If epsilon changes, start doing a search for the intersection point
          if (eps(1) > 1.0_dp) then
             ! Perform bisection to locate first cell inside dielectric
             xyz_gas = xyz - dvec
             xyz_nogas = xyz

             do k = 1, n_bisect
                xyz_middle = 0.5_dp * (xyz_gas + xyz_nogas)
                eps = af_interp0(tree, xyz_middle, [i_eps], success)
                if (.not. success) error stop "photon unexpectedly outside domain"

                if (eps(1) > 1.0_dp) then
                   xyz_nogas = xyz_middle
                else
                   xyz_gas = xyz_middle
                end if
             end do

             if (i <= n_steps) then
                ! The photon was absorbed within its normal travel path
                xyz_end(1:n_dim, n) = -1e50_dp
                call add_to_surface_photons(tree, xyz_nogas, photon_weight, &
                     gamma_se_ph_highenergy)
             end if
             call add_to_surface_photons(tree, xyz_nogas, photon_weight, &
                  gamma_se_ph_lowenergy)
             exit
          end if
       end do
    end do
  end subroutine dielectric_photon_absorption

  subroutine add_to_surface_photons(tree, xyz, w, frac)
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: xyz(NDIM)
    real(dp), intent(in)   :: w
    real(dp), intent(in)   :: frac
    integer                :: ix, ix_cell(NDIM-1)
    real(dp)               :: area

    call surface_get_surface_cell(tree, diel, xyz, ix, ix_cell)
    area = product(diel%surfaces(ix)%dr)

#if NDIM == 2
    diel%surfaces(ix)%sd(ix_cell(1), i_photon_flux) = &
         diel%surfaces(ix)%sd(ix_cell(1), i_photon_flux) &
         + frac * w / area
#else
    error stop
#endif
  end subroutine add_to_surface_photons

  subroutine dielectric_reset_photons()
    integer :: n

    do n = 1, diel%max_ix
       if (diel%surfaces(n)%in_use) then
#if NDIM == 1
          diel%surfaces(n)%sd(i_photon_flux) = 0.0_dp
#elif NDIM == 2
          diel%surfaces(n)%sd(:, i_photon_flux) = 0.0_dp
#elif NDIM == 3
          diel%surfaces(n)%sd(:, :, i_photon_flux) = 0.0_dp
#endif
       end if
    end do
  end subroutine dielectric_reset_photons

end module m_dielectric
