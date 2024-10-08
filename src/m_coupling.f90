!> Module for the coupling between the gas dynamics and the fluid model
module m_coupling
#include "../afivo/src/cpp_macros.h"
  use m_types
  use m_af_all
  use m_units_constants
  use m_gas
  use m_streamer
  use m_chemistry
  use m_field

  implicit none
  private

  public :: coupling_add_fluid_source
  public :: coupling_update_gas_density

contains

  !> Add source terms form the fluid model to the Euler equations
  subroutine coupling_add_fluid_source(tree, dt)
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt

    call af_loop_box_arg(tree, add_heating_box, [dt], .true.)
  end subroutine coupling_add_fluid_source

  subroutine add_heating_box(box, dt_vec)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: dt_vec(:)
    integer                    :: IJK, nc
    real(dp)                   :: J_dot_E, ehd_force(NDIM)
    real(dp)                   :: E_vt_release, tmp, eff_fast, eff_slow
    real(dp)                   :: E_vector(DTIMES(1:box%n_cell), NDIM)

    nc = box%n_cell
    do KJI_DO(1, nc)
       ! Compute inner product flux * field over the cell faces
       J_dot_E = 0.5_dp * sum(box%fc(IJK, :, flux_elec) * box%fc(IJK, :, electric_fld))
#if NDIM == 1
       J_dot_E = J_dot_E + 0.5_dp * (&
            box%fc(i+1, 1, flux_elec) * box%fc(i+1, 1, electric_fld))
#elif NDIM == 2
       J_dot_E = J_dot_E + 0.5_dp * (&
            box%fc(i+1, j, 1, flux_elec) * box%fc(i+1, j, 1, electric_fld) + &
            box%fc(i, j+1, 2, flux_elec) * box%fc(i, j+1, 2, electric_fld))
#elif NDIM == 3
       J_dot_E = J_dot_E + 0.5_dp * (&
            box%fc(i+1, j, k, 1, flux_elec) * box%fc(i+1, j, k, 1, electric_fld) + &
            box%fc(i, j+1, k, 2, flux_elec) * box%fc(i, j+1, k, 2, electric_fld) + &
            box%fc(i, j, k+1, 3, flux_elec) * box%fc(i, j, k+1, 3, electric_fld))
#endif

       tmp = J_dot_E * UC_elec_charge * dt_vec(1)

       if (gas_fraction_slow_heating > 0) then
          eff_fast = gas_heating_efficiency * (1 - gas_fraction_slow_heating)
          eff_slow = gas_heating_efficiency * gas_fraction_slow_heating

          E_vt_release = box%cc(IJK, i_vibration_energy)/gas_vt_time * dt_vec(1)
          box%cc(IJK, i_vibration_energy) = box%cc(IJK, i_vibration_energy) + &
               eff_slow * tmp - E_vt_release
          box%cc(IJK, gas_vars(i_e)) = box%cc(IJK, gas_vars(i_e)) + &
               eff_fast * tmp + E_vt_release
       else
          box%cc(IJK, gas_vars(i_e)) = box%cc(IJK, gas_vars(i_e)) + &
               gas_heating_efficiency * tmp
       end if
    end do; CLOSE_DO

    E_vector = field_get_E_vector(box)

    ! EHD force
    do KJI_DO(1, nc)
       ehd_force = UC_elem_charge * sum(box%cc(IJK, charged_species_itree) * &
            charged_species_charge) * E_vector(IJK, :)

       box%cc(IJK, gas_vars(i_mom)) = box%cc(IJK, gas_vars(i_mom)) + &
            gas_EHD_factor * ehd_force * dt_vec(1)
    end do; CLOSE_DO

  end subroutine add_heating_box

  !> Update gas number density in the fluid model
  subroutine coupling_update_gas_density(tree)
    type(af_t), intent(inout) :: tree

    call af_loop_box(tree, update_gas_density, .true.)
    call af_gc_tree(tree, [i_gas_dens])
  end subroutine coupling_update_gas_density

  subroutine update_gas_density(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc
    real(dp)                   :: inv_weight

    nc         = box%n_cell
    inv_weight = 1/gas_molecular_weight

    box%cc(DTIMES(1:nc), i_gas_dens) = &
         box%cc(DTIMES(1:nc), gas_vars(i_rho)) * inv_weight
  end subroutine update_gas_density

end module m_coupling
