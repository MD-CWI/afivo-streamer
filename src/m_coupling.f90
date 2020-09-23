#include "../afivo/src/cpp_macros.h"
!> Module for the coupling between the gas dynamics and the fluid model
module m_coupling
  use m_types
  use m_af_all
  use m_units_constants
  use m_gas
  use m_streamer

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
    real(dp)                   :: J_dot_E

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

       box%cc(IJK, gas_vars(i_e)) = box%cc(IJK, gas_vars(i_e)) + &
           gas_heating_efficiency *  J_dot_E * UC_elec_charge * dt_vec(1)
    end do; CLOSE_DO
  end subroutine add_heating_box

  !> Update gas number density in the fluid model
  subroutine coupling_update_gas_density(tree)
    type(af_t), intent(inout) :: tree

    call af_loop_box(tree, update_gas_density, .true.)
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
