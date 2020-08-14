#include "../afivo/src/cpp_macros.h"
!> Module to perform simulations with an external circuit
module m_circuit
  use m_types
  use m_af_all

  implicit none
  private

  !> Whether an external circuit is used
  logical, public, protected :: circuit_used = .false.

  !> Type of external circuit
  character(len=name_len) :: circuit_type = undefined_str

  !> Resistances (Ohm)
  real(dp), allocatable :: resistors(:)

  !> Capacities (farad)
  real(dp), allocatable :: capacitors(:)

  !> Initial charge on capacitors (Coulomb)
  real(dp), allocatable :: capacitors_q(:)

  public :: circuit_initialize
  public :: circuit_update

contains

  !> Initialize this module
  subroutine circuit_initialize(tree, cfg, restart)
    use m_config
    use m_units_constants
    use m_field
    use m_streamer
    type(af_t), intent(in)     :: tree
    type(CFG_t), intent(inout) :: cfg
    logical, intent(in)        :: restart
    real(dp)                   :: voltage
    logical                    :: V0_from_field = .true.

    call CFG_add_get(cfg, "circuit%type", circuit_type, &
         "Type of external circuit")
    call CFG_add(cfg, "circuit%resistors", [3.0e2_dp], &
         "Resistances (Ohm)", .true.)
    call CFG_add(cfg, "circuit%capacitors", [0.2e-9_dp], &
         "Capacities (farad)", .true.)
    call CFG_add(cfg, "circuit%capacitors_q", [0.0_dp], &
         "Initial charge on capacitors (Coulomb)", .true.)
    call CFG_add_get(cfg, "circuit%V0_from_field", V0_from_field, &
         "Get initial voltage (and capacitor charge) from applied field")

    if (restart .and. circuit_type /= undefined_str) then
       error stop "TODO: Circuit does not support restarting"
    end if

    select case (circuit_type)
    case (undefined_str)
       return
    case ("RC")
       circuit_used = .true.
       allocate(resistors(1), capacitors(1), capacitors_q(1))
       call CFG_get(cfg, "circuit%resistors", resistors)
       call CFG_get(cfg, "circuit%capacitors", capacitors)
       call CFG_get(cfg, "circuit%capacitors_q", capacitors_q)

       if (V0_from_field) then
          voltage = -ST_domain_len(NDIM) * &
               field_get_amplitude(tree, 0.0_dp)
          ! Set initial charge from voltage
          capacitors_q(1) = capacitors(1) * voltage
       else
          ! Set initial voltage from capacitor charge
          voltage = capacitors_q(1) / capacitors(1)
       end if
       print *, "initial voltage", voltage
       call field_set_voltage_externally(voltage)
    case default
       error stop "Unknown circuit type"
    end select

  end subroutine circuit_initialize

    !> Add source terms form the fluid model to the Euler equations
  subroutine circuit_update(tree, dt)
    use m_units_constants
    use m_field
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt
    real(dp)                  :: sum_J_E0, current, new_voltage

    ! Compute the volume integral of J dot E0, where J is the electron flux and
    ! E0 is the electric field without space charge.
    ! We then have V * I = -e * sum_J_E0 [Sato, J. Phys. D., 1980]
    call af_reduction(tree, get_J_dot_E0_box, my_sum, 0.0_dp, sum_J_E0)

    select case (circuit_type)
    case ("RC")
       ! Compute current for the old voltage. With small time steps, this simple
       ! approximation is fine. Energy is not strictly conserved, however.
       current = (UC_elec_charge * sum_J_E0)/field_voltage

       ! Adjust charge on capacitor
       capacitors_q(1) = capacitors_q(1) - dt * current

       ! Assume the current through the resistor stays the same at the next time
       ! step. Iterating to get the response from the discharge to the new
       ! voltage is unnecessary, since we use small time steps anyway.
       new_voltage = capacitors_q(1) / capacitors(1) - current * resistors(1)
    case default
       error stop "Unknown circuit type"
    end select

    call field_set_voltage_externally(new_voltage)
  end subroutine circuit_update

  !> Compute sum of J dot E0, where J is the electron flux and E0 is the
  !> electric field without space charge
  function get_J_dot_E0_box(box) result(J_dot_E)
    use m_field, only: field_voltage
    use m_streamer
    type(box_t), intent(in) :: box
    integer                 :: IJK, nc
    real(dp)                :: fac, J_dot_E, field(NDIM)

    ! Assume we have plate-plate electrodes and a homogeneous field without
    ! space charge. TODO: generalize this!
    field(:)    = 0.0_dp
    field(NDIM) = -field_voltage / ST_domain_len(NDIM)

    J_dot_E     = 0.0_dp
    fac         = product(box%dr)
    nc          = box%n_cell

    ! Sum inner product flux * field over the cell faces
#if NDIM == 1
    do KJI_DO(1, nc)
       J_dot_E = J_dot_E + fac * 0.5_dp * (&
            sum(box%fc(IJK, :, flux_elec) * field) + &
            box%fc(i+1, 1, flux_elec) * field(1))
    end do; CLOSE_DO
#elif NDIM == 2
    if (ST_cylindrical) then
       ! Multiply values with 2 * pi * r
       fac = fac * 2 * acos(-1.0_dp)
       do KJI_DO(1, nc)
          J_dot_E = J_dot_E + af_cyl_radius_cc(box, i) * &
               fac * 0.5_dp * (&
               sum(box%fc(IJK, :, flux_elec) * field) + &
               box%fc(i+1, j, 1, flux_elec) * field(1) + &
               box%fc(i, j+1, 2, flux_elec) * field(2))
       end do; CLOSE_DO
    else
       do KJI_DO(1, nc)
          J_dot_E = J_dot_E + fac * 0.5_dp * (&
               sum(box%fc(IJK, :, flux_elec) * field) + &
               box%fc(i+1, j, 1, flux_elec) * field(1) + &
               box%fc(i, j+1, 2, flux_elec) * field(2))
       end do; CLOSE_DO
    end if
#elif NDIM == 3
    do KJI_DO(1, nc)
       J_dot_E = J_dot_E + fac * 0.5_dp * (&
            sum(box%fc(IJK, :, flux_elec) * field) + &
            box%fc(i+1, j, k, 1, flux_elec) * field(1) + &
            box%fc(i, j+1, k, 2, flux_elec) * field(2) + &
            box%fc(i, j, k+1, 3, flux_elec) * field(3))
    end do; CLOSE_DO
#endif

  end function get_J_dot_E0_box

  real(dp) function my_sum(a, b)
    real(dp), intent(in) :: a, b
    my_sum = a + b
  end function my_sum

end module m_circuit
