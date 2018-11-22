module m_gas
  use m_types

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  !> Whether the gas has a constant density
  logical, public, protected :: gas_constant_density = .true.

  ! Pressure of the gas in bar
  real(dp), public, protected :: gas_pressure = 1.0_dp

  ! Gas temperature in Kelvin
  real(dp), public, protected :: gas_temperature = 300.0_dp

  real(dp), allocatable, public, protected :: gas_fractions(:)
  real(dp), allocatable, public, protected :: gas_densities(:)
  character(len=comp_len), allocatable, public, protected :: gas_components(:)

  ! Gas number density (1/m3)
  real(dp), public, protected :: gas_number_density

  ! Convert V/m to Townsend
  real(dp), public, protected :: SI_to_Townsend

  ! Convert Townsend to V/m
  real(dp), public, protected :: Townsend_to_SI

  public :: gas_init
  public :: gas_index

contains

  subroutine gas_init(cfg)
    use m_config
    use m_units_constants
    type(CFG_t), intent(inout) :: cfg
    integer                    :: n

    call CFG_add_get(cfg, "gas%constant_density", gas_constant_density, &
         "Whether the density is constant")
    call CFG_add_get(cfg, "gas%pressure", gas_pressure, &
         "The gas pressure (bar)")
    call CFG_add_get(cfg, "gas%temperature", gas_temperature, &
         "The gas temperature (Kelvin)")

    ! Ideal gas law
    gas_number_density = 1e5_dp * gas_pressure / &
         (UC_boltzmann_const * gas_temperature)
    SI_to_Townsend = 1e21_dp / gas_number_density
    Townsend_to_SI = 1.0_dp / SI_to_Townsend

    call CFG_add(cfg, "gas%components", ["N2", "O2"], &
         "Gas component names", .true.)
    call CFG_add(cfg, "gas%fractions", [0.8_dp, 0.2_dp], &
         "Gas component fractions", .true.)

    call CFG_get_size(cfg, "gas%components", n)
    allocate(gas_components(n))
    allocate(gas_fractions(n))
    call CFG_get(cfg, "gas%components", gas_components)
    call CFG_get(cfg, "gas%fractions", gas_fractions)

    if (any(gas_fractions < 0.0_dp)) &
         error stop "gas%fractions has negative value"
    if (abs(sum(gas_fractions) - 1.0_dp) > 1e-4_dp) &
         error stop "gas%fractions not normalized"

    gas_densities = gas_fractions * gas_number_density

  end subroutine gas_init

  !> Find index of a gas component, return -1 if not found
  elemental integer function gas_index(name)
    character(len=*), intent(in) :: name
    do gas_index = 1, size(gas_components)
       if (gas_components(gas_index) == name) exit
    end do
    if (gas_index == size(gas_components)+1) gas_index = -1
  end function gas_index

end module m_gas
