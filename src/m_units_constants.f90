!> Module that contains physical and numerical constants
module m_units_constants
  implicit none
  public

  integer, parameter, private :: dp          = kind(0.0d0)

  ! Numerical constants
  real(dp), parameter :: UC_pi               = acos(-1.0_dp)

  ! Permitivity of vacuum (SI)
  real(dp), parameter :: UC_eps0             = 8.8541878128e-12_dp

  ! The elementary charge in Coulombs
  real(dp), parameter :: UC_elem_charge      = 1.602176634e-19_dp

  ! The electron charge in Coulombs
  real(dp), parameter :: UC_elec_charge      = -UC_elem_charge

  ! The eV in joules
  real(dp), parameter :: UC_elec_volt        = UC_elem_charge

  ! The electron mass in kg
  real(dp), parameter :: UC_elec_mass        = 9.1093837015e-31_dp

  ! The atomic mass unit in kg
  real(dp), parameter :: UC_atomic_mass      = 1.6605390666e-27_dp

  ! The mass of a N2 molecule
  real(dp), parameter :: UC_N2_mass          = 28.0D0 * UC_atomic_mass

  ! The mass of an O2 molecule
  real(dp), parameter :: UC_O2_mass          = 32.0D0 * UC_atomic_mass

  ! The speed of light in m/s
  real(dp), parameter :: UC_lightspeed       = 299792458.0_dp

  ! The Boltzmann constant in J/K
  real(dp), parameter :: UC_boltzmann_const  = 1.380649e-23_dp

  ! The Bohr radius (m)
  real(dp), parameter :: UC_bohr_radius      = 5.29177210903e-11_dp

  ! One Torr in units of bar
  real(dp), parameter :: UC_torr_to_bar      = 133.32236842105263_dp * 1.0e-5_dp

  ! Electron charge over epsilon0
  real(dp), parameter :: UC_elec_q_over_eps0 = UC_elec_charge / UC_eps0

  ! Electron charge over electron mass
  real(dp), parameter :: UC_elec_q_over_m    = UC_elec_charge / UC_elec_mass

end module m_units_constants
