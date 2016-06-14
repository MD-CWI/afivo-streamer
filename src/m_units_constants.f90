!> Module that contains physical and numerical constants
module m_units_constants
  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  ! Numerical constants
  real(dp), parameter :: UC_pi          = acos(-1.0_dp)

  ! Physical constants
  ! TODO: add sign. digits
  real(dp), parameter :: UC_eps0             = 8.8541878176d-12        ! permitivity of vacuum (SI)
  real(dp), parameter :: UC_elem_charge      = 1.6022d-19              ! the elementary charge in Coulombs
  real(dp), parameter :: UC_elec_charge      = -1.6022d-19             ! the electron charge in Coulombs
  real(dp), parameter :: UC_elec_volt        = 1.6022d-19              ! the eV in joules
  real(dp), parameter :: UC_elec_mass        = 9.10938189d-31          ! the electron mass in kg
  real(dp), parameter :: UC_atomic_mass      = 1.66053886D-27          ! the atomic mass unit in kg
  real(dp), parameter :: UC_N2_mass          = 28.0D0 * UC_atomic_mass ! The mass of a N2 molecule
  real(dp), parameter :: UC_lightspeed       = 299792458d0             ! the speed of light in m/s
  real(dp), parameter :: UC_boltzmann_const  = 1.3806503d-23           ! the Boltzmann constant
  real(dp), parameter :: UC_bohr_radius      = 5.29d-11                ! the Bohr radius (m)
  real(dp), parameter :: UC_torr_to_bar      = 133.322368 * 1.0D-5     ! one Torr in units of bar
  real(dp), parameter :: UC_elec_q_over_eps0 = UC_elec_charge / UC_eps0
  real(dp), parameter :: UC_elec_q_over_m    = UC_elec_charge / UC_elec_mass

end module m_units_constants
