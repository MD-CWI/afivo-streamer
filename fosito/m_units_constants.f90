! Copyright 2005-2012, Chao Li, Margreet Nool, Anbang Sun, Jannis Teunissen
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

!> Module that contains physical and numerical constants
module m_units_constants
  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  ! Numerical constants
  real(dp), parameter :: UC_pi              = acos(-1.0_dp)
  real(dp), parameter :: UC_sqrt2           = sqrt(2.0_dp)
  real(dp), parameter :: UC_iSqrt2          = 1 / UC_sqrt2
  real(dp), parameter :: UC_one_third       = 1 / 3.0_dp
  real(dp), parameter :: UC_one_sixth       = 1 / 6.0_dp
  real(dp), parameter :: UC_one_seventh     = 1 / 7.0_dp
  real(dp), parameter :: UC_one_ninth       = 1 / 9.0_dp

  ! Physical constants
  ! TODO: add sign. digits
  real(dp), parameter :: UC_eps0            = 8.8541878176d-12        ! permitivity of vacuum (SI)
  real(dp), parameter :: UC_elem_charge     = 1.6022d-19             ! the elementary charge in Coulombs
  real(dp), parameter :: UC_elec_charge     = -1.6022d-19             ! the electron charge in Coulombs
  real(dp), parameter :: UC_elec_volt       = 1.6022d-19              ! the eV in joules
  real(dp), parameter :: UC_elec_mass       = 9.10938189d-31          ! the electron mass in kg
  real(dp), parameter :: UC_atomic_mass     = 1.66053886D-27          ! the atomic mass unit in kg
  real(dp), parameter :: UC_N2_mass         = 28.0D0 * UC_atomic_mass ! The mass of a N2 molecule
  real(dp), parameter :: UC_lightspeed      = 299792458d0             ! the speed of light in m/s
  real(dp), parameter :: UC_boltzmann_const = 1.3806503d-23           ! the Boltzmann constant
  real(dp), parameter :: UC_bohr_radius     = 5.29d-11                ! the Bohr radius (m)
  real(dp), parameter :: UC_torr_to_bar     = 133.322368 * 1.0D-5     ! one Torr in units of bar
  real(dp), parameter :: UC_elec_q_over_eps0 = UC_elec_charge / UC_eps0
  real(dp), parameter :: UC_elec_q_over_m = UC_elec_charge / UC_elec_mass

  ! Small and large numbers
  real(dp), parameter :: UC_tiny            = epsilon(1.0_dp)
  real(dp), parameter :: UC_huge            = huge(1.0_dp)

contains

                                                                        ! Convert (classical) kinetic energy to velocity
  real(dp) function UC_en_to_vel(en, mass)
    real(dp), intent(in) :: en, mass
    UC_en_to_vel                            = sqrt(2 * en) / mass
  end function UC_en_to_vel

  ! Convert velocity to energy
  real(dp) function UC_vel_to_en(vel, mass)
    real(dp), intent(in) :: vel, mass
    UC_vel_to_en = 0.5_dp * mass * vel**2
  end function UC_vel_to_en

  subroutine UC_xyz_to_spherical(xyz, radius, theta, phi)
    real(dp), intent(in) :: xyz(3)
    real(dp), intent(out) :: radius, theta, phi

    if (xyz(1) == 0.0_dp .and. xyz(2) == 0.0_dp) then
       phi = 0.0_dp
    else
       phi = atan2(xyz(2), xyz(1))
    end if

    radius = sqrt(sum(xyz**2))
    if (radius == 0.0_dp) then      ! Undefined
       theta = acos(xyz(3) / radius)
    else
       theta = 0.0_dp
    end if

 end subroutine UC_xyz_to_spherical

end module m_units_constants
