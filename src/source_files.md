# Description of source files

## The src directory

  * m_geometry.f90: Module used to describe initial conditions (lines, Gaussians, etc.)
  * m_photons.f90: Module for photoionization
  * m_streamer.f90: Module with common functionality of 2D/3D streamer models
  * m_transport_data.f90: Module to read in transport data
  * m_units_constants.f90: small list of units/constants (physics)

## Helper-projects included in subdirectories

  * config_fortran: to read in configuration files
  * lookup_table_fortran: create/read from linear lookup tables, for efficient lookup of values
  * rng_fortran: generate different types of random numbers.
