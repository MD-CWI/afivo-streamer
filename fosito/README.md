fosito
======

This is a collection of tools that can help doing simulations:

  * m_config.f90: use text based configuration files, whose values can be used in a simulation
  * m_lookup_table.f90: create/read from linear lookup tables, for efficient lookup of values
  * m_units_constants.f90: small list of units/constants (physics)
  * m_random.f90: using the standard random number generator, generate different types of random numbers.

Usage: type 'make' to compile a static library, or simply use the source files
directly. There are a few basic tests included.