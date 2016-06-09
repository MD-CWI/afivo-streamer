streamer tools
==============

originally the files mentioned below were included in the fosito directory.
Now they are moved to 'src' and together with the files

  * m_geometry.f90:
  * m_photons.f90:
  * m_streamer.f90:
  * m_find_index.f90:
  * m_linked_list.f90:
  * m_mrgrnk.f90:
  * m_transport_data.f90:

This is a collection of tools that can help doing simulations:

  * m_config.f90: use text based configuration files, whose values can be used in a simulation
  * m_lookup_table.f90: create/read from linear lookup tables, for efficient lookup of values
  * m_units_constants.f90: small list of units/constants (physics)
  * m_random.f90: using the standard random number generator, generate different types of random numbers.

Both collections are included in libstreamer.a

Usage: type 'make' to compile a static library, or simply use the source files directly.
There are a few basic tests included.
