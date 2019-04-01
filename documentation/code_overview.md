# Code overview

# Features

* Adaptive mesh refinement
* OpenMP parallelization
* Multigrid field solvers
* Plasma fluid model for streamer simulations
* Photoionization using a Helmholtz or Monte Carlo approach, see @ref md_documentation_photoionization
* Support for tabulated transport data, see @ref md_documentation_transport_data
* Support for user-defined chemical reactions, see @ref md_documentation_chemistry
* Flexible usage through configuration files and command line arguments, see @ref md_documentation_simulation_options

For more details see the publications describing afivo-streamer at @ref md_documentation_publications.

# Source code structure

Folder | Contents
---|---
src | The source code of afivo-streamer
src/config_fortran | Module for configuration files
src/lookup_table_fortran | Module for lookup tables
src/rng_fortran | Module for random numbers
afivo | The afivo library
documentation | The documentation (mostly Markdown files)
lib_2d | The afivo-streamer library (2D version)
lib_3d | The afivo-streamer library (3D version)
makefiles | Makefiles for building the library and programs
programs | Streamer simulation programs
tools | Python scripts for e.g. converting data
transport_data | Input data files for simulations





