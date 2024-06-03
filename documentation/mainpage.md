# Afivo-streamer

[TOC]

Afivo-streamer is a tool for the simulation of streamer discharges in 1D, 2D and 3D,
using plasma fluid models. Afivo-streamer makes uses of the [afivo
framework](https://github.com/MD-CWI/afivo), which for example provides
adaptive mesh refinement and parallelization.

# Recommended reading for new users

* [Installation instructions](documentation/installation.md)
* [Config files and command-line arguments](documentation/simulation_options.md)
* [Saving output and visualization](documentation/output_and_visualization.md)
* [Electrodes and boundary conditions](documentation/electrodes_bc.md)
* [Initial conditions](documentation/initial_conditions.md)
* Any of the tutorials listed on <a href="pages.html">the overview of all documentation</a>

# Git repository links

* [git repository](https://github.com/MD-CWI/afivo-streamer)
* [Recent changes](https://github.com/MD-CWI/afivo-streamer/activity)

# Features

* Adaptive mesh refinement
* OpenMP parallelization
* Multigrid field solvers
* Plasma fluid model for streamer simulations
* Photoionization using a Helmholtz or Monte Carlo approach
* Support for tabulated transport data
* Support for user-defined chemical reactions
* Flexible usage through configuration files and command line arguments

More details about the implementation are given in the [publications](documentation/publications.md) describing afivo-streamer.

# Source code structure

Folder | Contents
---|---
src | The source code of afivo-streamer
src/config_fortran | Module for configuration files
src/lookup_table_fortran | Module for lookup tables
src/rng_fortran | Module for random numbers
afivo | The afivo library
documentation | The documentation (mostly Markdown files)
lib_1d | The afivo-streamer library (1D version)
lib_2d | The afivo-streamer library (2D version)
lib_3d | The afivo-streamer library (3D version)
makefiles | Makefiles for building the library and programs
programs | Streamer simulation programs
tools | Python scripts for e.g. converting data
transport_data | Input data files for simulations

# Opiniated comparison against typical alternatives (advantages and disadvantages)

## Strong points

* **Computational efficiency**: relative to typical other codes, afivo-streamer has a high computational efficiency, see e.g. [this paper](https://doi.org/10.1088/1361-6595/aad768). This usually makes it possible to do 2D simulations on a normal computer rather interactively.
* **Possibility of doing 3D simulations**: due to the computational efficiency it is possible to do relevant 3D simulations in a reasonable time on typical hardware
* **Usage**: the code has frequently been used for streamer simulations in various publications, which also means that things like plasma chemisty and different photoionization models are included
* **Open source**: anyone can in principle change the code (although doing so might not always be easy)

## Neutral points

* **Linux/unix**: basic familiarity with such systems (and the command line) is required to compile and run the code
* **Fortran**: experience with Fortran is required in order to change the code

## Weak points

* **Complexity**: it can be rather difficult to change the code or to add new physics, in particular due to the use of adaptive mesh refinement and parallelization
* **Research type code**: which means there is not that much documentation and there is no graphical user interface
* **Complex geometries**: although curved electrodes can be included, curved dielectrics and complex geometries are generally not supported

# Links and related software

* [Lxcat website](https://lxcat.net), where input data can be obtained
* [Bolsig+](http://www.bolsig.laplace.univ-tlse.fr/), which can be used to compute transport data from electron-neutral cross sections
* [Particle_swarm](https://github.com/MD-CWI/particle_swarm), which can also be used to compute transport data using a Monte Carlo approach
* [Chombo discharge](https://github.com/chombo-discharge/chombo-discharge), “A multiphysics code which uses Chombo for discharge simulations with adaptive mesh refinement (AMR) on embedded boundary grids”, developed by Robert Marskar
* [PASSKey](http://www.plasma-tech.net/parser/passkey/)"Parallel Streamer Solver with Kinetics" by Yifei Zhu and others at LPP
