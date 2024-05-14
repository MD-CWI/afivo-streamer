# Afivo-streamer

Afivo-streamer is a code for fluid simulations of streamer discharges. It is based on
the [afivo](https://github.com/MD-CWI/afivo) framework, which provides adaptive mesh refinement (AMR) and a multigrid solver for Poisson's equation.

A brief summary of features:

* 1D, 2D, 3D and a cylindrical fluid model
* Electrodes
* Chemistry
* Gas dynamics
* OpenMP parallelization

## Documentation

Documentation is available at http://teunissen.net/afivo_streamer

## Requirements

* A unix-like system (e.g., GNU/Linux)
* Gfortran 5 or newer

## News

* The Hypre library has been updated. It has to be recompiled using the `build_hypre.sh` script in the `afivo/external_libraries` folder.
