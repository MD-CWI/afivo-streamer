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

* A new version 0.3.0 has been pushed, which contains the following changes:

- Improved flux computation: the code now computes flux along a line, which is automatically generalized to 2D/3D simulations
- Inclusion of fluid model based on an energy equation (so m_fluid_lfa has been renamed to m_fluid)
- Cleaning up of unused functionality
- Improvement of boundary conditions for gas dynamics
- Improved time step computation
- The code now avoids tiny negative rates when using cubic spline interpolation of input data
- And more ...

* The code has been moved to Github due to [these changes at gitlab](https://about.gitlab.com/blog/2022/03/24/efficient-free-tier/)


