# Geometric multigrid

# Introduction

[Multigrid methods](https://en.wikipedia.org/wiki/Multigrid_method) can be used to efficiently solve elliptic partial differential
equations, such as Poisson's equation. The error in the solution is iteratively
damped on a hierarchy of grids, with the coarse grids reducing the low frequency
(i.e., long wavelength) error components, and the fine grids the high frequency
components.

For an introduction to multigrid methods, consult for example \cite Brandt_2011
or \cite Trottenberg_2000_multigrid. The multigrid implementation in Afivo is
described in \cite afivo_paper. Here is a brief summary:

* Afivo uses FAS (Full Approximation Scheme) multigrid. This means that the (approximate) solution is available on all refinement levels
* A basic V-cycle and FMG-cycle are implemented
* A procedure for the conservative filling of ghost cells is implemented
* There is built-in support for Poisson's/Laplace's equation in 2D, 3D and cylindrical coordinates
* Gauss-Seidel red-black smoothers are implemented
* The user can customize the elliptic operator as well as the prolongation,
  restriction, and correction method.

# How to use multigrid in Afivo

TODO

# Solving different equations

TODO
