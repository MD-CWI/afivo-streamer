# Afivo's documentation

# Introduction

Afivo (which stands for *Adaptive Finite Volume Octree*) is a framework for
simulations on adaptively refined quadtree and octree grids. It was designed
with a focus on simplicity. Some of the key features are:

* Adaptively refined quadtree and octree grids
* OpenMP parallelization
* FAS multigrid solver (v-cycle and FMG)
* Flexible handling of refinement boundaries and physical boundaries
* Written in modern Fortran
* Silo and VTK unstructured output, which can be visualized with e.g.
  [Visit](https://wci.llnl.gov/simulation/computer-codes/visit)

# Documentation

Afivo's Doxygen-based documentation is
available [here](http://cwimd.nl/other_files/afivo_doc/html/index.html). With a
recent version of Doxygen, it can also be generated locally using `make doc`.

## References / how to cite

A draft of a paper about the Afivo framework is available upon request.

The first paper in which the code was used is:

[The role of free electrons in the guiding of positive streamers](http://dx.doi.org/10.1088/0963-0252/25/4/044001),
S. Nijdam, J. Teunissen, E. Takahashi, U. Ebert, Plasma Sources Sci. Technol.
(2016)

Afivo was first presented in chapter 10 of:

[3D Simulations and Analysis of Pulsed Discharges](http://repository.tue.nl/801516),
J. Teunissen, PhD Thesis, Eindhoven University of Technology, 2015.
