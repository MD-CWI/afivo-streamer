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

* \cite afivo_arXiv Paper describing Afivo
* \cite Nijdam_Teunissen_2016 Paper in which Afivo was first used
