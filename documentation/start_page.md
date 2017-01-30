# Afivo

## Introduction

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

## Documentation quick links

* [Installation instructions](documentation/installation.md)
* [About the source code](documentation/source_code.md)
* [Why Afivo?](documentation/why_afivo.md)
* [Data structures](documentation/data_structures.md)
* [Multigrid tutorial](documentation/multigrid_tutorial.md)
* [Related projects](documentation/other_projects.md)
* [FAQ](documentation/faq.md)
* [Authors](documentation/authors.md)

## References / how to cite

* \cite Nijdam_Teunissen_2016 Paper in which Afivo was first used
* \cite afivo_arXiv Paper describing the design and functionality of Afivo
