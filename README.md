# Afivo's documentation

# Introduction

Afivo (which stands for *Adaptive Finite Volume Octree*) is a framework for
simulations on adaptively refined quadtree and octree grids. It was designed
with a focus on simplicity.

Some of the key features:

* Adaptively refined quadtree and octree grids
* OpenMP parallelization
* FAS multigrid solver (v-cycle and FMG)
* Flexible handling of refinement boundaries and physical boundaries
* Written in modern Fortran
* Silo and VTK unstructured output, which can be visualized with e.g.
  [Visit](https://wci.llnl.gov/simulation/computer-codes/visit)

In order to keep the code simple, the following choices were made:

* The refinement ratio is always 2
* Quantities are either cell-centered or face-centered
* Parallellization using OpenMP only, no MPI
* There is one layer of ghost cells (but you can get more)
* No corner ghost cells

# Documentation contents

The best way to read Afivo's Doxygen-based documentation is
through [this link](http://cwimd.nl/other_files/afivo_doc/html/index.html).

* [Why Afivo?](documentation/why_afivo.md)
* [Installation instructions](documentation/installation.md)
* [About the documentation](documentation/documentation.md)
* [Discussion of design](documentation/design.md)
* [Geometric multigrid](documentation/multigrid.md)
* [Tutorial multigrid](documentation/tutorial_mg.md)
* [Related projects & alternatives](documentation/other_projects.md)
* [FAQ](documentation/faq.md)
* [Authors](documentation/authors.md)
