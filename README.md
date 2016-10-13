# Afivo

Afivo (<b>A</b>daptive <b>Fi</b>nite <b>V</b>olume <b>O</b>ctree) is a framework
for simulations on adaptively refined quadtree and octree grids. It was designed
with a focus on simplicity.

## Features

* Adaptively refined quadtree and octree grids
* OpenMP parallelization
* FAS multigrid solver (v-cycle and FMG)
* Flexible handling of refinement boundaries and physical boundaries
* Written in modern Fortran
* Silo and VTK unstructured output, which can be visualized with e.g.
  [Visit](https://wci.llnl.gov/simulation/computer-codes/visit).

## Design choices

Some of the design choices are liste below. A more detailed discussion can be
found on [this page](documentation/design.md).

* The refinement ratio is always 2.
* Quantities are either cell-centered or face-centered.
* Parallellization using OpenMP only, no MPI.
* There is always one layer of ghost cells (but of course you can get
  more data from neighbors).
* No 'corner' ghost cells

## Installation & compilation

Instructions can be found on [this page](documentation/compiling.md).

## Authors

* Jannis Teunissen (started the project to do discharge simulations)
* Margreet Nool (worked on documentation, examples, testing)

## Some alternatives

* [Boxlib](https://ccse.lbl.gov/BoxLib/index.html)
* [Paramesh](https://sourceforge.net/projects/paramesh/)
* [Chombo](https://commons.lbl.gov/display/chombo/Chombo+-+Software+for+Adaptive+Solutions+of+Partial+Differential+Equations)
