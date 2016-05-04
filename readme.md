Afivo
=====

Afivo stands for <b>A</b>daptive <b>Fi</b>nite <b>V</b>olume <b>O</b>ctree. It
is a framework for adaptive finite volume simulations on quadtree / octree
grids. Some alternatives are listed below, but what makes Afivo different is the
focus on **simplicity**. It intentionally supports less features than most
alternatives, which keeps the code relatively simple.

### Features

* All the basic grid functionality for finite volume simulations with
  adaptive mesh refinement.
* Relatively simple.
* OpenMP implementation.
* FAS multigrid solver (v-cycle and FMG) with point relaxation.
* Silo and VTK unstructured output, which can be visualized with e.g.
  [Visit](https://wci.llnl.gov/simulation/computer-codes/visit).

### Design choices

* Afivo uses an octree-like grid. The basic blocks in this grid contain NxN (2D) or
  NxNxN (3D) cells, where N is an even number. This is similar to
  [Paramesh](http://www.physics.drexel.edu/~olson/paramesh-doc/Users_manual/amr.html).
* Parallellization using OpenMP only, no MPI.
* The refinement ratio is always 2.
* There is always a layer of width 1 of ghost cells (you can of course get
  more data from neighbors, this is easy with OpenMP).
* Corner ghost cells are not automatically filled.
* Quantities are either cell-centered or face-centered.
* Boundary conditions are provided by the user.

For some elaboration on the design choices, see
[design considerations](documentation/design.md).

### Installation & compilation

Have a look here: [Compiling the code](documentation/compiling.md).

### Author

Jannis Teunissen, jannis@teunissen.net

### Todo

* Look for "pretty" examples, generate some animations, and put them online.
* Fix Doxygen output

### Some alternatives

* Boxlib
* Paramesh
* Chombo
