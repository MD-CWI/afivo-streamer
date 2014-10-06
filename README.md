afivo
=====

<b>A</b>daptive <b>Fi</b>nite <b>V</b>olume <b>O</b>ctree

The goal of this project is to provide a framework for adaptive finite volume
simulations, with a focus on **simplicity** and **adaptive grid refinement**.

The framework should:

* Provide all the basic grid functionality for finite volume simulations with
  adaptive mesh refinement.
* Be relatively easy to use.
* Be easy to modify.
* Efficiently make use of a multicore machine.
* Provide a FAS multigrid solver (v-cycle and fmg) with point relaxation.
* Provide routines for generating output that can be visualized with
  [Visit](https://wci.llnl.gov/simulation/computer-codes/visit).

### Design choices

To keep the code lean & fast, we restrict the framework in the following way:

* We use an octree-like grid. The basic blocks in this grid contain NxN (2D) or
  NxNxN (3D) cells, where N is an even number. This is similar to
  [Paramesh](http://www.physics.drexel.edu/~olson/paramesh-doc/Users_manual/amr.html).
* The refinement ratio is always 2.
* There is always a layer of width 1 of ghost cells (however you can always get
  more data from neighbors)
* Corner ghost cells can be optionally filled.
* Quantities are either cell-centered or face-centered.
* We use one type of conservative restriction.
* We implement bi/tri-linear prolongation and injection.
* Boundary conditions are provided by the user.
* Parallellization is provided only for shared memory systems (using OpenMP).

### Todo
* Reorganize 2D code and start writing 3D code

### Features
* I'll write about them soon
