afivo
=====

<b>A</b>daptive <b>Fi</b>nite <b>V</b>olume <b>O</b>ctree

The goal of this project is to provide a framework for adaptive finite volume
simulations, with a focus on **simplicity** and **adaptive grid refinement**.

The framework should:

* Provide all the basic grid functionality for finite volume simulations with
  adaptive mesh refinement.
* Be relatively easy to use.
* Efficiently make use of a multicore machine.
* Provide a simple multigrid solver (with point relaxation).
* Provide routines for generating output that can be visualized with
  [Visit](https://wci.llnl.gov/simulation/computer-codes/visit).

### Design choices

To keep the code lean & fast, we restrict the framework in the following way:

* We use an octree-like grid. The basic blocks in this grid have size
  2<sup>kD</sup>, where k is an integer > 1 and D is the dimension used. This is
  similar to
  [Paramesh](http://www.physics.drexel.edu/~olson/paramesh-doc/Users_manual/amr.html).
* The refinement ratio is always 2.
* There is always a layer of width 1 of ghost cells.
* Corner ghost cells are used (although optionally).
* Quantities are either cell-centered or face-centered.
* We use one type of conservative restriction.
* We use bi/tri-linear prolongation.
* Boundary conditions are provided by the user.
* Parallellization is provided only for shared memory systems (using OpenMP).

### Todo
* Use morton order for enhancing data locality
* (Idem) Reorder memory
* Fill ghost cells from boundary condition

### Status
* Using vtk unstructured for output
* Use 1d array with indexing as data storage
