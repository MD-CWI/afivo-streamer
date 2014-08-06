afivo
=====

Adaptive Finite Volume Octree

The goal of this project is to provide a framework for adaptive finite volume
simulations, with a focus on **simplicity** and **adaptive grid refinement**.

The framework should:

* Be easy to use for many types of simulations.
* Efficiently make use of a multicore machine.
* Provide a simple multigrid solver (with point relaxation).
* Provide routines for generating output that can be visualized with

### Design choices

To keep the code lean & fast, we restrict the framework in the following way:

* We use an octree-like grid. The basic blocks in this grid have size
  2<sup>kD</sup>, where k is an integer > 1 and D is the dimension used. This is
  similar as
  [Paramesh](http://www.physics.drexel.edu/~olson/paramesh-doc/Users_manual/amr.html).
* The refinement ratio is always 2.
* There is always a layer of widht 1 of ghost cells.
* Corner ghost cells are **not** used.
* Quantities are either cell-centered or face-centered.
* We use one type of conservative restriction.
* We use 2-1-1 prolongation (because there are no corner ghost cells).
* We support only periodic and dirichlet boundary conditions at first. Neumann
  may also be included if it isn't too much work with the MG solver.
* Parallellization is provided only for shared memory systems (using OpenMP).

### Current status

Nothing is done yet.
