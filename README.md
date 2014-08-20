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
* Corner ghost cells are **not** used.
* Quantities are either cell-centered or face-centered.
* We use one type of conservative restriction.
* We use 2-1-1 prolongation (because there are no corner ghost cells).
* We support only periodic and dirichlet boundary conditions at first. Neumann
  may also be included if it isn't too much work with the MG solver.
* Parallellization is provided only for shared memory systems (using OpenMP).

### Todo

* Select output format / library from
  [this list](http://www.visitusers.org/index.php?title=Detailed_list_of_file_formats_VisIt_supports).
  This will probably be Silo...
* Create basic datatypes (now wip)
* Decide: dimension independent code or not (e.g., each block gets a variable
n_dim indicating its dimension). Pro: same code for 2D/3D, Con: overhead, many
checks/ifs
* Decide: use tree structure, or index into an array, or index into a "pool" of
  arrays?
* Decide on programming language, started in F90 but that can change
* For morton order / "bit" stuff, C would be nice, especially since it has
  unsigned integers
