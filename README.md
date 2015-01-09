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
* Corner ghost cells are not automatically filled.
* Quantities are either cell-centered or face-centered.
* We use one type of conservative restriction.
* We implement linear prolongation and injection.
* Boundary conditions are provided by the user.
* Parallellization is provided only for shared memory systems (using OpenMP).

### Todo
* Add more streamer example for 3D code
* Look for "pretty" examples, generate some animations, and put them online.

### Features
* I'll write about them soon

### Design considerations

#### Why don't you fill corner ghost cells by default?

Filling these is relatively easy in 2D, but quite a pain in 3D. There you have to
consider 8 corner points and 12 edges between these corners. These edges can be
shared by multiple neighbors, and filling them in a consistent way is quite
difficult.

#### Why use Fortran (2003+)?

Because it is one of the more convenient languages for scientific computing.

#### Why don't you use MPI?

There are a couple of reasons for this:

* There are already many frameworks out there aimed at "big" simulations.
* I think most people dont need "big" simulations, running on more than say 10
  cores.
* Most of the complexity of the current frameworks is in the communication, this
  is much simpler for AFiVO. There is much less code, and it is probably easier
  to make changes in a project if you can read all the data from each core, so
  that you don't have to think about MPI. (Although getting good OpenMP
  performance is also quite tricky, I admit).
* When your code is more efficient, you can use a smaller system to do the same
  type of simulations. This is what I hope to achieve.
* Most parallel codes don't scale so well, especially if there is a lot of grid
  refinement. The work is then harder to distribute, and more communication is
  required.
* If your simulation fits in memory, then you can also consider running 5
  different test cases on a big system instead of one bigger one. The bigger one
  will almost always be less efficient.

#### Why use one fixed ghost cell?

I have considered a couple of options, which are listed below with some remarks:

* Variable number of ghost cells (depending on the variable)

	* Suitable for all ghost cell requirements, flexible.
	* Having more than 1 ghost cells is not very memory efficient. For example,
	  for a 2D 8x8 block with 2 ghost cells per side, you would have to store
	  12x12 cells (144). So the flexibility of having more than one is not really
	  all that useful.
	* It is harder to write code for a variable number of ghost cells, for
      example, when copying data, should we copy the ghost cells? And writing
      code for filling more than one ghost cell is also hard.
	* Storing variables is annoying, because they cannot be stored in the same
      array (since they have different shapes). This makes indexing harder.

* One ghost cell

	* Restricted, not flexible.
	* Simple to implement because all variables are the same.
	* One ghost cells does not cost too much memory.
	* One ghost cells is convenient for 2nd order schemes, such as the multigrid
      examples, which are quite common.
	* When you need more than one ghost cell, you can simply access the data on
      a neighbor directly. See for example the drift-diffusion test.
	* Perhaps I will add variables without ghost cells in the future, to not
      waste memory for them.

* No ghost cells

	* Perhaps the most general / elegant idea: don't waste memory on ghost cells
      but just look the values up at the neighbors.
	* Hard to write efficient code: typically you would work on an enlarged copy
      of the current box that includes neighbor data. Copying data takes time,
      and it is hard to write elegant routines for this. For example, to get a
      corner ghost cell, you typically want to use the "side" ghost cell of a
      neighbor, but if these are not stored, they have to be recomputed each
      time.
	* If you don't work on an enlarged copy of the box, indexing is really
      annoying.
