# Frequently Asked Questions

[TOC]

# What does Afivo stand for?

Adaptive Finite Volume Octree. Note that these names do not describe the full
functionality of the framework, since you can also use e.g., quadtrees or finite
difference methods.

# Why no MPI?

There are a couple of reasons for this:

* There are already many frameworks out there aimed at "big" simulations.
* Many people do not need 'big' simulations, running on more than say 10 cores.
* A lot of the complexity of the current frameworks is in the communication,
  this is much simpler for Afivo. There is much less code, and it is probably
  easier to make changes in a project if you can read all the data from each
  core, so that you do not have to think about MPI. (Although getting good
  OpenMP performance can be quite tricky).
* It quite challenging to write a distributed memory code with an efficient
  multigrid solver that scales well to 100 cores or more, in particular when
  there is a lot of grid refinement.
* If your simulation fits in memory, you can also consider running 5 different
  test cases on a big system instead of one bigger one. The bigger one will
  almost always be less efficient.

# Why use one fixed ghost cell?

We have considered a couple of options, which are listed below with some
remarks:

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

	* Perhaps the most general / elegant idea: do not waste memory on ghost cells
      but just look the values up at the neighbors.
	* Hard to write efficient code: typically you would work on an enlarged copy
      of the current box that includes neighbor data. Copying data takes time,
      and it is hard to write elegant routines for this. For example, to get a
      corner ghost cell, you typically want to use the "side" ghost cell of a
      neighbor, but if these are not stored, they have to be recomputed each
      time.
	* If you do not work on an enlarged copy of the box, indexing is really
      annoying.


