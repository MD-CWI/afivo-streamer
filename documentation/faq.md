# Frequently Asked Questions

[TOC]

# What does Afivo stand for? {#afivo-name}

"Adaptive Finite Volume Octree". Note that these names do not
describe the full functionality of the framework, since you can also use e.g.,
quadtrees or finite difference methods.

# Why no MPI? {#why-no-mpi}

There are a couple of reasons for this:

* Simplicity
* There are already frameworks aimed at 'big' simulations.
* It is quite challenging to write a distributed memory code with an efficient
  multigrid solver that scales well to 100 cores or more, in particular when
  there is a lot of grid refinement.
