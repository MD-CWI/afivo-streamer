# Introduction to Afivo

[TOC]

# What is Afivo? {#what-is-afivo}

Afivo is a framework for simulations on adaptively
refined [quadtree and octree](@ref quadtree_octree.md) grids. Because Afivo has
no built-in support for specific physics applications a user has to write
his/her own numerical methods. Some key features/characteristics of the
framework are:

* Adaptively refined quadtree and octree grids
* OpenMP parallelization
* Geometric multigrid routines
* Flexible handling of refinement and physical boundaries
* Written in modern Fortran
* Fully open source
* Application-independent
* Silo and VTK unstructured output

The motivation for developing Afivo is discussed in \cite afivo_paper. In
summary, the main reason was to provide a relatively simple framework that can
easily be modified.

# For which problems can Afivo be used? {#which-problems}

Afivo can be used to simulate physical systems exhibiting *multiscale* features,
e.g. features that appear at different spatial and temporal scales. Numerical
simulations of such systems benefit from Afivo's adaptive mesh refinement (AMR),
especially if a high-resolution mesh is only required in a small fraction of the
total volume.

# What is included? {#included-functionality}

Afivo provides general functionality for parallel simulations with adaptive mesh
refinement:

* It can adjust the refinement according to user-supplied information
* It stores cell-centered and face-centered variables
* It provides routines to perform restriction and prolongation (to convert fine
  grid values to coarse ones and vice versa)
* It can fill so-called *ghost cells*, which allow the user to perform
  computations as on a uniform grid
* It can (help) solve elliptic partial differential equations with the built-in
  multigrid methods
* It can write output in Silo and VTK format, which can directly be visualized
  with tools such
  as [Visit](https://wci.llnl.gov/simulation/computer-codes/visit/)

# What is not included? {#excluded-functionality}

Besides the multigrid methods, there are no built-in solvers. To use Afivo for
e.g. hyperbolic problems, a user will have to implement both a spatial and
temporal discretization. Some of the included @ref examples-page show how this
can be done.
