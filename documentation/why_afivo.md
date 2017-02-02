# Why Afivo?

Many physical systems exhibit *multiscale* features, which appear at different
spatial and temporal scales. Numerical simulations of such systems can be speed
up with adaptive mesh refinement (AMR), especially if a high-resolution mesh is
only required in a small fraction of the total volume. Afivo is a framework for
simulations on adaptively refined [quadtree and octree](@ref quadtree_octree.md)
grids. Because Afivo has no built-in support for specific physics applications a
user has to write his/her own numerical methods. Some key
features/characteristics of the framework are:

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




