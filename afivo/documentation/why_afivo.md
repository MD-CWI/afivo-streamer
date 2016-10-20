# Why Afivo?

Many physical systems exhibit *multiscale* features, which appear at different
spatial scales. Numerical simulations of such systems can be speed up with
adaptive mesh refinement (AMR), especially if a high-resolution mesh is only
required in a small fraction of the total volume. Afivo is a framework for
simulations on adaptively refined quadtree (2D) and octree (3D) grids. These
grids have a simple structure and connectivity, which allows for efficient
numerical computations. At the same time, they allow for strong local
refinements, as illustrated below.

![Projection of an octree mesh around a streamer discharge. Note that a fine grid is only required in a small region.](branch_view.png)

One of the main reasons to develop Afivo was to be able to experiment
with [geometric multigrid](documentation/multigrid.md) algorithms. Geometric
multigrid can be used to efficiently solve *elliptic partial differential
equations* such as Poisson's equation. There already exist several frameworks
with multigrid solvers. The main reason to develop yet another one was that a
relatively *simple* framework seemed to be missing. Simulations with adaptive
mesh refinement often require significant experimentation, for example to:

* determine a suitable refinement criterion
* compare multigrid algorithms
* investigate different discretizations near refinement boundaries

Afivo was designed to facilitate such experiments, by keeping the implementation
relative simple, as discussed in @ref design:

* Only shared-memory parallelism is supported, so no parallel communication or load balancing is required. (Most frameworks rely on MPI)
* Quadtree and octree grids are used, which are probably the simplest grids that
support adaptive refinement.
* Only cell-centered and face-centered variables are supported.
* Afivo is application-independent, i.e., it includes no code or
algorithms for specific applications.

Because of these simplifications we expect that Afivo can more easily be
modified, thus providing an option in between some of the other `advanced'
frameworks and uniform grid computations.

