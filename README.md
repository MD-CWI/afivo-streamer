# Afivo

Afivo (which stands for *Adaptive Finite Volume Octree*) is a framework for
simulations on adaptively refined quadtree and octree grids. It was designed
with a focus on simplicity. Some of the key features are:

* Adaptively refined quadtree and octree grids
* OpenMP parallelization
* FAS multigrid solver (v-cycle and FMG)
* Flexible handling of refinement boundaries and physical boundaries
* Written in modern Fortran
* Silo and VTK unstructured output, which can be visualized with e.g.
  [Visit](https://wci.llnl.gov/simulation/computer-codes/visit)

Documentation is available at
[this link](http://cwimd.nl/other_files/afivo_doc/html/index.html).
[This paper](http://arxiv.org/abs/1701.04329) gives an overview of Afivo, and
[this paper](http://dx.doi.org/10.1088/0963-0252/25/4/044001) contains
simulations performed with Afivo.

## Quick getting started

The instructions below are for the not-so patient, and assume that a recent
version of `gfortran` and `gcc` are installed. See
the [manual](http://cwimd.nl/other_files/afivo_doc/html/index.html) for more
details.

1. Type `make` in the afivo folder. This will first download and compile Silo
   and then compile the Afivo framework itself.
2. Go to the `examples` directory, which will now contain several executables.
3. Run some of the tests by executing e.g. `./computational_domain_2d`,
   `./poisson_basic_2d` etc.
4. The output files can be visualized
   using [Visit](https://wci.llnl.gov/simulation/computer-codes/visit)




