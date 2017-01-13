# Afivo's documentation

# Introduction

Afivo (which stands for *Adaptive Finite Volume Octree*) is a framework for
simulations on adaptively refined quadtree and octree grids. It was designed
with a focus on simplicity.

Some of the key features:

* Adaptively refined quadtree and octree grids
* OpenMP parallelization
* FAS multigrid solver (v-cycle and FMG)
* Flexible handling of refinement boundaries and physical boundaries
* Written in modern Fortran
* Silo and VTK unstructured output, which can be visualized with e.g.
  [Visit](https://wci.llnl.gov/simulation/computer-codes/visit)

In order to keep the code simple, the following choices were made:

* The refinement ratio is always 2
* Quantities are either cell-centered or face-centered
* Parallellization using OpenMP only, no MPI
* There is one layer of ghost cells (but you can get more)
* No corner ghost cells

# Documentation contents

The best way to read Afivo's Doxygen-based documentation is
through [this link](http://cwimd.nl/other_files/afivo_doc/html/index.html).

* [Why Afivo?](http://cwimd.nl/other_files/afivo_doc/html/md_documentation_why_afivo.html)
* [Installation instructions](http://cwimd.nl/cwimd.nl/other_files/afivo_doc/html/md_documentation_installation.html)
* [About the documentation](http://cwimd.nl/cwimd.nl/other_files/afivo_doc/html/md_documentation_installation.html)
* [Discussion of design](http://cwimd.nl/cwimd.nl/other_files/afivo_doc/html/md_documentation_design.html)
* [Geometric multigrid](http://cwimd.nl/cwimd.nl/other_files/afivo_doc/html/md_documentation_multigrid.html)
* [Tutorial multigrid](http://cwimd.nl/cwimd.nl/other_files/afivo_doc/html/md_documentation_tutorial_mg.html)
* [Related projects & alternatives](http://cwimd.nl/cwimd.nl/other_files/afivo_doc/html/md_documentation_other_projects.html)
* [FAQ](http://cwimd.nl/cwimd.nl/other_files/afivo_doc/html/md_documentation_faq.html)
* [Authors](http://cwimd.nl/cwimd.nl/other_files/afivo_doc/html/md_documentation_authors.html)

## References / how to cite

A draft of a paper about the Afivo framework is available upon request.

The first paper in which the code was used is:

[The role of free electrons in the guiding of positive streamers](http://dx.doi.org/10.1088/0963-0252/25/4/044001),
S. Nijdam, J. Teunissen, E. Takahashi, U. Ebert, Plasma Sources Sci. Technol.
(2016)

Afivo was first presented in chapter 10 of:

[3D Simulations and Analysis of Pulsed Discharges](http://repository.tue.nl/801516),
J. Teunissen, PhD Thesis, Eindhoven University of Technology, 2015.
