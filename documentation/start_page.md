# Introduction to Afivo

[TOC]

![Snapshot of a simulation performed with Afivo](branch_view.png)

# What is Afivo? {#what-is-afivo}

Afivo is a generic framework for simulations on adaptively refined quadtree and octree grids \cite afivo_paper. 
Some key features/characteristics of the framework are:

* Adaptively refined [quadtree and octree grids](@ref documentation/quadtree_octree.md)
* OpenMP [parallelization](@ref documentation/parallelization.md)
* Geometric [multigrid routines](@ref documentation/multigrid.md)
* Flexible handling of [refinement](@ref documentation/grid_refinement.md) and [physical boundaries](@ref documentation/boundary_conditions.md)
* Silo and VTK unstructured [output](@ref documentation/writing_viewing_output.md)
* Written in modern Fortran
* Open source

# Installation

See the [installation instructions](@ref documentation/installation.md), afterwards have a look at the <a href="examples.html">examples</a>.

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

See the <a href="namespaces.html">Modules</a> page for all included modules.
Have a look at the @ref doc-contents to get started with Afivo.

# Relevant publications {#publications}

* \cite afivo_paper Paper describing Afivo
* \cite afivo_streamer_paper Paper describing a streamer simulation code based on Afivo
* \cite Nijdam_Teunissen_2016 Paper in which Afivo was first used
* Papers in which Afivo was used: \cite Bagheri_2019 \cite Bagheri_2018 \cite Malag_n_Romero_2020 \cite Li_2020 \cite Li_2020b
