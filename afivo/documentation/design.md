# Design of Afivo

# Introduction {#design-intro}

Given the fact that there are already several frameworks available, why develop
another one? The main reason was that a relatively simple framework seemed to be
missing. Therefore, much of the design choices were made to keep the framework
simple, e.g.:

* The refinement ratio is always 2
* Quantities are either cell-centered or face-centered
* Parallellization using OpenMP only, no MPI
* There is one layer of ghost cells (but you can get more)
* No corner ghost cells

![](quadtree_cex4.png)

# One ghost cell {#design-one-ghost-cell}

In Afivo, we have opted for a simple approach: there is always one layer of
ghost cells for cell-centered variables. When additional ghost cells are
required, these can of course still be computed, there is just no default
storage for them.

# No corner ghost cells {#design-no-corner-ghost}

In Afivo, corner ghost cells are not used. The reason for this is that in three
dimensions, the situation is quite complicated: there are eight corners, twelve
edges and six sides. It is hard to write an elegant routine to fill all these
ghost cells, especially because the corners and edges have multiple neighbors.
Therefore, only the sides of a box are considered in Afivo. This means that
Afivo is not suitable for stencils with diagonal terms.
