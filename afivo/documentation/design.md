# Design of Afivo

# Introduction {#design-intro}

Given the fact that there are already several frameworks available, why develop
another one? The main reason was that a relatively simple framework seemed to be
missing. Therefore, much of the design choices were made to keep the framework
simple.

# Grids {#design-grid}

Two types of structured meshes are commonly used: *block-structured* meshes and
*orthtree* meshes, see the figure below.

Block-structured meshes are more general than orthtree meshes: any orthtree mesh
is also a block-structured mesh, whereas the opposite is not true.
Some of the advantages and disadvantages of these approaches are:

- In a block-structured mesh, blocks can have a flexible size. Computations on
  larger blocks are typically more efficient, especially when ghost cells are
  required (virtual cells on the boundary of a block), which is typically the
  case.

- For an orthtree grid, there is a trade-off: larger block sizes allow for more
  efficient computations, but reduce the adaptivity of the mesh. For a
  block-structured grid, there is a similar trade-off: in principle it can be
  refined in a more flexible way, but adding many refined blocks increases the
  overhead.

- The connectivity of the mesh is simpler for an orthtree mesh, because each
  block has the same number of cells, and blocks are refined in the same way.
  This also ensures a simple relation between fine and coarse meshes. These
  properties make operations such as prolongation and restriction easier to
  implement, especially in parallel.

Block-structured grid | Quadtree grid
---|---
![a](block_structured.png) | ![](quadtree_cex4.png)

# One ghost cell {#design-one-ghost-cell}

There are essentially two ways to implement ghost cells in a framework such as Afivo.

1. Ghost cells are not stored for boxes. When a computation has to be performed
   on a box, there are typically two options: algorithms can be made aware of
   the mesh structure, or a box can be temporarily copied to an enlarged box on
   which ghost cells are filled.
2. Each box includes ghost cells, either a fixed number for all variables or a
   variable-dependent number.

Storing ghost cells can be quite costly. For example, adding two layers of ghost
cells to a box of \f$8^3\f$ cells requires \f$(12/8)^3 = 3.375\f$ times as much
storage. With one layer, about two times as much storage is required. Not
storing ghost cells prevents this extra memory consumption. However, some
operations can become more complicated to program, for example when some type of
interpolation depends on coarse-grid ghost cells. Furthermore, one has to take
care not to unnecessarily recompute ghost values, and parallelization becomes
slightly harder.

If ghost cells are stored for each box, then there are still two options: store
a fixed number of them for each variable, or let the number of ghost cells vary
per variable. In Afivo, we have opted for the simplest approach: there is always
one layer of ghost cells for cell-centered variables. For numerical operations
that depend on the nearest neighbors, such as computing a second order
Laplacian, one ghost cell is enough. When additional ghost cells are required,
these can of course still be computed, there is just no default storage for
them.

# No corner ghost cells {#design-no-corner-ghost}

In Afivo, corner ghost cells are not used.
The reason for this is that in three dimensions, the situation is quite
complicated: there are eight corners, twelve edges and six sides.
It is hard to write an elegant routine to fill all these ghost cells, especially
because the corners and edges have multiple neighbors.
Therefore, only the sides of a box are considered in Afivo.
This means that Afivo is not suitable for stencils with diagonal terms.

# OpenMP for parallelism {#design-openmp}

The two conventional methods for parallel computing are OpenMP (shared memory)
and MPI (communicating tasks).
Afivo was designed for small scale parallelism, for example using at most 16
cores, and therefore only supports OpenMP.
Compared to an MPI implementation, the main advantage of OpenMP is simplicity:
data can always be accessed, sequential (user) code can easily be included,
there is no need for load balancing and no communication between processes needs
to be set up.

Most operations in Afivo loop over a number of boxes, for example the leaves at
a certain refinement level.
All such loops have been parallelized by adding OpenMP statements around them,
for example as shown below:

    do lvl = 1, tree%max_lvl
       !$omp parallel do private(id)
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call my_method(tree%boxes(id))
       end do
       !$omp end parallel do
    end do

The parallel speedup that one can get depends on the cost of the algorithm that
one is using. The communication cost (updating ghost cells) is always about the
same, so that an expensive algorithm will show a better speedup. Furthermore, on
a shared memory system, it is not unlikely for an algorithm to be memory-bound
instead of CPU-bound.
