# Overview of data structures {#main-data-structures}

We now start with the description of the implementation of Afivo.
First, the properties of *orthtree* meshes are discussed.
Three data types are used to store these meshes: boxes, levels and trees. These
data types are described below.

|       |        | 
| ----: |  :---- |
<img src="../../documentation/figures/quadtree_cex1.png" width=100px /> | <img src="../../documentation/figures/quadtree_cex2.png" width=100px /> 
<img src="../../documentation/figures/quadtree_cex3.png" width=100px /> | <img src="../../documentation/figures/quadtree_cex4.png" width=100px /> |
<a name="fig_example-quadtree" />
<img src="../../documentation/figures/box_indices.png" width=200px style="right" /> 

**Figure 4**. **TODO**: refine in one corner Upper: Example of a quadtree mesh
that gets refined. Here boxes contain \f$2 \times 2\f$ cells, and different
boxes have different colors.
Below: The spatial indices of the boxes. When a box with indices \f$(i,j)\f$ is refined,
its children have indices \f$(2i-1,2j-1)\f$ up to \f$(2i,2j)\f$.}

## Orthtree meshes {#main-quad-octree}

In Afivo, quadtree (2D) and octree (3D) meshes are used, which can be described
by the following rules:

	- The mesh is constructed from boxes that each contain \f$N^D\f$ cells, where
	\f$D\f$ is the number of coordinates, or the dimension of the problem,
	 and \f$N\f$ is an even number.

	- When a box with a grid spacing \f$h\f$ is refined, it is subdivided in \f$2^D\f$
	**children** with grid spacing \f$h/2\f$.

	- The difference in refinement level for adjacent boxes is at most one.
	This is called **2:1 balance**, or, proper nesting.

<a href="#fig_example-quadtree">Figure 4</a> shows an example of a quadtree that gets
refined.
All the boxes are stored in a single one-dimensional array, so that an integer
index can be used to point to a box, see section \ref main-tree-type below.


## Box data type {#main-box-type}

A *box* is the basic mesh unit in Afivo. Each box consists of \f$N^D\f$ grid
cells, where \f$N\f$ has to be an even number and \f$D\f$ is the spatial dimension. In
<a href="#fig_example-quadtree">Figure 4</a>, there is for example a box with \f$2 \times 2\f$
cells at coordinate \f$(1,1)\f$. Each box stores its parent, an array of \f$2^D\f$
children and an array of \f$2D\f$ neighbors. In <a href="#fig_location-box-indices">Figure 5a</a>,
 the indices of the children and the neighbors are shown. 
A special value
<a class="el" href="namespacem__af__t.html#aa0395be4a37c035fae9b5663b2d41698">af_no_box</a>
(which is zero) is
used to indicate that a parent, child or neighbor does not exist. In the case of
neighbors, boundary conditions are specified by negative numbers.

Two types of cell data are supported by default: cell-centered data and
face-centered data, see <a href="#fig_location-box-indices">Figure 5b</a>.
These are stored in \f$D+1\f$-dimensional arrays, so that multiple variables can be
stored per location.
Furthermore, boxes contain some `convenience` information, such at their
refinement level, minimum coordinate and spatial index.

<a name="fig_location-box-indices" />
a) | b) 
------------- | -------------
<img src="../../documentation/figures/children_neighbors.png" width=200px /> | <img src="../../documentation/figures/location_cc_fx.png" width=250px />


**Figure 5**. a) Each box contains an array of children and neighbors.
The ordering of these arrays is shown here for a 2D box.
Red: children, blue: neighbors.
b) Location and indices of the cell-centered variables (black dots) and the
face-centered variables in the x-direction (red dots) for a box of
\f$2\times 2\f$ cells.}

## Level data type {#main-level-type}

The *level* data type contains three lists:

	- A list with all the boxes at refinement level \f$l\f$
	- A list with the *parents* (boxes that are refined) at level \f$l\f$
	- A list with the *leaves* (boxes that are not refined) at level \f$l\f$

This separation is often convenient, because some algorithms operate only on
leaves while others operate on parents or on all boxes.
These lists contain the integer indices of the boxes in the tree data structure
described below.


## Tree data type {#main-tree-type}

The tree data type contains all the data of the mesh.
Most importantly, it stores two arrays: one that contains all the boxes and one
that contains all the levels.

\verbatim
Since Afivo is implemented in Fortran, these arrays start at index one.
\endverbatim

Some other information is also stored: the current maximum refinement level, the
number of cells per box-dimension \f$N\f$, the number of face and cell-centered
variables and the grid spacing \f$\Delta x\f$ at the coarsest level.


# Methods {#main-methods}

In this section we give a brief overview of the most important methods in Afivo.
The names of the methods for two-dimensional meshes are used, which have the
prefix `a2_`. The three-dimensional analogs have, not surprisingly, a prefix `a3_`.

## Creating the initial mesh {#main-init-mesh}

In Afivo the coarsest mesh, which covers the full computational domain, is not
supposed to change. To create this mesh there is a routine
<a class="el" href="namespacem__a2__core.html#ab7007734c1625a6057a63f4676ce7244">a2_set_base</a>,
which takes as input the spatial indices of
the coarse boxes and their neighbors. In <a href="#fig_set-base">Figure 6</a>, a 2D
example is shown for creating a single coarse box at index \f$(1,1)\f$. This box is
its own neighbor in all four directions, or in other words, there are periodic
boundary conditions. Physical (non-periodic) boundaries can be indicated by a
negative index for the neighbor. By adjusting the neighbors one can specify
different geometries, the possibilities include meshes that contain a hole, or
meshes that consist of two isolated parts. The treatment of boundary conditions
is discussed in section \ref main-fill-ghost-cell.

<a name="fig_set-base" />
\verbatim
! Initialize tree
call a2_init(tree, & ! Tree to initialize
     box_size, &     ! Number of cells per coordinate in a box
     n_var_cell, &   ! Number of face-centered variables
     n_var_face, &   ! Number of cell-centered variables
     dr)             ! Distance between cells on base level

! Set the spatial index and neighbors
ix_list(1:2, 1) = 1  ! One box at (1,1)
nb_list(1:4, 1) = 1  ! Periodic box is its own neighbor

! Create the base mesh
call a2_set_base(tree, ix_list, nb_list)
\endverbatim
**Figure 6**. Fortran code fragment that shows how a base mesh can be constructed.
In this case, there is one box at \f$(1,1)\f$, with periodic boundary conditions.


## The refinement procedure {#main-ref-procedure}

Mesh refinement can be performed by calling the
<a class="el" href="namespacem__a2__core.html#a3f0dd0f48f710e79b6b2085177b19dd6">a2_adjust_refinement</a>
routine, which should be passed a user-defined refinement function.
This function is then called for each box, and should set the refinement flag of
the box to one of three values: refine (add children), derefine (remove this
box) or keep refinement.
It is also possible to set the refinement flags for other boxes than the current
one, which can for example be useful to extend refinements to a neighbor.

In Afivo, a box is either fully refined (with \f$2^D\f$ children) or not refined.
Furthermore, 2:1 balance is ensured, so that there is never a jump of more than
one refinement level between neighboring boxes.
These constraints are automatically handled, so that the user-defined refinement
function does not need to impose them.

Each call to
<a class="el" href="namespacem__a2__core.html#a3f0dd0f48f710e79b6b2085177b19dd6">a2_adjust_refinement</a>
changes the mesh by at most one level.
To introduce larger changes one should call the routine multiple times.
A number of rules is used to make the user-supplied refinement consistent:

	- Only leaves can be removed (because the grid changes by at most one
	level at a time)
	- A box flagged for refinement will always be refined, including neighbors
	that are required for 2:1 balance
	- Boxes cannot be removed if that would violate 2:1 balance
	- If all the \f$2^D\f$ children of a box are flagged for removal, and the box
	itself not for refinement, then the children are removed
	- Boxes at level one cannot be removed
	- Boxes cannot be refined above the maximum allowed refinement level


The <a class="el" href="namespacem__a2__core.html#a3f0dd0f48f710e79b6b2085177b19dd6">a2_adjust_refinement</a>
routine returns information on the added and removed boxes per level, so that a user can set values on the new
boxes or clean up data on the removed ones.

When boxes are added or removed in the refinement procedure, their connectivity
is automatically updated.
References to a removed box are removed from its parent and neighbors.
When a new box is added, its neighbors are found through its parent.
Three scenarios can occur: the neighbor can be one of the other children of the
parent, the neighbor can be a child from the neighbor of the parent, or the
neighbor does not exist.
In the latter case, there is a refinement boundary, which is indicated by the
special value
<a class="el" href="namespacem__af__t.html#aa0395be4a37c035fae9b5663b2d41698">af_no_box</a>.

## Filling of ghost cells {#main-fill-ghost-cell}

When working with numerical grids that are divided in multiple parts, it is
often convenient to use *ghost cells*.
The usage of ghost cells has two main advantages: algorithms can operate on the
different parts without special care for the boundaries, and algorithms can
straightforwardly operate in parallel.

In Afivo each box has one layer of ghost cells for its cell-centered variables,
as illustrated in <a href="#fig_location-box-indices">Figure 5b</a>.
The built-in routines only fill the ghost cells on the sides of boxes, not those
on the corners.
The reasons for implementing ghost cells in this way are discussed in sections
\ref main-one-ghost-cell and \ref main-no-corner-ghost .

For each side of a box, ghost cells can be filled in three ways
	- If there is a neighboring box, then the ghost cells are simply copied
	from this box
	- If there is a physical boundary, then the box is passed to a user-defined
	routine for boundary conditions
	- If there is a refinement boundary, then the box is passed to another
	user-defined routine

Physical boundaries are indicated by negative values for the neighbor index, and
these values are passed on to the user-defined routine.
In this way, one can set up different types of boundary conditions.


## Interpolation and restriction {#main-interp-restrict}

Because corner ghost cells are not filled by Afivo, interpolation schemes cannot
use diagonal elements. Instead of the standard bilinear and trilinear
interpolation schemes, a \f$2-1-1\f$ and \f$1-1-1-1\f$ scheme is used in 2D and 3D,
respectively. These interpolation schemes use information from the closest and
second-closest neighbors; the 2D case is illustrated in
<a href="#fig_interp-2d">Figure  7</a>. Zeroth-order interpolation is also included, in
which the coarse values are simply copied without any interpolation. As a
restriction method (going from fine to coarse) Afivo just includes averaging, in
which the parent gets the average value of its children.

A user can of course implement higher order interpolation and restriction
methods, by using information from additional grid locations.
It is generally quite complicated to do this consistently near refinement
boundaries.

<a name="fig_interp-2d" />
<img src="../../documentation/figures/interp_2d.png" width=200px />
**Figure 7**. Schematic drawing of \f$2-1-1\f$ interpolation.
The three nearest coarse grid values are used to interpolate to the center
of a fine grid cell.
Note that the same interpolation scheme can be used for all fine grid cells,
because of the symmetry in a cell-centered discretization.


## The list of boxes {#main-tidy-up}

In Afivo, all the boxes are stored in a single array. New boxes are always added
to the end of the array, which means that removed boxes leave a `hole`. There is
a route
<a class="el" href="namespacem__a2__core.html#a84a6486bbb34876a27c9a7955e8a17ff">a2_tidy_up</a>
to tidy up the array: all the unused
boxes are moved to the end of the array, and the boxes that are still in use are
sorted by their refinement level. Furthermore, for each level, the boxes are
ranked according to their Morton index Xcite{Morton_1966}.

Sometimes, extra storage is required when
<a class="el" href="namespacem__a2__core.html#a3f0dd0f48f710e79b6b2085177b19dd6">a2_adjust_refinement</a>
has to add new boxes to the mesh.
In such a case, the array of boxes is simply resized so that there is enough
space.

## Producing output {#main-output}

It is important that one is able to quickly and conveniently visualize the
results of a simulation. Afivo supports two output formats: VTK unstructured
files and Silo files.

For VTK files, Afivo relies on the unstructured format, which support much more
general grids than quadtree and octree meshes.
This format should probably only be used for grids of moderate sizes (e.g.,
\f$10^5\f$ or \f$10^6\f$ cells), because visualizing larger grids can be computationally
expensive.
Although there is some support for octrees in VTK, this support does not yet
extend to data visualization programs such as Paraview Xcite{www_paraview} and
Visit Xcite{HPV:VisIt}.

Afivo also supports writing Silo files.
These files contain a number of Cartesian blocks ( **quadmeshes** in Silo''s terminology)
that can each contain multiple boxes.
This is done by starting with a region \f$R\f$ that contains a single box.
If all the neighbors to the left of \f$R\f$ exist, have no children and are not yet
included in the output, then these neighbors are added to \f$R\f$.
The procedure is repeated in all directions, until \f$R\f$ can no longer grow.
Then \f$R\f$ represents a rectangular collection of boxes which can be added to the
output, and the procedure start again from a new box that is not yet included.
This merging of boxes is done because writing and reading a large number
separate meshes can be quite costly with the Silo library.
