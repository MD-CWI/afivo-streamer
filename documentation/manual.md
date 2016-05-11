\section Contents
<ul>
	<li> \ref sect_intro </li>
	<li> \ref sect_design-discussion
		<ul>
			<li> \ref sect_grid </li>
			<li> \ref sect_one-ghost-cell </li>
			<li> \ref sect_no-corner-ghost </li>
			<li> \ref sect_openmp </li>
		</ul>
	</li>
	<li> \ref sect_data-structures
		<ul>
			<li> \ref sect_quad-octree </li>
			<li> \ref sect_box-type </li>
			<li> \ref sect_level-type </li>
			<li> \ref sect_tree-type </li>
		</ul>
	</li>
	<li> \ref sect_methods
		<ul>
			<li> \ref sect_init-mesh </li>
			<li> \ref sect_ref-procedure </li>
			<li> \ref sect_fill-ghost-cell </li>
			<li> \ref sect_interp-restrict </li>
			<li> \ref sect_tidy-up </li>
			<li> \ref sect_output </li>
		</ul>
	</li>
	<li> \ref sect_afivo-multigrid
		<ul>
			<li> \ref sect_mg-v-cycle </li>
			<li> \ref sect_mg-fmg-cycle </li>
			<li> \ref sect_gsrb</li>
			<li> \ref sect_mg-ghost-cells
				<ul>
					<li> \ref sect_3d_case </li>
					<li> \ref sect_mg-varepsilon </li>
					<li> \ref sect_mg-cyl </li>
				</ul>
			</li>
			<li> \ref sect_mg-examples </li>
		</ul>
	</li>
	<li>
		\ref sect_example-fluid
		<ul>
			<li> \ref sect_model-formulation </li>
			<li> \ref sect_flux-calc-time-stepping </li>
			<li> \ref sect_refinement-crit </li>
			<li> \ref sect_results </li>
		</ul>
	</li>
	<li> \ref sect_References </li>
</ul>

\section sect_intro Introduction

Many physical systems have *multiscale* properties, i.e., their features
appear at different spatial scales. Numerical simulations of such systems can be
speed up with adaptive mesh refinement (AMR), especially if a high-resolution
mesh is only required in a small fraction of the total volume. Here, we present
Afivo, a framework for simulations on adaptively refined quadtree (2D) and
octree (3D) grids. These grids have a simple structure and connectivity, which
allows for efficient numerical computations. At the same time, they allow for
strong local refinements, as illustrated in <a href="#fig_ex-orthtree" >Figure 1</a>.

One of the main reasons to develope Afivo (which stands for Adaptive Finite
Volume Octree) was to be able to experiment with geometric multigrid algorithms.
This has led to the development of a geometric multigrid solver for quadtrees /
octrees that handles refinement boundaries consistently. The framework with the
multigrid routines is available under the <code>GNU GPLv3</code> license; both were
implemented in Fortran 2008.

<a name="fig_ex-orthtree" />
<img src="../../figures/two_boxes.png" width=300px />
**Figure 1**. In \f$D\f$ dimensions, a quadtree or octree grid consist of boxes that each contain \f$N^D\f$ cells. A box with a grid spacing \f$h\f$ can be subdivided into \f$2^D\f$ smaller boxes with grid spacing \f$h/2\f$. By refining boxes up to different levels, local refinements can be created. In Afivo the **2:1 balance** is ensured, which means that between neighboring cells the refinement level differs by at most one.

Below, the basic data
structures and methods are described in sections \ref sect_data-structures
and \ref sect_methods, after which the multigrid implementation is described in
section \ref sect_afivo-multigrid . Finally, simulation examples are presented in section \ref sect_afivo-examples.


\section sect_design-discussion Design discussion

\subsection sect_grid Grids

Two types of structured meshes are commonly used: *block-structured* meshes and *orthtree*
meshes. Examples of these meshes are shown in <a href="#fig_block-orth-example" >Figure 2</a>.

Note that we here refer to quadtrees and octrees as *orthtrees*, because the
general name for quadrants and octants is orthants Wong \cite www_orthtree_wong.


Block-structured meshes are more general than orthtree meshes: any orthtree mesh
is also a block-structured mesh, whereas the opposite is not true.
Some of the advantages and disadvantages of these approaches are:

	- In a block-structured mesh, blocks can have a flexible size.
	Computations on larger blocks are typically more efficient, especially when
	ghost cells are required (virtual cells on the boundary of a block), which is
	typically the case.

	- For an orthtree grid, there is a trade-off: larger block sizes allow for
	more efficient computations, but reduce the adaptivity of the mesh. For a
	block-structured grid, there is a similar trade-off: in principle it can be
	refined in a more flexible way, but adding many refined blocks increases the
	overhead.

	- The connectivity of the mesh is simpler for an orthtree mesh, because
	each block has the same number of cells, and blocks are refined in the same
	way. This also ensures a simple relation between fine and coarse meshes. These
	properties make operations such as prolongation and restriction easier to
	implement, especially in parallel.


<a name="fig_block-orth-example" />
a) A block-structured grid | b) A quadtree grid
------------- | -------------
<img src="../../figures/block_structured.png" height=200px /> | <img src="../../figures/quadtree_cex4.png" height=200px />

**Figure 2**. a) Example of a block-structured grid, taken from Xcite{www_bell_block_amr}.
b) A quadtree grid consisting of boxes of \f$2 \times 2\f$ cells.

\subsection sect_one-ghost-cell One ghost cell

There are essentially two ways to implement ghost cells in a framework such as Afivo.

	1. Ghost cells are not stored for boxes.
	When a computation has to be performed on a box, there are typically two
	options: algorithms can be made aware of the mesh structure, or a box can be
	temporarily copied to an enlarged box on which ghost cells are filled.
	2. Each box includes ghost cells, either a fixed number for all variables
	or a variable-dependent number.


Storing ghost cells can be quite costly. For example, adding two layers of ghost
cells to a box of \f$8^3\f$ cells requires \f$(12/8)^3 = 3.375\f$ times as much storage.
With one layer, about two times as much storage is required. Not storing ghost
cells prevents this extra memory consumption. However, some operations can
become more complicated to program, for example when some type of interpolation
depends on coarse-grid ghost cells. Furthermore, one has to take care not to
unnecessarily recompute ghost values, and parallelization becomes slightly
harder.

If ghost cells are stored for each box, then there are still two options: store
a fixed number of them for each variable, or let the number of ghost cells vary
per variable. In Afivo, we have opted for the simplest approach: there is always
one layer of ghost cells for cell-centered variables. For numerical operations
that depend on the nearest neighbors, such as computing a second order
Laplacian, one ghost cell is enough. When additional ghost cells are required,
these can of course still be computed, there is just no default storage for
them.

\subsection sect_no-corner-ghost No corner ghost cells

In Afivo, corner ghost cells are not used.
The reason for this is that in three dimensions, the situation is quite
complicated: there are eight corners, twelve edges and six sides.
It is hard to write an elegant routine to fill all these ghost cells, especially
because the corners and edges have multiple neighbors.
Therefore, only the sides of a box are considered in Afivo.
This means that Afivo is not suitable for stencils with diagonal terms.

\subsection sect_openmp OpenMP for parallelism


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
for example as in <a href="#fig_omp-example">Figure 3</a>.

\verbatim
do lvl = 1, tree%max_lvl
   !$omp parallel do private(id)
   do i = 1, size(tree%lvls(lvl)%ids)
      id = tree%lvls(lvl)%ids(i)
      call my_method(tree%boxes(id))
   end do
   !$omp end parallel do
end do
\endverbatim

<a name="fig_omp-example" />
**Figure 3**. Fortran code fragment that shows how to call
my_method for all the boxes in a tree, from level 1 to
the maximum level.
Within each level, the routine is called in parallel using OpenMP.}

The parallel speedup that one can get depends on the cost of the algorithm that
one is using.
The communication cost (updating ghost cells) is always about the same, so that
an expensive algorithm will show a better speedup.
Furthermore, on a shared memory system, it is not unlikely for an algorithm to
be memory-bound instead of CPU-bound.


\section sect_data-structures Overview of data structures

We now start with the description of the implementation of Afivo.
First, the properties of *orthtree* meshes are discussed.
Three data types are used to store these meshes: boxes, levels and trees. These
data types are described below.

|       |        | 
| ----: |  :---- |
<img src="../../figures/quadtree_cex1.png" width=100px /> | <img src="../../figures/quadtree_cex2.png" width=100px /> 
<img src="../../figures/quadtree_cex3.png" width=100px /> | <img src="../../figures/quadtree_cex4.png" width=100px /> |

<a name="fig_example-quadtree" />
<img src="../../figures/box_indices.png" width=200px /> 
**Figure 4**. **TODO**: refine in one corner Upper: Example of a quadtree mesh
that gets refined. Here boxes contain \f$2 \times 2\f$ cells, and different
boxes have different colors.
Below: The spatial indices of the boxes. When a box with indices \f$(i,j)\f$ is refined,
its children have indices \f$(2i-1,2j-1)\f$ up to \f$(2i,2j)\f$.}

\subsection sect_quad-octree Orthtree meshes

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
index can be used to point to a box, see section \ref sect_tree-type below.


\subsection sect_box-type Box data type

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
<img src="../../figures/children_neighbors.png" width=200px /> | <img src="../../figures/location_cc_fx.png" width=250px />


**Figure 5**. a) Each box contains an array of children and neighbors.
The ordering of these arrays is shown here for a 2D box.
Red: children, blue: neighbors.
b) Location and indices of the cell-centered variables (black dots) and the
face-centered variables in the x-direction (red dots) for a box of
\f$2\times 2\f$ cells.}

\subsection sect_level-type Level data type

The *level* data type contains three lists:

	- A list with all the boxes at refinement level \f$l\f$
	- A list with the *parents* (boxes that are refined) at level \f$l\f$
	- A list with the *leaves* (boxes that are not refined) at level \f$l\f$

This separation is often convenient, because some algorithms operate only on
leaves while others operate on parents or on all boxes.
These lists contain the integer indices of the boxes in the tree data structure
described below.


\subsection sect_tree-type Tree data type

The tree data type contains all the data of the mesh.
Most importantly, it stores two arrays: one that contains all the boxes and one
that contains all the levels.

\verbatim
Since Afivo is implemented in Fortran, these arrays start at index one.
\endverbatim

Some other information is also stored: the current maximum refinement level, the
number of cells per box-dimension \f$N\f$, the number of face and cell-centered
variables and the grid spacing \f$\Delta x\f$ at the coarsest level.


\section sect_methods Methods

In this section we give a brief overview of the most important methods in Afivo.
The names of the methods for two-dimensional meshes are used, which have the
prefix <code>a2_</code>. The three-dimensional analogs have, not surprisingly, a prefix <code>a3_</code>.

\subsection sect_init-mesh Creating the initial mesh

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
is discussed in section \ref sect_fill-ghost-cell.

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


\subsection sect_ref-procedure The refinement procedure

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

\subsection sect_fill-ghost-cell Filling of ghost cells

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
\ref sect_one-ghost-cell and \ref sect_no-corner-ghost .

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


\subsection sect_interp-restrict Interpolation and restriction

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
<img src="../../figures/interp_2d.png" width=200px />
**Figure 7**. Schematic drawing of \f$2-1-1\f$ interpolation.
The three nearest coarse grid values are used to interpolate to the center
of a fine grid cell.
Note that the same interpolation scheme can be used for all fine grid cells,
because of the symmetry in a cell-centered discretization.


\subsection sect_tidy-up The list of boxes

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

\subsection sect_output Producing output

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


\section sect_afivo-multigrid Multigrid

Multigrid can be seen as a technique to improve the convergence of a relaxation
method, by using a hierarchy of grids.
Afivo comes with a built-in geometric multigrid solver
<a name="eq_mg-linear-equation">(1)</a>, to solve problems of the form
\f[
  A(u) = \rho,
\f]
where \f$A\f$ is a (nearly) elliptic operator, \f$\rho\f$ the right-hand side and
\f$u\f$ the solution to be computed.
In discretized form, we write the 
<a href="#eq_mg-linear-equation)" >equation(1)</a> above as 
<a name="eq_mg-discr-equation">(2)</a>
\f[
  A_h(u_h) = \rho_h,
\f]
where \f$h\f$ denotes the mesh spacing at which the equation is discretized.

There already exists numerous sources on the foundations of multigrid, the
different cycles and relaxation methods, convergence behaviour and other
aspects, see for example
Xcite{Brandt_2011,trottenberg2000multigrid,Briggs_2000,Hackbusch_1985}.
Here we will not provide a general introduction to multigrid.
Instead we briefly summarize the main ingredients, and focus on one particular
topic: how to implement multigrid on an adaptively refined quadtree or octree
mesh.
On such a mesh, the solution has to be specified on all levels.
Therefore we use FAS multigrid, which stands for Full Approximation Scheme.
Below, the implementation of the various multigrid components in Afivo are
described.


\subsection sect_mg-v-cycle The V-cycle

Suppose there are levels \f$l = l_\mathrm{min}, l_\mathrm{min}+1, \ldots, l_\mathrm{max}\f$, then the FAS
V-cycle can be described as

	1. For \f$l\f$ from \f$l_\mathrm{max}\f$ to \f$l_\mathrm{min}+1\f$, perform
	\f$N_\mathrm{down}\f$ relaxation steps on level \f$l\f$, then update level \f$l-1\f$ (see
	below)
	2. Perform \f$N_\mathrm{base}\f$ relaxation steps on level \f$l_\mathrm{min}\f$, or apply a direct solver
	3. For \f$l\f$ from \f$l_\mathrm{min}+1\f$ to \f$l_\mathrm{max}\f$, perform a
	correction using the data from level \f$l-1\f$
	<a name="eq_coarse-corr" >(3)</a>
	\f[
	  u_h \leftarrow u_h + I_H^h(v_H - v'_H),
	\f]
	then perform \f$N_\mathrm{up}\f$ relaxation steps on level \f$l\f$.
	(See below for the notation)

The first two steps require some extra explanation.
Let us denote the level \f$l-1\f$ grid by \f$H\f$ and the level \f$l\f$ grid by \f$h\f$, and let
\f$v\f$ denote the current approximation to the solution \f$u\f$.
Furthermore, let \f$I_H^h\f$ be an interpolation operator to go from coarse to fine
and \f$I_h^H\f$ a restriction operator to go from fine to coarse.
For these operators, the schemes described in section
\ref sect_interp-restrict are used.

In the first step, the coarse grid is updated in the following way

	1. Set \f$v_H \leftarrow I_h^H v_h\f$, and store a copy \f$v'_H\f$ of \f$v_H\f$
	2. Compute the fine grid residual \f$r_h = \rho_h - A_h(v_h)\f$
	3. Update the coarse grid right-hand side
	<a name="eq_coarse-rhs" >(4)</a>
	\f[
	  \rho_H \leftarrow I_h^H r_h + A_H(v_H).
	\f]

This last equation can also be written as
\f[
  \rho_H \leftarrow I_h^H \rho_h + \tau_h^H,
\f]
where \f$\tau_h^H\f$ is given
by Xcite{Brandt_2011,trottenberg2000multigrid,Briggs_2000,Bai_1987,Hackbusch_1985}
<a name="eq_mg-tau">(5)</a>
\f[
  \tau_h^H = A_H(I_h^H v_h) - I_h^H A_h(v_h).
\f]
This term can be seen as a correction to \f$\rho\f$ on the coarse grid.
When a solution \f$u_h\f$ is found such that \f$A_h(u_h) = \rho_h\f$, then
\f$u_H = I_h^H u^h\f$ will be a solution to \f$A_H(u_H) = \rho_H\f$.

In the second step, relaxation takes place on the coarsest grid. In order to
quickly converge to the solution with a relaxation method, this grid should
contain very few points (e.g., \f$2\times 2\f$ or \f$4\times 4\f$ in 2D). Alternatively,
a direct solver can be used on the coarsest grid, in which case it can contain more points.
Such a direct method has not yet been built into Afivo, although this is
planned for the future. As a temporary solution, additional coarse grids can be
constructed below the coarsest quadtree/octree level. For example, if a quadtree
has boxes of \f$16\times 16\f$ cells, then three levels can be added below it
(\f$8\times 8\f$, \f$4\times 4\f$ and \f$2\times 2\f$), which can then be used in the
multigrid routines.


\subsection sect_mg-fmg-cycle The FMG-cycle

The full multigrid (FMG) cycle that is implemented in Afivo works in the
following way:

	1. If there is no approximation of the solution yet, then set the initial
	guess to zero on all levels, and restrict \f$\rho\f$ down to the coarsest grid
	using \f$I_h^H\f$.
	If there is already an approximation \f$v\f$ to the solution, then restrict \f$v\f$
	down to the coarsest level.
	Use <a href="#eq_coarse-rhs"> equation(4)</a> to set \f$\rho\f$ on coarse grids.
	2. For \f$l = l_\mathrm{min}, l_\mathrm{min}+1, \ldots, l_\mathrm{max}\f$

		* Store the current approximation \f$v_h\f$ as \f$v`_h\f$
		* If \f$l > l_\mathrm{min}\f$, perform a coarse grid correction using
		<a href="#eq_coarse-corr" >equation(3)</a>
		* Perform a V-cycle starting at level \f$l\f$, as described in the previous
		section

\subsection sect_gsrb Gauss Seidel red-black

In Afivo, we have implemented Gauss Seidel red-black or GS-RB as a relaxation
method.
This method is probably described in almost all textbooks on multigrid, such as
Xcite{Brandt_2011,trottenberg2000multigrid,Briggs_2000,Hackbusch_1985}, so we
just give a very brief description.

The *red-black* refers to the fact that points are relaxed in an
alternating manner, using a checkerboard-like pattern. For example, in two
dimensions with indices \f$(i,j)\f$ points can be labeled *red* when \f$i+j\f$ is
even and *black* when \f$i+j\f$ is odd. Now consider 
<a href="#eq_mg-discr-equation">equation(2)</a>, which typically relates a value \f$u_h^{(i,j)}\f$ to
neighboring values and the source term \f$\rho\f$. If we keep the values of the
neighbors fixed, then we can determine the value \f$u_h^{(i,j)}\f$ that locally
solves the linear equation. This is precisely what is done in GS-RB: the linear
equations are solved for all the red points while keeping the old black values,
and then vice-versa.


\subsection sect_mg-ghost-cells Conservative filling of ghost cells

The finer levels will typically not cover the complete grid in Afivo, so that
ghost cells have to be used near refinement boundaries.
These ghost cells can be filled in multiple ways, which will affect the
multigrid solution and convergence behavior.
Here we consider *conservative* schemes for filling ghost cells
Xcite{Bai_1987,trottenberg2000multigrid}.
A conservative scheme ensures that the coarse flux across a refinement boundary
equals the average of the fine fluxes, see <a href="#fig_mg-ref-bound" >Figure 8</a>.

Ensuring consistent fluxes near refinement boundaries helps in obtaining a
*consistent* solution. For example, if we consider a general equation of
the form \f$\nabla \cdot \vec{F} = \rho\f$, then the divergence theorem gives
\f[
  \int_V \rho \, dV = \int_V \nabla \cdot \vec{F} \, dV = \int \vec{F} \cdot
  \vec{n} \, dS,
\f]
where the last integral runs over the surface of the volume \f$V\f$, and \f$\vec{n}\f$
is the normal vector to this surface.
This means that when fine and coarse fluxes are consistent, the integral over
\f$\rho\f$ will be same on the fine and the coarse grid.

<a name="fig_mg-ref-bound" />
<img src="../../figures/mg_refinement_boundary.png" width=300px />
**Figure 8**. Two coarse cells, of which the right one is refined.
The cell centers are indicated by dots.
There are two ghost values (red dots) on the left of the refinement boundary.
Fluxes across the refinement boundary are indicated by arrows.

The construction of a conservative scheme for filling ghost cells is perhaps
best explained with an example.
Consider a 2D Poisson problem
\f[
  \nabla^2 u = \nabla \cdot (\nabla u) = \rho,
\f]
with a standard 5-point stencil for the Laplace operator
<a name="eq_5-point-stencil" >(6)</a>
\f[
  L_h = h^{-2}
  \begin{bmatrix}
    & \;1 &\\
    1 & -4 & 1\\
    & \;1 &
  \end{bmatrix}.
\f]
With this stencil, the coarse flux \f$f_H\f$ across the refinement boundary in
 <a href="#fig_mg-ref-bound" >Figure 8</a> is given by
\f[
  f_H = [u_H^{(2,1)} - u_H^{(1,1)}]/H,
\f]
and on the fine grid, the two fluxes are given by
\begin{align}
  f_{h,1} &= [u_h^{(3,1)} - g_h^{(2,1)}]/h,\\
  f_{h,2} &= [u_h^{(3,2)} - g_h^{(2,2)}]/h.
\end{align}
The task is now to fill the ghost cells \f$g_h^{(2,1)}\f$ and \f$g_h^{(2,2)}\f$ in such
a way that the coarse flux equals the average of the fine fluxes, i.e., such
that
<a name="eq_mg-flux-cons" >(7)</a>
\f[
  f_H = (f_{h,1} + f_{h,2})/2
\f]
To relate \f$u_H^{(2,1)}\f$ to the refined values \f$u_h\f$, the restriction operator
\f$I_h^H\f$ needs to be specified.
In our implementation, this operator does averaging over the children, which can
be represented as
\f[
  I_h^H = \frac{1}{4}
  \begin{bmatrix}
    1 & 1\\
    1 & 1
  \end{bmatrix}.
\f]
The constraint from <a href="#eq_mg-flux-cons">equation(7)</a> can then be written as
\f[
  g_h^{(2,1)} + g_h^{(2,2)} = u_H^{(1,1)} + \frac{3}{4} \left(u_h^{(3,1)} + u_h^{(3,2)}\right)
  - \frac{1}{4} \left(u_h^{(4,1)} + u_h^{(4,2)}\right).
\f]
Any scheme for the ghost cells that satisfies this constraint will be a
conservative discretization.

Bilinear *extrapolation* (similar to standard bilinear interpolation) gives
the following scheme for \f$g_h^{(2,1)}\f$
\f[
  g_h^{(2,1)} = \frac{1}{2} u_H^{(1,1)} + \frac{9}{8} u_h^{(3,1)} -
  \frac{3}{8} \left (u_h^{(3,2)} + u_h^{(4,1)} \right)
  + \frac{1}{8} u_h^{(4,2)}.
\f]
(The scheme for \f$g_h^{(2,2)}\f$ should then be obvious.)
Another option is to use only the closest two neighbors for the extrapolation,
which gives the following expression for \f$g_h^{(2,1)}\f$
<a name="eq_ghost-cell-standard-2d" >(8)</a>
\f[
  g_h^{(2,1)} = \frac{1}{2} u_H^{(1,1)} + u_h^{(3,1)} -
  \frac{1}{4} \left (u_h^{(3,2)} + u_h^{(4,1)} \right).
\f]
This last scheme is how refinement-boundary ghost cells are filled by default in
Afivo.


\subsubsection sect_3d_case Three-dimensional case

In three spatial dimensions, the 5-point stencil of
<a href="#eq_5-point-stencil">equation(6)</a> becomes a 7-point stencil with \f$-6\f$ at the center, and
the restriction operator has eight entries of 1/8.
The analog of <a href="#eq_ghost-cell-standard-2d">equation(8)</a> then becomes
<a name="eq_ghost-cell-standard-3d" >(9)</a>
\f[
  g_h^{(2,1,1)} = \frac{1}{2} u_H^{(1,1,1)} + \frac{5}{4} u_h^{(3,1,1)} -
  \frac{1}{4} \left (u_h^{(4,1,1)} + u_h^{(3,2,1)} + u_h^{(3,1,2)}\right).
\f]


\subsubsection sect_mg-varepsilon Change in \f$\varepsilon\f$ at cell face

**TODO**: math symbol \varepsilon is not allowed in (subsub)section heading.

For the more general equation with a coefficient \f$\varepsilon\f$
\f[
  \nabla \cdot (\varepsilon \nabla u) = \rho,
\f]
we consider a special case: \f$\varepsilon\f$ jumps from \f$\varepsilon_1\f$ to
\f$\varepsilon_2\f$ at a cell face.
Local reconstruction of the solution shows that a gradient \f$(\phi_{i+1} -
\phi_i) / h\f$ has to be replaced by
\f[
  \frac{2 \, \varepsilon_1 \varepsilon_2}{\varepsilon_1 \varepsilon_2} \, \frac{\phi_{i+1} -
    \phi_i} {h},
\f]
or in other words, the gradient is multiplied by the harmonic mean of the
\f$\varepsilon\f$`s (see for example chapter 7.7 of
Xcite{trottenberg2000multigrid}).
The 5-point stencil for the Laplacian can be modified accordingly.

When a jump in \f$\varepsilon\f$ occurs on a coarse cell face, it will also be
located on a fine cell face, see <a href="#fig_mg-ref-bound" >Figure 8</a>.
In this case, the ghost cell schemes described above for constant \f$\varepsilon\f$
still ensure flux conservation.
The reason is that the coarse and fine flux are both weighted by a factor
\f$2 \, \varepsilon_1 \varepsilon_2 / (\varepsilon_1 \varepsilon_2)\f$.


\subsubsection sect_mg-cyl Cylindrical case

In cylindrical coordinates, the Laplace operator can be written as
\f[
  \nabla^2 u = \frac{1}{r} \partial_r (r \partial_r u) + \partial^2_z u
  = \partial_r^2 u + \frac{1}{r} \partial_r u + \partial_z^2 u,
\f]
where we have assumed cylindrical symmetry (no \f$\phi\f$ dependence).
At a radius \f$r \neq 0\f$, the 5-point stencil is
<a name="eq_5-point-stencil-cyl" >(10)</a>
\f[
  L_h = h^{-2}
  \begin{bmatrix}
    & \;1 &\\
    1-\frac{h}{2 r} & -4 & 1+\frac{h}{2 r}\\
    & \;1 &
  \end{bmatrix}.
\f]
With the cell-centered grids in Afivo, radial grid points are located at
\f$(i - \frac{1}{2}) h\f$ for \f$i = 1, 2, 3, \ldots\f$, which means we do not have to
consider the special case \f$r = 0\f$.
For this type of grid indexing, the 5-point stencil can also be written as
\f[
  L_h = h^{-2}
  \begin{bmatrix}
    & \;1 &\\
    \frac{2i-2}{2i-1} & -4 & \frac{2i}{2i-1}\\
    & \;1 &
  \end{bmatrix}.
\f]

If we do not modify the restriction operator, then the ghost cells can still be
filled with the schemes from <a href="#eq_ghost-cell-standard-2d">equations(8)</a> and
<a href="#eq_ghost-cell-standard-3d">(9)</a>.
One way to interpret this is that fluxes are computed in the same way in
cylindrical coordinates, although their divergence is weighted by the radius:
\f[
  \nabla \cdot \vec{F} = \frac{1}{r} \partial_r (r F_r) + \ldots
\f]
From <a href="#fig_mg-ref-bound">Figure 8</a>, we can see that for refinement in the
\f$r\f$-direction, the coarse and fine flux are `weighted` by the same radius.
For the fluxes in the \f$z\f$-direction, the computations are the same as for the
Cartesian case.

\verbatim
Note that when the restriction operator is changed to include radial weighting,
these arguments are no longer valid.
\endverbatim

\subsection sect_mg-examples Multigrid test problems

In this section we present several Poisson test problems
to demonstrate the multigrid behavior on a partially refined mesh.
We use the **method of manufactured solutions**: from an analytic solution the right-hand side and boundary
conditions are computed. Two test problems are considered, a
constant-coefficient
<a class="code" href="poisson__basic__2d_8f90.html#aef02b53cac21b72a47afc1c34286c443">two-</a> and
<a class="code" href="poisson__basic__3d_8f90.html#a1b0ee6adabd66bb15df1120e11776ead">three-dimensional</a>
Poisson equation shown in <a href="examples.html"><span>Examples</span></a>
<a name="eq_mg-example-lpl-1" >(11)</a>
\f[
  \nabla^2 u = \nabla \cdot (\nabla u) = \rho
\f]
and also a 
<a class="code" href="poisson__cyl_8f90.html#aaa62a5d2b70da45c7561014a92ccb9f0">cylindrical</a>
version with a coefficient \f$\varepsilon\f$
<a name="eq_mg-example-lpl-2" >(12)</a>
\f[
  \frac{1}{r} \partial_r (r \varepsilon \partial_r u) + \partial_z (\varepsilon \partial_z u) = \rho,
\f]
both on a two-dimensional rectangular domain \f$[0,1] \times [0,1]\f$.
For the cylindrical case, \f$\varepsilon\f$ has a value of \f$100\f$ in the lower left
quadrant \f$[0,0.25] \times [0,0.25]\f$, and a value \f$1\f$ in the rest of the domain.
In all cases, we pick the following solution for \f$u\f$
\f[
  u(r) = \exp(|{\vec{r}-\vec{r}_1}|/\sigma) + \exp(|{\vec{r}-\vec{r}_2}|/\sigma),
\f]
where \f$\vec{r_1} = (0.25, 0.25)\f$, \f$\vec{r_2} = (0.75, 0.75)\f$ and
\f$\sigma = 0.04\f$.
An analytic expression for the right-hand side \f$\rho\f$ is obtained by plugging
the solution in <a href="#eq_mg-example-lpl-1">equation(11)</a> and
<a href="#eq_mg-example-lpl-2">equation(12)</a>.
Note that jumps in \f$\varepsilon\f$ also contribute to the source term \f$\rho\f$,
and the solution itself is used to set boundary conditions.

The two different problems can now be solved numerically.
For the cylindrical case with the varying \f$\varepsilon\f$, a modified Laplacian
operator is used, as described in section \ref sect_mg-ghost-cells.
The Gauss Seidel red-black relaxation methods are also modified, because they
depend on the applied operator, see section \ref sect_gsrb.
For these examples, we have used
\f$N_\mathrm{down} = N_\mathrm{up} = N_\mathrm{base} = 2\f$ (number of down/up/base
smoothing steps), and a coarsest grid of \f$2\times 2\f$ cells.

It is possible to do adaptive mesh refinement in multigrid, for example by using
an estimate of the local truncation error based on <a href="#eq_mg-tau">equation(5)</a>
(see also chapter 9 of Xcite{trottenberg2000multigrid}).
Such a technique is not used here, instead the refinement criterion is based on
the right-hand side: refine if \f$\Delta x^2 |\rho| > 5.0 \times 10^{-4}\f$.
The resulting mesh spacing is shown in <a href="#fig_mg-ex1">Figure 9a</a>.

In both cases, one FMG (full multigrid) cycle is enough to achieve convergence
up to the discretization error, which was approximately \f$10^{-4}\f$ for the mesh
of <a href="#fig_mg-ex1">Figure 9a</a>.
Consecutive FMG cycles have a negligible effect on the absolute error, although
they do reduce the residual \f$r = \rho - \nabla^2 u\f$.
The maximum value of \f$|r|\f$ is shown versus iteration number in
<a href="#fig_mg-ex1">Figure 9b</a>.
The convergence behaviour is similar for both cases, with each iteration
reducing the residual by a factor of about \f$0.056\f$.
The offset between the lines is caused by the \f$\varepsilon = 100\f$ region, which
locally amplifies the source term by a factor of 100.

<a name="fig_mg-ex1" />
a) Mesh spacing | b) Maximum residual versus FMG
------------- | -------------
<img src="../../figures/mesh_mg_ex_1.png" width=300px /> | <img src="../../figures/ex1_mg_res_pdf.tex" width="300px" />

**Figure 9**. a) mesh spacing used for the multigrid examples, in a
\f$[0,1] \times [0,1]\f$ domain.
Each step in color is a factor two in refinement, with red indicating
\f$\Delta x = 2^{-5}\f$ and the darkest blue indicating \f$\Delta x = 2^{-12}\f$.
b) the maximum residual versus FMG iteration, case 1 corresponds to
<a href="#eq_mg-example-lpl-1">equation(11)</a> and case 2 to 
<a href="#eq_mg-example-lpl-2">equation(12)</a>.

\section sect_example-fluid Implementing a plasma fluid model


To illustrate how Afivo can be used, we describe the implementation of a simple
2D/3D *plasma fluid model* for streamer discharges below. For simplicity,
photoionization is not included in this example. A review of fluid models for
streamer discharges can be found in Xcite{Luque_2012}.


\subsection sect_model-formulation Model formulation

We use the so-called drift-diffusion-reaction approximation:
<a name="eq_fluid-model" >(13)</a>
\f[
\begin{align}
  \partial_t n_e &= -\nabla \cdot \vec{j}_e + \bar{\alpha} |{\vec{j}_e}|,\\
  \partial_t n_i &= \bar{\alpha} |{\vec{j}_e}|,\\
  \vec{j}_e &= -\mu_e n_e \vec{E} - D_e \nabla n_e,
\end{align}
\f]

where \f$n_e\f$ is the electron density, \f$n_i\f$ the positive ion density, \f$\vec{j}_e\f$
the electron flux, \f$\bar{\alpha}\f$ the effective ionization coefficient, \f$\mu_e\f$
the electron mobility, \f$D_e\f$ the electron diffusion coefficient and \f$\vec{E}\f$
the electric field. The above equations are coupled to the electric field, which
we compute in the electrostatic approximation:
\f[
\begin{align}
  \vec{E} &= -\nabla \phi,\\
  \nabla^2 \phi &= -\rho/\varepsilon_0\\
  \rho &= e (n_i - n_e),
\end{align}
\f]
where \f$\phi\f$ is the electric potential, \f$\varepsilon_0\f$ the permittivity of
vacuum and \f$e\f$ the elementary charge. The electric potential is computed with
the multigrid routines described in section \ref sect_afivo-multigrid.

We make use of the *local field approximation* Xcite{Li_2007}, so that
\f$\mu_e\f$, \f$D_e\f$ and \f$\bar{\alpha}\f$ are all functions of the local electric field
strength \f$E = |{\vec{E}}|\f$. These coefficients can be obtained
experimentally, or they can be computed with a Boltzmann solver
Xcite{Bolsighagelaar011,Dujko_2011} or particle swarms Xcite{Li_hybrid_i_2010}.


\subsection sect_flux-calc-time-stepping Flux calculation and time stepping

The electron flux is computed as in Montijn \cite Montijn_2006. For the diffusive part,
we use central differences. The advective part is computed using the Koren
limiter Xcite{koren_limiter}. The Koren limiter was not designed to include
refinement boundaries, and we use linear interpolation to obtain fine-grid ghost
values. These ghost cells lie inside a coarse-grid neighbor cell, and we limit
them to twice the coarse values to preserve positivity. (We would like to
improve this in the future.)

Time stepping is also performed as in Montijn \cite Montijn_2006, using the explicit
trapezoidal rule, also known as the modified Euler`s method. The time step is
determined by a CFL condition for the electron flux and the dielectric
relaxation time, as in Montijn \cite Montijn_2006.

\subsection sect_refinement-crit Refinement criterion

Our refinement criterion contains two components: a *curvature monitor*
\f$c_\phi\f$ for the electric potential and a monitor \f$\bar{\alpha} \Delta x\f$ which
gives information on how well the ionization length (\f$1/\bar{\alpha}\f$) is
resolved. For both, we use the maximum value found in a box in order to decide
whether to (de)refine it.

Since \f$\nabla^2 \phi = -\rho / \varepsilon_0\f$, the curvature monitor can be
computed as \f$c_\phi = \Delta x^2 |\rho| / \varepsilon_0\f$. The quantity
\f$\bar{\alpha} \Delta x\f$ is computed by locating the highest electric field in
the box, and looking up the corresponding value of \f$\bar{\alpha}\f$. The combined
refinement criterion is then as follows, where later rules can override earlier
ones:
	- If \f$\bar{\alpha} \Delta x < 0.1\f$ and \f$\Delta x < 25 \, \mu\textrm{m}\f$,
	derefine.
	- If \f$t < 2.5 \, \textrm{ns}\f$, ensure that there is enough refinement
	around the initial seed to resolve it.
	- If \f$\bar{\alpha} \Delta x > 1.0\f$ and \f$c_\phi > 0.1 \, \textrm{Volt}\f$,
	refine.

\subsection sect_results Simulation conditions and results

<a name="fig_ex-streamer-seed" />
<img src="../../figures/streamer_seed.png" width=300 />
**Figure 10**. Cross section through the center of the three-dimensional simulation
domain. The ionized seeds with a density of \f$10^{20} \, \textrm{m}^{-3}\f$
electrons and ions are indicated in red. There is a background density of
\f$5 \times 10^{15} \, \textrm{m}^{-3}\f$ electrons and ions, and the background
electric field points down with a magnitude \f$E_0 = 2.5 \, \textrm{MV/m}\f$.

A cross section through the computational domain of \f$(32 \, \textrm{mm})^3\f$ is
shown in <a href="#fig_ex-streamer-seed" >Figure 10</a>. The background field points down,
with a magnitude \f$E_0 = 2.5 \, \textrm{MV/m}\f$, which is about
\f$5/6\f$\textsuperscript{th} of the critical field. The background field is applied
by grounding the bottom boundary of the domain, and applying a voltage at the
top. At the other sides of the domain we use Neumann boundary conditions for the
potential. We use transport coefficients (e.g., \f$\bar{\alpha}\f$, \f$\mu_e\f$) for
atmospheric air, but for simplicity photoionization has not been included.
Instead a background density of \f$5 \times 10^{15} \, \textrm{m}^{-3}\f$ electrons
and positive ions is present.

Two seeds of electrons and ions locally enhance the background electric field,
see <a href="#fig_ex-streamer-seed" >Figure 10</a>. These seeds have a density of
\f$10^{20} \, \textrm{m}^{-3}\f$, a width of about \f$0.3 \, \textrm{mm}\f$ and a length
of \f$1.6 \, \textrm{mm}\f$. The electrons from these seeds will drift upwards,
enhancing the field at the bottom of the seed where a positive streamer can
form.

In <a href="#fig_ex-streamer-dns">Figure 11</a>, the time evolution of the electron
density is shown, and in <a href="#fig_ex-streamer-fld">Figure 12</a>. the electric field
is shown. Two positive streamers grow downwards from the ionized seeds. The
upper one is attracted to the negatively charged end of the lower one, and
connects to it at around \f$9.5 \, \textrm{ns}\f$. The three-dimensional simulation
took about 3.5 hours on a 16-core machine, and eventually used about
\f$1.3 \times 10^7\f$ grid cells.

<a name="fig_ex-streamer-dns" />
| 1 ns  | 3 ns   | 5 ns   | 7 ns   | 9 ns   | 11 ns  | 13 ns  | 15 ns  | 17 ns  | |
| ----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :---- |
<img src="../../figures/visit/crop_two_str0010.png" width=60px /> | <img src="../../figures/visit/crop_two_str0011.png" width=60px /> | <img src="../../figures/visit/crop_two_str0012.png" width=60px /> | <img src="../../figures/visit/crop_two_str0013.png" width=60px /> | <img src="../../figures/visit/crop_two_str0014.png" width=60px /> | <img src="../../figures/visit/crop_two_str0015.png" width=60px /> | <img src="../../figures/visit/crop_two_str0016.png" width=60px /> | <img src="../../figures/visit/crop_two_str0017.png" width=60px /> | <img src="../../figures/visit/crop_two_str0018.png" width=60px /> | <img src="../../figures/legend_elec_3d.png" width=80px /> 
**Figure 11**. A three-dimensional simulation showing two positive streamers
propagating downwards. The upper one connects to the back of the lower one.
The electron density is shown using volume rendering, for which the opacity
is indicated in the legend; low densities are transparent

<a name="fig_ex-streamer-fld" />
| 1 ns  | 3 ns   | 5 ns   | 7 ns   | 9 ns   | 11 ns  | 13 ns  | 15 ns  | 17 ns  | |
| ----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :----: | :---- |
<img src="../../figures/visit/crop_two_str0000.png" width=60px /> | <img src="../../figures/visit/crop_two_str0001.png" width=60px /> | <img src="../../figures/visit/crop_two_str0002.png" width=60px /> | <img src="../../figures/visit/crop_two_str0003.png" width=60px /> | <img src="../../figures/visit/crop_two_str0004.png" width=60px /> | <img src="../../figures/visit/crop_two_str0005.png" width=60px /> | <img src="../../figures/visit/crop_two_str0006.png" width=60px /> | <img src="../../figures/visit/crop_two_str0007.png" width=60px /> | <img src="../../figures/visit/crop_two_str0008.png" width=60px /> | <img src="../../figures/legend_2d.pdf" width=80px /> 
**Figure 12**. Cross section through the three-dimensional domain showing the time
evolution of the electric field. The full height of the domain is shown
(\f$32 \, \textrm{mm}\f$), but only \f$6 \, \textrm{mm}\f$ of the width.

\section sect_References References

\bibliography{../afivo}{}
