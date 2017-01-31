# Important data structures

* Box
* Level
* Tree

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


