# Important data structures

# Box data type {#main-box-type}

A *box* is the basic mesh unit in Afivo, see `m_a2_types::box2_t` and
`m_a3_types::box3_t`. Each box consists of \f$N^D\f$ grid cells, where \f$N\f$
has to be an even number and \f$D\f$ is the spatial dimension.

Two types of cell data are supported: cell-centered (cc) data and face-centered
(fc) data. For example, in 2D `cc(i, j, n)` selects the nth variable at cell
(i,j), whereas `fc(i, j, dim, n)` selects the nth face-centered variable in
direction dim. A box contains one layer of ghost cells for its cell-centered
variables, are illustrated in the figure below. How these ghost cells are filled
is discussed in @ref documentation/filling_ghost_cells.md.

![Location and indices of the cell-centered variables (black dots) and the face-centered variables in the x-direction (red dots) for a box of 2x2 cells.](location_cc_fx.png)

Boxes store their neighbors, their children and their parent. A special value
`m_afivo_types::af_no_box` (which is zero) is used to indicate that a parent,
child or neighbor does not exist. In the case of neighbors, physical boundary
conditions are specified by `m_afivo_types::af_phys_boundary` (which is -1).

Furthermore, boxes contain some `convenience` information, such at their
refinement level, minimum coordinate and spatial index.


# Level data type {#main-level-type}

The *level* data type (see `m_afivo_types::lvl_t`) contains three lists:

- A list with all the boxes at refinement level \f$l\f$
- A list with the *parents* (boxes that are refined) at level \f$l\f$
- A list with the *leaves* (boxes that are not refined) at level \f$l\f$

This separation is often convenient, because some algorithms operate only on
leaves while others operate on parents or on all boxes. These lists contain the
integer indices of the boxes in the tree data structure described below.

# Tree data type {#main-tree-type}

The tree data type contains all the data of the mesh, see `m_a2_types::a2_t` and
`m_a3_types::a3_t`. Most importantly, it stores two arrays: one that contains
all the boxes and one that contains all the levels.

Some other information is also stored: the current maximum refinement level, the
number of cells per box-dimension \f$N\f$, the number of face and cell-centered
variables and the grid spacing \f$\Delta x\f$ at the coarsest level.


