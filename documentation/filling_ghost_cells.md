# Filling ghost cells

Each box has storage for a single layer of ghost cells for its cell-centered
variables, which can be filled with one the following routines (and their
equivalents in 3D):

* m_af_ghostcell::af_gc_tree: Fill ghost cells on the complete mesh (all boxes)
* m_af_ghostcell::af_gc_ids: Fill ghost cells on some boxes
* m_af_ghostcell::af_gc_box: Fill ghost cells on a single box

To understand in which order ghost cells are filled, have a look at the implementation of m_af_ghostcell::af_gc_box. Basically, the procedure is as follows:

1. Loop over all the sides of a box. If there is a neighbor, copy ghost cells
   from that neighbors. If there is a refinement boundary, use a
   method for refinement boundaries, and if there is a physical boundary, use
   a method for physical boundaries.
2. Only in 3D: Fill the *edge* ghost cells. Note that a box in 3D has 12 edges.
   If available, copy data from a neighboring box at the same refinement level.
   Otherwise, use extrapolation to fill edge ghost cells.
3. Fill the corner ghost cells. If a neighbor is available at the same
   refinement level, copy data from the neighbor. Otherwise, use extrapolation
   to fill the ghost cells.

# How corners/edges are filled by extrapolation

When data cannot be copied from a neighbor at the same refinement level, we use
linear extrapolation to fill ghost cells. Suppose that there is a corner like this:

    c | d
    -----
    a | b

where a lies inside the box, c and d are ghost cells on the sides, and d is a
corner ghost cell. The recipe for d is then: \f$d = (b+c) - a\f$.

This linear extrapolation has a nice property: if we later perform bilinear
prolongation using the values a, b, c, and d, then we get the same result as
when we linearly interpolated using only a, b, and c.

# Filling a second layer of ghost cells

TODO
