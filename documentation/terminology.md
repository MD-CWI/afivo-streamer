# Terminology and variable names

Here the most important terminology and variable names are collected.

# Terminology

Term | Description
---|---
AMR | Adaptive Mesh Refinement
quadtree/octree | Afivo uses adaptively refined grids of the [quadtree and octree](@ref documentation/quadtree_octree.md) type
box | A box contains \f$N^D\f$ grid cells, where \f$D\f$ is the problem dimension and \f$N\f$ is e.g. 8 or 16, see @ref documentation/data_structures.md
tree | The full grid, containing all the boxes at all refinement levels
neighbors | A box can have neighbors (adjacent boxes at the same refinement level) in the x, y, and z-direction
parent | When a box is refined, it becomes a parent (since it now has refined 'children')
children | When a box is refined, it is covered by its children, which 4 (2D) or 8 (3D) boxes with half the grid spacing
leaf | A box that is not refined / which does not have children

# Variable names

name | Typically used for
---|---
tree | Data structure containing the full grid
boxes | List of all the boxes, at all refinement levels
id | Integer indicating the index of a boxes in the boxes array
nc | The number of cells along each dimension of a box (a box contains \f$nc^D\f$ cells)
lvl | Integer indicating the refinement level (1, 2, ...)
neighbors | List of \f$2D\f$ integers to indicate the neighbors of a box. They can be positive (index of neighbor in boxes array), zero (refinement boundary) or negative (physical boundary)
parent | An integer with the index of the parent, or zero if there is no parent
children | List of \f$2^D\f$ integers with the indices of the children, and zero if there are no children
`af_no_box` | The integer value zero


