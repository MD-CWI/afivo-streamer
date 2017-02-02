# Quadtree and octree grids

# Introduction

In Afivo so-called [quadtree](https://en.wikipedia.org/wiki/Quadtree) (2D)
or [octree](https://en.wikipedia.org/wiki/Octree) (3D) grids are used.

# Quadtree grids

A quadtree grid in Afvio consists of boxes (i.e., blocks) of \f$N \times N\f$
cells, with \f$N\f$ an even number. An example of a quadtree grid with boxes of
\f$4 \times 4\f$ cells is shown below.

![](mesh_example.png)

A box in a quadtree grid can be refined by covering it with four refined boxes
('children'). These children contain the same number of cells as their parent,
but half the grid spacing. Each of the children can again be refined, as
illustrated above. In Afivo so-called *proper nesting* or *2:1 balance* is
ensured, which means that neighboring boxes differ by at most one refinement
level.

# Octree grids

Octrees are the 3D equivalent of quadtrees. When a box is refined, it is covered
by eight children instead of four in 2D. Some of the other differences are
summarized in the table below.

Property | Quadtree | Octree
---|---|---
Children of a box | 4 | 8
Number of corners | 4 | 8
Number of faces | 0 | 6
Number of edges | 4 | 12


