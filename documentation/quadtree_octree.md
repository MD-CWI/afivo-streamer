# Quadtree and octree grids

Afivo uses so-called quadtree (2D) or octree (3D) grids. A quadtree grid in
Afvio consists of boxes (i.e., blocks) of $N \times N$ cells, with $N$ an
even number. A user can for example select to use boxes of $4 \times 4$ cells, see @ref documentation/coarse_grid.md.

A box in a quadtree grid can be refined by adding four \emph{child} boxes.
These children contain the same number of cells but half the grid spacing, so
that they together have the same area as their parent. Each of the children can
again be refined, and so on, as illustrated in figure
\ref{fig:quadtree-comp-domain}. There can be up to 30 refinement levels in
Afivo. So-called \emph{proper nesting} or \emph{2:1 balance} is ensured, which
means that neighboring boxes differ by at most one refinement level.

The coarse grid, which defines the computational domain, can then be constructed
from one or more of these boxes, see figure~\ref{fig:quadtree-comp-domain} and
section \ref{sec:coarse-grid}.

Two types of variables are stored: cell-centered variables and face-centered
variables. When initializing Afivo, the user has to set the number of these
variables. For the cell-centered variables, each box has one layer of ghost
cells, as discussed in section~\ref{sec:ghost-cell}.

  \includegraphics[width=0.9\textwidth]{visit/mesh_example.pdf}
  \caption{Left: example of a coarse grid consisting of two boxes of
    $4 \times 4$ cells. The middle and right figure show how the boxes can be
    refined, by covering them with four `child' boxes.}

