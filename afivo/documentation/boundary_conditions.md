# Boundary conditions

# Introduction

Unless a computational domain is fully periodic, there will be physical
boundaries. Afivo implements boundary conditions using ghost cells, see @ref
documentation/filling_ghost_cells.md. Currently, Dirichlet and Neumann type
boundary conditions are supported. An example of how to fill boundary conditions
is given in @ref boundary_conditions_2d.f90 and @ref boundary_conditions_3d.f90.

# Writing a boundary condition routine

When filling ghost cells, see @ref documentation/filling_ghost_cells.md, the user has to provide a routine that will be called near physical boundaries. The interface of this routine is given by a2_types::a2_subr_bc and a3_types::a3_subr_bc. An example is shown below:

\snippet boundary_conditions_2d.f90 boundary_method

Note that the user method fills ghost cells on the sides of a box, but not
corners. For example, for a boundary condition in the lower-x direction, the
cells with indices (0, 1:nc) need to be filled.
