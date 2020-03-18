# Boundary conditions

# Introduction

Unless a computational domain is fully periodic, there will be physical
boundaries. Afivo implements boundary conditions using ghost cells, see @ref
documentation/filling_ghost_cells.md. Currently, Dirichlet and Neumann type
boundary conditions are supported. An example of how to fill boundary conditions
is given in @ref boundary_conditions.f90.

# Writing a boundary condition routine

The user has to provide routines that will be called near physical boundaries when adding variables through `m_af_core::af_add_cc_variable()`. The interface of these routines is given by `m_af_types::af_subr_bc`. An example is shown below:

\snippet boundary_conditions.f90 boundary_method
