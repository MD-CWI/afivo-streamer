# Defining the coarse grid

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
