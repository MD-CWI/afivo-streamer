# Grid refinement

In order to adjust the refinement of the mesh, the user can call
m_a2_core::a2_adjust_refinement() (m_a3_core::a3_adjust_refinement() in 3D).
Each call will increase or decrease the refinement level at any location by at
most one level. An example from @ref advection_2d.f90 is shown below:

\snippet advection_2d.f90 adjust_refinement

The last argument specifies that refinement within a distance of 2 cells of the
box's boundary should also trigger refinement in neighboring boxes. After
refinement, information about the changes is returned through the `refine_info`
argument. The argument `ref_routine` is a user-supplied subroutine that will be
called for eac *relevant* box (boxes covered with multiple layers of refinement
don't matter). For each cell of a box, the subroutine should set a refinement
flag. The options are:

* `af_do_ref`: Extra refinement is required
* `af_rm_ref`: The current refinement can be removed
* `af_keep_ref`: The current refinement should stay

An example from @ref advection_2d.f90 is shown below.

\snippet advection_2d.f90 refine_routine

If any of the cells is marked for refinement, all the cells will be refined.
Refinement will only be removed if all cells are marked for removal, and if the
parent is not marked for refinement.

# Prolongation after refinement

After the mesh has been refined, there will typically be new boxes. It is the user's responsibility to set meaningful values on these new boxes. An example from @ref advection_2d.f90 is shown below.

\snippet advection_2d.f90 prolong_to_new_children

