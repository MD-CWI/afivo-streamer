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

After the mesh has been refined, there will typically be new boxes, and some boxes will have been removed. A user can specify prolongation, restriction and boundary conditions routines for a cell-centered variable with the `m_a2_core::a2_set_cc_methods` routine. For such variables, restriction and prolongation is automatically performed. For other variables, the user is responsible for performing restriction (before the mesh update) and prolongation (after the mesh update). An example is show below:

```{f90}
! Set values on newly added boxes
subroutine prolong_to_new_children(tree, ref_info)
  use m_a2_prolong
  type(a2_t), intent(inout)    :: tree
  type(ref_info_t), intent(in) :: ref_info
  integer                      :: lvl, i, id, p_id

  do lvl = 1, tree%highest_lvl
     ! For each newly added box ...
     do i = 1, size(ref_info%lvls(lvl)%add)
        ! Use prolongation to set the values of variable i_phi
        id = ref_info%lvls(lvl)%add(i)
        p_id = tree%boxes(id)%parent

        call a2_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_phi)
     end do

     ! After values have been set on this level, fill ghost cells
     call a2_gc_ids(tree%boxes, ref_info%lvls(lvl)%add, i_phi, &
          a2_gc_interp, a2_bc_dirichlet_zero)
  end do
end subroutine prolong_to_new_children
```

\snippet advection_2d.f90 prolong_to_new_children

