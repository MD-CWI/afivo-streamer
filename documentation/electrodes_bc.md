# Electrodes and boundary conditions {#electrodes-bc}

[TOC]

# Boundary conditions

## Electric potential boundary conditions {#bc-phi}

There is a parameter `field_bc_type`, that can be set to:

* 'homogeneous': Neumann zero on the sides, zero voltage at the bottom, and a fixed applied voltage on top.
* 'neumann': Neumann zero on the sides, zero voltage at the bottom, and a Neumann boundary at the top corresponding to the applied electric field

To have more flexibility, a custom routine for boundary conditions can be defined in the `m_user.f90` file, which will override the above setting

## Species boundary conditions {#bc-species}

There is a parameter `species_boundary_condition`, that can be set to:

* 'neumann_zero': Neumann zero for all species on all domain boundaries. This means that electrons can flow out of electrodes.
* 'dirichlet_zero': Dirichlet zero for all species on all domain boundaries.

# Specifying an electrode

To enable an electrode, set `use_electrode = T`. The electrode can be controlled with the following parameters:

\snippet m_field.f90 electrode_settings

The parameter `field_electrode_type` can be set to:

* `rod` (default): a cylindrical electrode with a semi-spherical tip
* `rod_cone_top`: a cylindrical electrode with a conical tip

For the rod electrode, the relative coordinates `field_rod_r0` and `field_rod_r0` need to be defined. For the conical tip, `field_rod_r2` also needs to be specified.