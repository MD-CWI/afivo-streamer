# Initial conditions

# Using a standard seed

Initial seeds in the form of line segments are supported. These seeds can be defined with the following parameters, where `N` refers to the number of seeds the user wants to define:

name | value | description
---|---|---
`seed_density` | N * real | Maximum density for the seed (in 1/m3)
`seed_rel_r0` | N * NDIM reals | relative start position (between 0 and 1)
`seed_rel_r1` | N * NDIM reals | relative end position (between 0 and 1)
`seed_charge_type` | N * integer | 0 for neutral, 1 for positive ions, -1 for electrons
`seed_falloff` | N * string | Fall-off type for seed (sigmoid, gaussian, smoothstep, step, laser)
`seed_width` | N * real | Seed decay width (m)

Other settings that control the initial conditions are:

name | value | description
---|---|---
`background_density` | real | The background ion and electron density (1/m3)

# Defining custom initial conditions

To define custom initial conditions, specify the following in your `m_user.f90` file:

    user_initial_conditions => my_init_cond

where `my_init_cond` is a routine (see @ref
m_user_methods::user_initial_conditions) that can modify the initial state in
any way you like.

## Finding the index of a species

To set the value of a particular species, first find its index with

    i_species = af_find_cc_variable(tree, "species_name")

then you can for example set species density to zero (in 2D) with

    box%cc(:, :, i_species) = 0.0_dp

## Using a static grid for a specific region in the domain

The user can specify a specific region of the domain to have fixed grid spacing. The following flags are utilized for that purpose:

* `refine_regions_dr`: the minimum grid spacing in the region
* `refine_regions_rmin`: the starting boundary of the region
* `refine_regions_rmax`: the closing boundary of the region

The flag `refine_regions_rmin` should have two values if working in a two dimensional domain and three values if in a three dimensional domain. The same is the case for `refine_regions_rmax`. The values are separated by spaces. For example, having `refine_regions_rmin = 0.0 0.0` and  `refine_regions_rmax = 6.0e-4 5.0e-2` in a cylindrically symmetric simulation corresponds to the region bounded by the lines `r = 0.0`, `z = 0.0`, `r = 6.0e-4`, and `z = 5.0e-2`.
