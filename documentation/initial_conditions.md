# Initial conditions {#initial-conditions}

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
