# Tutorial: including electrodes {#tutorial-electrodes}

## A single rod with a hemispherical cap

In most experiments, a pointed electrode is used to generate streamer discharges. As explained at [electrodes and boundary conditions](documentation/electrodes_bc.md), there are a a few pre-defined electrode shapes that can easily be used. Here we will include a rod with a hemi-spherical cap, placed at the bottom of the domain. The electrode should therefore be grounded (as the bottom of the domain also is). We will perform an axisymmetric simulation, so the electrode should be on the axis. The can be specified with the following parameters:

    # Whether to include an electrode:
    use_electrode = T

    # Whether electrode 1 is grounded or at the applied voltage:
    field_electrode_grounded = T

    # Type of electrode (sphere, rod, rod_cone_top, rod_rod, user):
    field_electrode_type = 'rod'

    # Electrode 1: first relative coordinate:
    field_rod_r0 = 0. 0.

    # Electrode 1: second relative coordinate:
    field_rod_r1 = 0.2 0.

    # Electrode 1 radius (in m):
    field_rod_radius = 0.5e-3

This defines an electrode with a length of 0.2 times the domain height and a radius of 0.5 mm. Another important setting is how fine the mesh around the electrode should *at least* be (it can be refined further according to other refinement criteria). In this case, we specify a maximum grid spacing of 0.025 mm:

    # Ensure grid spacing around electrode is less than this value:
    refine_electrode_dx = 2.5e-5

An example configuration file is provided at `tutorials/including_electrodes/rod_electrode.cfg`, and the resulting initial electric field is shown below.

![Initial field around rod electrode](rod_electrode.png){html: width=50%}

## Two rod electrodes

One of the other built-in options is to include two rod electrodes, for example as follows:

     # Whether to include an electrode:
    use_electrode = T

    # Type of electrode (sphere, rod, rod_cone_top, rod_rod, user):
    field_electrode_type = 'rod_rod'

     # Whether electrode 1 is grounded or at the applied voltage:
    field_electrode_grounded = T

    # Whether electrode 2 is grounded or at the applied voltage:
    field_electrode2_grounded = F

    # Electrode 1: first relative coordinate:
    field_rod_r0 = 0.0 0.0

    # Electrode 1: second relative coordinate:
    field_rod_r1 = 0.0 0.2

    # Electrode 2: first relative coordinate:
    field_rod2_r0 = 0.0 0.8

    # Electrode 2: second relative coordinate:
    field_rod2_r1 = 0.0 1.0

    # Electrode 1 radius (in m):
    field_rod_radius = 0.5e-3

    # Electrode 2 radius (in m):
    field_rod2_radius = 0.5e-3

An example configuration file is provided at `tutorials/including_electrodes/two_rod_electrodes.cfg`, and the initial field is shown below:

![Initial field due to two rod electrodes](two_rod_electrodes.png){html: width=50%}

