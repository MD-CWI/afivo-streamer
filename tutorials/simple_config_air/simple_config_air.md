# A simple config file for an axisymmetric streamer in air  {#tutorial-air-simple}

Here we will briefly show how to construct a simple config file, leaving most parameters to their default values. An axisymmetric positive streamer will be simulated in air. Below, it is assumed that a file `tutorial_air_simple.cfg` is created in the `programds/standard_2d` folder. The following settings should then be specified in the config file.

## The computational domain

This defines axisymmetric domain measuring 20 mm in the r and z directions:

    # Whether cylindrical coordinates are used (only in 2D):
    cylindrical = T

    # The length of the domain (m):
    domain_len = 20e-3 20e-3

## The simulation name, output time step and end time

Output files with start with `output%%name`. In this case, output will be written every 0.5 ns, and the simulation will run until 10 ns:

    # Name for the output files (e.g. output/my_sim):
    output%name = output/tutorial_air_simple

    # The timestep for writing output (s):
    output%dt = 0.5e-9

    # The desired endtime (s) of the simulation:
    end_time = 10e-9

## The gas

The simulations will be performed in artificial air (80% N2, 20% O2) at 1 bar and 300 K (which are actually the default values):

    # The gas pressure (bar):
    gas%pressure = 1.0

    # Gas component names:
    gas%components = N2 O2

    # Gas component fractions:
    gas%fractions = 0.8 0.2

    # The gas temperature (Kelvin):
    gas%temperature = 300.0

## The electron transport data and reactions to use

There are several transport data and reaction files included with the code, in the `transport_data` folder. In this case, we will use a simple air chemistry:

    # Input file with transport (and reaction) data:
    input_data%file = ../../transport_data/air_light_example_v0.txt

Note that the path should point to the `afivo_streamer/transport_data` directory.

## The background electric field

In this example, the background electric field will be 2.0 MV/m, which is about half of the critical field of air at 300 K and 1 bar:

    # How the electric field or voltage is specified:
    field_given_by = field 2.0e6

The geometry will be plate-to-plate (the default), with a voltage difference applied between the plates.
For other ways of specifying the background field or the applied voltage, see [Electrodes and boundary conditions](documentation/electrodes_bc.md).

## The initial conditions

The formation of a streamer requires:

* A region wherhe the electric field exceeds the critical field
* Some initial electrons

In this example, the background field is below the critical field. To enhance the background field, we will place an elongated conducting channel in the domain (a "seed"), see [Initial conditions](documentation/initial_conditions.md). This seed will have an electron density and positive ion density of `5e19 / m3`, and will be 2 mm long and 0.25 mm wide. Furthermore, we will a so-called "smoothstep" profile (see @ref m_geometry). After some time, the electric field in the seed will be partially screened (i.e., have a lower value), which will enhance the electric field at the endpoints of the seed.

    # Type of seed: neutral (0), ions (1) or electrons (-1)
    seed_charge_type = 0

    # Initial density of the seed (1/m3):
    seed_density = 5e19

    # The relative start position of the initial seed:
    seed_rel_r0 = 0.0 0.45

    # The relative end position of the initial seed:
    seed_rel_r1 = 0.0 0.55

    # Seed width:
    seed_width = 0.25e-3

    # Fallof type for seed, see m_geom.f90:
    seed_falloff = smoothstep

Finally, a small background ionization level of `10^10/m3` electrons and positive ions is included. Such background ionization can help to start the discharge.

    # The background ion and electron density (1/m3):
    background_density = 1e11

## Photoionization

In air, photoionization is an important process, especially for positive streamer discharges. Here we simply use the default Helmholtz photoionization model, see [Photoionization](documentation/photoionization.md) for more details.

    # Whether photoionization is enabled:
    photoi%enabled = T

    # Which photoionization method to use (helmholtz, montecarlo):
    photoi%method = helmholtz

## Running the example and analyzing the output

After putting all the above parameters in a file `tutorial_air_simple.cfg`, the simulation can be performed using the command:

    ./streamer tutorial_air_simple.cfg
