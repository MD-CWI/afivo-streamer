# Simulation options

Afivo-streamer uses the
[config_fortran](https://github.com/jannisteunissen/config_fortran) module to
allow users to customize simulations with different options. Each option has a
name given by `[group%]name`, where the group can also be empty. Options can be
specified in configuration files like this:

    group%name = value

Multiple settings from the same group can be given like this:

    [group]
        name = value
        name_2 = value_2

Options can be of type integer, real, logical or string. Some options are arrays,
and they can also have a variable length.

# Examples of a few important options

name | example | meaning
---|---|---
`output%%name` | output/my_sim | filename base for output files
`output%%dt` | 0.25e-9 | time step for writing output
`input_data%%file` | [filename] | input file with transport data and reactions
`end_time` | 10e-9 | end time of the simulation (s)
`domain_len` | 32e-3 32e-3 | length of the domain (m)
`gas%%components` | N2 O2 | names of the gas components
`gas%%fractions` | 0.8 0.2 | gas fractions
`gas%%pressure` | 1.0 | pressure (bar)
`field_amplitude` | 2e6 | amplitude of the background field (V/m)

# Passing configuration files via the command line

One or more configuration files can be specified when running a simulation:

    ./streamer file_1.cfg file_2.cfg ...

Options from later files will override those from earlier files.

# Setting variables from the command line

Individual options can be specified via the command line:

    ./streamer file_1.cfg -var=value

If necessary, you can use quotes to specify for example an array:

    ./streamer file_1.cfg -var="a b c"

# A list of all options

To obtain a list with all simulation options, first run a simulation. If the
value for `output%name` is for example `output/streamer_2d`, then a file
`output/streamer_2d_out.cfg` should be generated with all options and their
values.
