# Writing output

# Silo and VTK output

TODO

# Binary output

Afivo-streamer allows a simulation to be continued at a later time through the use of binary output files (DAT files). To do this, the simulation is initially run and set up to write DAT files, which could be achieved by using the following `.cfg` file options:

 * `datfile%%write` - set to true if you want to write DAT files
 * `datfile%%per_outputs` - control when DAT files are written by equating this to an integer (e.g. datfile%%per_outputs = 20 means that a DAT file is created every after 20 output files are written)

The created DAT files would have the same name as the SILO files written by the program. To begin another simulation run from a particular DAT file, set the flag restart_from_file equal to a string with the name (and location, if in a different directory) of the DAT file you want to use as the starting point of the other simulation run.

# The log file

The log file is a text file containing information about the physics and numerical properties of the simulation. It is possible to write extra variables to the log file by defining the routine `user_log_variables()`, see `m_user_methods`.

Which variables the log file contains depends on the dimensionality of the run. The meaning of these variables is described below

name | meaning
---|---
`it` | simulation iteration
`time` | simulation time
`dt` | time step
`v` | estimate of streamer velocity (compared to last output)
`sum(n_e)` | integral of electron density
`sum(n_i)` | integral of first positive ion species
`sum(charge)` | integral of charge density (considering all species)
`max(E) x y` | maximum electric field + location
`max(n_e) x y` | maximum electron density + location
`max(E_r) x y` | maximum radial field + location
`min(E_r)` | minimum radial field
`voltage` | current applied voltage
`wc_time` | simulation wall-clock time (how long it has been running)
`n_cells` | number of grid cells used in simulation
`min(dx)` | minimum grid spacing
`highest(lvl)` | highest refinement level

# 2D planes

TODO

# 1D lines

TODO
