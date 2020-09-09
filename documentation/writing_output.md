# Writing output

# Silo and VTK output

TODO

# Binary output

Afivo-streamer allows a simulation to be continued at a later time through the use of binary output files (DAT files). To do this, the simulation is initially run and set up to write DAT files, which could be achieved by using the following `.cfg` file options:

 * `datfile%%write` - set to true if you want to write DAT files
 * `datfile%%per_outputs` - control when DAT files are written by equating this to an integer (e.g. datfile%%per_outputs = 20 means that a DAT file is created every after 20 output files are written)

The created DAT files would have the same name as the SILO files written by the program. To begin another simulation run from a particular DAT file, set the flag restart_from_file equal to a string with the name (and location, if in a different directory) of the DAT file you want to use as the starting point of the other simulation run.

# The log file

TODO

# 2D planes

TODO

# 1D lines

TODO
