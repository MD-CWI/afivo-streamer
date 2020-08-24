# Examples

## Standard 2D examples

To run the standard 2D examples, go to the folder `programs/standard_2d`:

    cd programs/standard_2d

To be sure that the simulation is up to date, you can compile it with:

    make

An executable `streamer` should be present. You can run it with a configuration
file like this:

    ./streamer streamer_2d.cfg

The 2D program can also be used for axisymmetric simulations:

    ./streamer streamer_cyl.cfg

For more information about the configuration files, see @ref md_documentation_simulation_options.

## Examples in 1D and 3D

Examples in 1D and 3D can be found in

    programs/standard_1d
    programs/standard_3d

## Visualization

Output can be written in the
[Silo](https://wci.llnl.gov/simulation/computer-codes/silo) (recommended) or VTK
unstructured format.
[Visit](https://wci.llnl.gov/simulation/computer-codes/visit/downloads) is the
recommended tool to visualize the output.

