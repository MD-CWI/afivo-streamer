# Tutorial: running built-in examples

To run the standard 2D examples, go to the folder `programs/standard_2d`:

    cd programs/standard_2d

To be sure that the simulation code is compiled, you can compile it with:

    make

An executable `streamer` should be present. You can run it with a configuration
file like this:

    ./streamer streamer_2d.cfg

This example writes its output to a subdirectory `output`. If this directory does not exist, the code will abort with a message

    Output name: output/streamer_2d_...
    ERROR STOP Directory not writable (does it exist?)

After creating the output directory with `mkdir output`, the code should run and print (among other things) the following information:

    af_write_silo: written output/streamer_2d_000001.silo
    af_write_silo: written output/streamer_2d_000002.silo
    af_write_silo: written output/streamer_2d_000003.silo
    ...

These Silo files can be visualized using for example Visit, as explained at [Saving output and visualization](documentation/output_and_visualization.md). An example is included below.

![2D Pseudocolor plot of electric field](images/visit-pseudocolor-example.png){html: width=80%}

The 2D version of the `streamer` executable can also be used for axisymmetric simulations:

    ./streamer streamer_cyl.cfg

There are also 1D and 3D examples, which can be found in

    programs/standard_1d
    programs/standard_3d (note that this example takes a long time!)
