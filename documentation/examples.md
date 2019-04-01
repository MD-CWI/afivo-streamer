# Examples

    cd programs/standard_2d
    ./streamer streamer_2d.cfg
    ./streamer streamer_cyl.cfg

    cd programs/standard_3d
    ./streamer_3d configs/streamer_3d.cfg

where the configuration files include the parameters that you want to use, see
the examples in the `configs` directory. You can also specify multiple
configuration files, like

    ./streamer cfg_1.txt cfg_2.txt ...

Options from later files will override those from earlier files.

Output can be written in the
[Silo](https://wci.llnl.gov/simulation/computer-codes/silo) or VTK unstructured
format. [Visit](https://wci.llnl.gov/simulation/computer-codes/visit/downloads)
is the recommended tool to visualize the output.

