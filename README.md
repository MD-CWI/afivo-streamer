# Afivo-streamer

This is a collection of streamer models based on
[Afivo](https://github.com/jannisteunissen/afivo). A 2D, a 3D and a cylindrical
model are included.

## Requirements

Gfortran 4.8 or newer. Afivo makes use of the
[Silo](https://wci.llnl.gov/simulation/computer-codes/silo/downloads) library,
which is automatically downloaded and compiled. It also supports writing VTK
unstructured grids in XML-binary format, for which no external libraries are
required.

## Getting the code

    git clone TODO [new url]

## Compiling the code

    cd afivo_streamer
    make

## Running the code

    ./streamer_2d  ConfigureFiles/example_2d.txt
    ./streamer_cyl ConfigureFiles/example_cyl.txt
    ./streamer_3d  ConfigureFiles/example_3d.txt

where the configuration files include the parameters that you want to use, see
the examples in the `configs` directory. You can also specify multiple
configuration files, like

    ./streamer_2d cfg_1.txt cfg_2.txt ...

Options from later files will override those from earlier files.

## Getting input data (transport coefficients)

Transport data can be computed from cross sections with for example:

* [Bolos](https://github.com/aluque/bolos) (Open source two-term Boltzmann
  solver)
* [Bolsig+](http://www.bolsig.laplace.univ-tlse.fr) (Popular two-term Boltzmann
  solver)
* [Magboltz](http://consult.cern.ch/writeup/magboltz/)
* [Particle_swarm](https://gitlab.com/MD-CWI-NL/particle_swarm) (Particle /
  Monte Carlo transport data computation)

Cross section can be obtained from http://lxcat.net. The format of the transport
data used for the streamer models is illustrated in the files in the `transport_data`
folder.

## Format of output data

Output can be written in the
[Silo](https://wci.llnl.gov/simulation/computer-codes/silo) or VTK unstructured
format. [Visit](https://wci.llnl.gov/simulation/computer-codes/visit/downloads)
is the recommended tool to visualize the output.
