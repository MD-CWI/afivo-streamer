# Afivo-streamer

This is a collection of streamer models based on
[Afivo](https://github.com/jannisteunissen/afivo). A 2D, a 3D and a cylindrical
model are included.

## How to use

### Step 1: Get the code

Using `git`:

    git clone https://gitlab.com/MD-CWI-NL/afivo-streamer.git

Alternatively,
a
[zip of the latest version](https://gitlab.com/MD-CWI-NL/afivo-streamer/repository/archive.zip?ref=master) can
be downloaded, but this is not recommended, since there is no update mechanism.

### Step 2: Compile

    cd afivo_streamer
    make

### Step 3: Run the examples

    ./streamer_2d  configs/example_2d.txt
    ./streamer_cyl configs/example_cyl.txt
    ./streamer_3d  configs/example_3d.txt

where the configuration files include the parameters that you want to use, see
the examples in the `configs` directory. You can also specify multiple
configuration files, like

    ./streamer_2d cfg_1.txt cfg_2.txt ...

Options from later files will override those from earlier files.

Output can be written in
the [Silo](https://wci.llnl.gov/simulation/computer-codes/silo) or VTK
unstructured
format. [Visit](https://wci.llnl.gov/simulation/computer-codes/visit/downloads)
is the recommended tool to visualize the output.

### Requirements

* A unix-like system, e.g. `GNU/Linux`
* Gfortran 4.8 or newer

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

## References / how to cite

A paper specifically showcasing `afivo-streamer` is in the works. The first
paper in which the code was used is:

[The role of free electrons in the guiding of positive streamers](http://dx.doi.org/10.1088/0963-0252/25/4/044001),
S. Nijdam, J. Teunissen, E. Takahashi, U. Ebert, Plasma Sources Sci. Technol.
(2016)

More information about the photoionization method can be found in chapter 11 of:

[3D Simulations and Analysis of Pulsed Discharges](http://repository.tue.nl/801516),
J. Teunissen, PhD Thesis, Eindhoven University of Technology, 2015.

