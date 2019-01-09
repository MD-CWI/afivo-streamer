# Afivo-streamer

This is a collection of streamer models based on
[Afivo](https://github.com/jannisteunissen/afivo). A 2D, a 3D and a cylindrical
model are included.

## News/updates

* In November/December 2018 I (Jannis) am planning to improve this code a bit so
  that it becomes easier to write custom user applications.

## How to use

### Step 1: Get the code

Using `git`:

    git clone https://gitlab.com/MD-CWI-NL/afivo-streamer.git

Alternatively, a [zip of the latest
version](https://gitlab.com/MD-CWI-NL/afivo-streamer/repository/archive.zip?ref=master)
can be downloaded, but this is not recommended, since there is no update
mechanism.

### Step 2: Compile

    cd afivo_streamer
    make

### Step 3: Run the examples

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

### Requirements

* A unix-like system (e.g., GNU/Linux)
* Gfortran 4.8 or newer
* git (recommended)

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

* Paper about Afivo-streamer: [Simulating streamer discharges in 3D with the
  parallel adaptive Afivo framework](https://doi.org/10.1088/1361-6463/aa8faf)
* First paper using Afivo-streamer: [The role of free electrons in the guiding
  of positive streamers](http://dx.doi.org/10.1088/0963-0252/25/4/044001)
* Paper about Afivo: [Afivo: a framework for quadtree/octree AMR with
  shared-memory parallelization and geometric multigrid
  methods](https://doi.org/10.1016/j.cpc.2018.06.018)
  [arXiv](https://arxiv.org/abs/1701.04329)
* More information about the photoionization method can be found in chapter 11
  of: [3D Simulations and Analysis of Pulsed
  Discharges](https://research.tue.nl/en/publications/3d-simulations-and-analysis-of-pulsed-discharges),
  J. Teunissen, PhD Thesis, Eindhoven University of Technology, 2015.
* [Comparison of six simulation codes for positive streamers in air](https://doi.org/10.1088/1361-6595/aad768)
