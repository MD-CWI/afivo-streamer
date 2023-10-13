# Afivo-streamer

Afivo-streamer is a tool for the simulation of streamer discharges in 1D, 2D and 3D,
using plasma fluid models. Afivo-streamer makes uses of the [afivo
framework](https://github.com/MD-CWI/afivo), which for example provides
adaptive mesh refinement and parallelization.

## Quick links

* [git repository](https://github.com/MD-CWI/afivo-streamer)
* [Installation instructions](documentation/installation.md)
* [Contents of the documentation](@ref doc-contents)
* [Recent changes](https://github.com/MD-CWI/afivo-streamer/activity)

## Opiniated comparison against typical alternatives (advantages and disadvantages)

### Strong points

* **Computational efficiency**: relative to typical other codes, afivo-streamer has a high computational efficiency, see e.g. [this paper](https://doi.org/10.1088/1361-6595/aad768). This usually makes it possible to do 2D simulations on a normal computer rather interactively.
* **Possibility of doing 3D simulations**: due to the computational efficiency it is possible to do relevant 3D simulations in a reasonable time on typical hardware
* **Usage**: the code has frequently been used for streamer simulations in various publications, which also means that things like plasma chemisty and different photoionization models are included
* **Open source**: anyone can in principle change the code (although doing so might not always be easy)

### Neutral points

* **Linux/unix**: basic familiarity with such systems (and the command line) is required to compile and run the code
* **Fortran**: experience with Fortran is required in order to change the code

### Weak points

* **Complexity**: it can be rather difficult to change the code or to add new physics, in particular due to the use of adaptive mesh refinement and parallelization
* **Research type code**: which means there is not much documentation, no extensive test suite, there are experimental features, etc.
* **Complex geometries**: although curved electrodes can be included, curved dielectrics and complex geometries are generally not supported





