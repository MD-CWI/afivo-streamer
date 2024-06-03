# Installation, compilation and requirements

[TOC]

# Requirements {#requirements}

* A UNIX-like operating system such as GNU/Linux, or the Windows subsystem for linux (see @ref wsl-installation)
* A recent Fortran compiler such as `gfortran`, which supports Fortran 2008
* A recent C compiler such as `gcc` to compile the included libraries
* `git` to download and update the source code
* `python3` with `numpy`, `pandas` and `matplotlib` for some of the included tools

# Obtaining the source code {#obtaining-src}

Run the following command in a folder to clone the complete repository:

    git clone https://github.com/MD-CWI/afivo-streamer.git

If you only want the `master` branch, then you can use

    git clone -b master --single-branch https://github.com/MD-CWI/afivo-streamer

The code can also be downloaded as a [zip file](https://github.com/MD-CWI/afivo-streamer/archive/refs/heads/master.zip).

# Compiling the source code {#compiling}

Go into the afivo-streamer folder and compile the code:

    cd afivo-streamer
    make

This will first compile the `afivo` library and its dependencies, and afterwards compile most of the programs included with `afivo-streamer`. To test whether the code is working as expected, you can run several tests using

    bash run_tests.sh

To run some of the examples, see [running examples](documentation/running_examples.md).

# Updating to the latest version {#updating}

If you want to update your previously downloaded code, go into your afivo-streamer folder and pull the new version:

    git pull

Afterwards, you can recompile the code by typing

    make

either in the `afivo-streamer` folder, or in a specific program folder with a `Makefile`.

# Compilation flags and recompilation {#compilation-flags}

A number of flags can be set to help with debugging or performance testing. To apply these, first remove all previously compiled files:

    make allclean

This is also useful when upgrading e.g. the compiler on your system. Afterwards, a number of flags can be applied, most importantly:

* `make DEBUG=1` leads to a much slower executable that does all kinds of error checking (array bounds, division by zero etc.)
* `make PROF=gprof` compile with `-pg` support for `gprof`
* `make PROF=gperftools` link with `libprofiler` to enable `gperftools`, see https://teunissen.net/afivo/md_documentation_profiling.html
* `make COMPILER=ifort` to use the `ifort` compiler

Afterwards, perform `make allclean` again to revert to the standard compilation settings.

# Windows subsystem for linux {#wsl-installation}

Follow the instructions on https://learn.microsoft.com/en-us/windows/wsl/install to install `wsl` (windows subsystem for linux). After installation, launch `wsl`. For Ubuntu/Debian (Ubuntu is the default), the required compilers and python tools can be installed with:

    sudo apt update
    sudo apt install gfortran gcc g++ make git python3 python3-numpy python3-pandas python3-matplotlib

Afterwards, the regular installation instructions can be followed. The Linux files are also accessible from Windows, so simulation output can be visualized from the Windows environment.

# Potential issues and solutions {#issues-solutions}

## Compilation errors after pulling a new version

1. Remove the old compiles with `make allclean`
2. Rarely, the Silo or Hypre libraries have to be recompiled. To do this, go to `afivo/external_libraries` and execute the `build_...` scripts.

## Specifying the compilers to use

Sometimes, no suitable compilers are found (or used). This can be overcome by manually specifying which compilers to use. If `gcc`, `g++` and `gfortran` can be found in `/usr/bin/`, then this can  for example be done with:
    * `export CC=/usr/bin/gcc`
    * `export CXX=/usr/bin/g++`
    * `export FORT=/usr/bin/gfortran`

## Missing libraries

* `/usr/bin/ld: cannot find -lsz`. **Solution** (on Fedora) `sudo dnf install libaec-devel`
