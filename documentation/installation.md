# Installation

# Requirements

* A UNIX-like operating system such as GNU/Linux
* A recent Fortran compiler such as gfortran 4.8 or newer
* A recent C compiler such as gcc to compile the [Silo](https://wci.llnl.gov/simulation/computer-codes/silo) library
* `git` to download and update the source code

# Installation

Run the following command in a folder to clone the repository:

    git clone https://github.com/MD-CWI/afivo-streamer.git

Then you can go into the folder and compile the code:

    cd afivo-streamer
    make

Afterwards, you can run some of the example, see @ref md_documentation_examples.

# Updating to the latest version

If you want to update your previously downloaded code, go into your afivo-streamer folder and pull the new version:

    git pull

Afterwards, you can recompile the code by typing

    make

either in the `afivo-streamer` folder, or in a specific program folder with a `Makefile`.

# Compilation flags and recompilation

A number of flags can be set to help with debugging or performance testing. To apply these, first remove all previously compiled files:

    make allclean

This is also useful when upgrading e.g. the compiler on your system. Afterwards, a number of flags can be applied, most importantly:

* `make DEBUG=1` leads to a much slower executable that does all kinds of error checking (array bounds, division by zero etc.)
* `make PROF=gprof` compile with `-pg` support for `gprof`
* `make PROF=gperftools` link with `libprofiler` to enable `gperftools`, see https://teunissen.net/afivo/md_documentation_profiling.html
* `make COMPILER=ifort` to use the `ifort` compiler

Afterwards, perform `make allclean` again to revert to the standard compilation settings.

# List of issues and solutions

## Installation 'Error 2':
* Simply typing `make` to compile the code in Fedora systems (30 and above) gives an error. A closer inspection indicates that the `g77` compiler is missing. This is because the makefile is unable to find the `gcc`, `gfortran` and `g++` compilers. This can be overcome by typing the following in the terminal:
    * `export CC=/usr/bin/gcc`
    * `export CXX=/usr/bin/g++`
    * `export FORT=/usr/bin/gfortran`
* *Attention*: The above locations `/usr/bin` is the default installation location for the GNU compilers. If you have them installed in a different directory (in case you have multiple versions of compilers), then make sure you use that particular location (instead of `/usr/bin`).

## Problems compiling Silo

* The script `build_silo.sh` can be executed manually from the `afivo/external_libraries' folder for easier debugging. 
* Make sure the default C, C++ and Fortran compiler are configured correctly. For example, use `export FC=/usr/bin/gfortran' in your shell to specify the Fortran compiler.
* `/usr/bin/ld: cannot find -lsz`. **Solution** (on Fedora) `sudo dnf install libaec-devel`

### The above solution gives a new error (encountered by some users using Fedora 30)

* Add the following lines to `build_silo.sh` just before the `#Configure` comment:
    `export CC=/usr/bin/gcc`
    `export CXX=/usr/bin/g++`
    `export FC=/usr/bin/gfortran`
    
