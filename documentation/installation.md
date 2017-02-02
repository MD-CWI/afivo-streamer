# Installation instructions

# Getting the source code

If you have `git` installed, you can download a copy of the repository with:

    $ git clone https://gitlab.com/MD-CWI-NL/afivo.git

To update the local copy, use the command

    $ git pull

in the `afivo` directory. You can also download
a
[zip file with the latest version](https://gitlab.com/MD-CWI-NL/afivo/repository/archive.zip?ref=master),
but then it harder to update.

# Compiling

On an up-to-date system with `gfortran`, simply do:

    $ cd afivo
    $ make -j 4

With `make -j 4` the build is done in parallel, using four processors. To build
sequentially, remove the `-j` argument. If you want to compile with runtime
checks (useful for debugging) and profiling options, then do:

    $ make clean
    $ make DEBUG=1

It is also possible to compile with the Intel Fortran compiler, using the
command

    $ make COMPILER=ifort

# Requirements

* Fortran 2008 compatible compiler (gfortran >= 4.8 has been tested)
* C/C++ compiler for compiling Silo (see below).

Afivo has one dependency: Silo. The script `external_libraries/build_silo.sh`
automatically downloads, configures and compiles Silo into your Afivo directory.
If this fails, then you should probably download, configure and install Silo
yourself. If you really do not want to use Silo (Afivo also supports .vtu
output), you could remove the Silo routines.
