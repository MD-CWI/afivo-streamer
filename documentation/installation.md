# Installation instructions

# Getting the source code

If you have `git` installed, you can download a copy of the repository with:

    $ git clone https://github.com/MD-CWI/afivo.git

To update the local copy, use the command

    $ git pull

in the `afivo` directory. You can also download
a
[zip file with the latest version](https://github.com/MD-CWI/afivo/repository/archive.zip?ref=master),
but then it harder to update.

# Requirements

* A Linux/Unix system
* A recent Fortran compiler (gfortran >= 5 is recommended)
* C/C++ compiler for compiling Silo (gcc is recommended, also see below).

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

# Dependencies

Afivo has two dependencies: Silo and Hypre. The scripts in `external_libraries/`
automatically download, configure and compile these libraries into your Afivo directory.
If this fails, some known issues and solutions can be found at https://teunissen.net/afivo_streamer/md_documentation_installation.html