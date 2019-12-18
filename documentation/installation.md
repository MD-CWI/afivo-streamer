# Installation

# Requirements

* A UNIX-like operating system such as GNU/Linux
* A recent Fortran compiler such as gfortran 4.8 or newer
* A recent C compiler such as gcc to compile the [Silo](https://wci.llnl.gov/simulation/computer-codes/silo) library
* `git` to download and update the source code

# Installation

Run the following command in a folder to clone the repository:

    git clone https://gitlab.com/MD-CWI-NL/afivo-streamer.git

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
