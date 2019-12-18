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

If you want to update your previously downloaded code, go into your afivo-streamer folder and pull the new version:

    cd afivo-streamer
    git clone https://gitlab.com/MD-CWI-NL/afivo-streamer.git

See @ref md_documentation_examples for a list of examples.
