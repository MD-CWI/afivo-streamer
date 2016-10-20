# Description of source code

This folder (src) contains the source of Afivo, split over several modules. The
code can be recompiled by typing

    $ make

The Silo library has to be available, which means that you need to type 'make'
at least once in the Afivo root folder (one directory up):

    $ cd ..
    $ make

The m_aX_[...].f90 files contain code for both the 2D and 3D version. The
corresponding m_a2_[...].f90 and m_a3_[...].f90 files are generated
automatically, so you should **not** edit those! Instead, modify the
m_aX_[...].f90 files!
