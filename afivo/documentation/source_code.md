# Source code structure

# Introduction

The source code of Afivo is located in different folders:

Folder | Description
---|---
src | Contains the source code of the framework
examples | Contains example programs
documentation | Contains the Doxygen documentation
tests | Contains tests of the framework
`external_libraries` | Contains required external libraries
`lib_2d` and `lib_3d` | The Afivo library in 2D/3D

# The src folder

The folder `src` contains the source of Afivo, split over several modules (see
the **Modules** page). These modules contain the code for both 2D and 3D
functionality. The translation is performed using a couple simple preprocessor
rules:

Pattern | Replacement 2D | Replacement 3D
---|---|---
`NDIM` | 2 | 3
`IJK` | `i, j` | `i, j, k`

To simplify writing loops, a few patterns are defined in
`cpp_macros_2d.h` and `cpp_macros_3d.h`:

    do KJI_DO(a,b)
    --->
    2D: do j = a, b; do i = a, b
    3D: do k = a, b; do j = a, b; do i = a, b

    end do; CLOSE_DO
    --->
    2D: end do; end do
    3D: end do; end do; end do

    DTIMES(xx)
    --->
    2D: xx, xx
    3D: xx, xx, xx

Blocks of code specific for 2D and 3D can be written as:

    #if NDIM == 2
        ...
    #elif NDIM == 3
        ...
    #endif

