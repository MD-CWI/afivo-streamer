# Source code structure

# Introduction

The source code of Afivo is located in different folders:

Folder | Description
---|---
src | Contains the source code of the framework
examples | Contains example programs
documentation | Contains the Doxygen documentation
tests | Contains some simple tests of the framework
`external_libraries` | After compilation, contains the Silo library
`misc_files` | Miscellaneous and old files (could be removed)

# The src folder

The folder `src` contains the source of Afivo, split over several modules. The
`m_aX_[...].f90` files contain code for both the 2D and 3D version. These files
are automatically translated to 2D and 3D version (`m_a2_[...].f90` and
`m_a3_[...].f90`), as described below. **Modifications should therefore be made
in the `m_aX_[...].f90` files**.

## List of modules and prefixes

The table below contains a brief description of Afivo's modules. The `a2_`
modules have 3D equivalents with the prefix `a3_`.

Module name | Prefix | Description
---|---|---
`m_a2_core` | `a2_` | Core routines to construct and refine the mesh
`m_a2_types` | `a2_` | Dimension-dependent types and routines related to indexing
`m_a2_utils` | `a2_` | All kinds of different 'helper' routines for Afivo
`m_a2_output` | `a2_write_` | Routines for writing output
`m_a2_ghostcell` | `a2_gc_` | Routines for filling ghost cells
`m_a2_restrict` | `a2_restrict_` | Routines for restriction (going from fine to coarse values)
`m_a2_prolong` | `a2_prolong_` | Routines for prolongation (going from coarse to fine values)
`m_a2_interp` | `a2_interp_` | Routines for interpolation on the grid
`m_a2_multigrid` | `mg2_` | Multigrid routines and data structures
`m_a2_all` | - | Module that contains all of Afivo

# Translating code to 2D/3D

The translation is performed using a couple simple preprocessor rules:

Pattern | Replacement 2D | Replacement 3D
---|---|---
`$D` | 2 | 3
`IJK` | `i, j` | `i, j, k`
`do KJI_DO(a,b)` | `do j=a,b; do i=a,b` | `do k=a,b; do j=a,b; do i=a,b`

To simplify writing loops, the following patterns are defined in
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

Furthermore, blocks specific for 2D and 3D can be written as:

    #if $D == 2
        ...
    #elif $D == 3
        ...
    #endif

