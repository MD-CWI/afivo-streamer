Afivo tests
==

Each file in this folder with the extension `.90` is considered a test, with the
exception of module files, which should have a prefixe `m_`. The output of the
tests is compared against a stored answer, located in the `answers` directory.

Running tests
===

All tests:

    make

which should (after compilation) produce output such as:

    PASSED test_refinement
    PASSED test_ghostcell
    PASSED test_init
    PASSED test_reduction_2d
    PASSED test_morton
    PASSED test_types_2d_3d

Perform an individual test

    make results/test_name

Adding tests
===

Simply add a new program to this directory.

Editing/adding answers
===

Adding or updating an *answer* (i.e., the correct output for a test):

    make answers/test_name
