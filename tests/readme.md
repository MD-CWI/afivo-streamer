Afivo tests
==

Each file in this folder with the extension `.90` is considered a test, with the
exception of module files, which should have a prefixe `m_`. The output of the
tests is compared against a stored answer, located in the `answers` directory.

Running tests
===

All tests:

    make

Perform individual test

    make results/test_name

Adding tests
===

Simply add a new program to this directory

Editing/adding answers
===

Adding or updating an *answer* (i.e., the correct output for a test):

    make answers/test_name
