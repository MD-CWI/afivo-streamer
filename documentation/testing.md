# Testing

# Introduction

Currently, there are so-called [regression tests](https://en.wikipedia.org/wiki/Regression_testing). The idea is to run a short simulation, and compare its output to stored results. This comparison is done with a certain tolerance, because results can be slightly different, for example when using a different CPU.

# Running tests

To run all tests, execute the following command in the `afivo-streamer` directory:

    ./run_test.sh

This should produce output like

    PASSED programs/standard_2d/tests/test_2d.cfg
    PASSED programs/standard_2d/tests/test_2d_photoi.cfg
    PASSED programs/standard_2d/tests/test_cyl.cfg
    ...

It is also possible to run individual tests, for example:

    ./run_test.sh programs/standard_2d/tests/test_2d.cfg

# Adding new tests

To add a new test called `my_test`:

1. Add a file `my_test.cfg` in `programs/.../tests`, so for example in `programs/standard_2d/tests`. The variable `output%name` should be set to `output/my_test`.
2. Run the simulation from the `tests` directory with `../streamer my_test.cfg`
3. Copy the log file `tests/output/my_test_rtest.log` to `tests/`

You should now be able to run the test, for example like this:

    # In the directory of my_test.cfg:
    <afivo-streamer-path>/run_tesh.sh my_test.cfg

Optionally, add the `tests/` directory to the variable `test_dirs` in `run_test.sh` if it is not yet present. The test will then run automatically, so make sure it doesn't take too long.

# Test configuration

The test configuration can be found in the file `.gitlab-ci.yml`.
