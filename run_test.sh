#!/bin/bash

# Usage:
# ./run_test.sh (runs all tests)
# ./run_test.sh my_test.cfg (run a single test)
#
# The parent directory of a cfg file should contain the "streamer" executable

set -eu                         # Exit on errors
top_dir=$(dirname "$0")         # Top directory

run_test() {
    # Arguments: config file, top directory
    # Run simulation (executable should be in parent folder)
    echo "Running $1"

    dir=$(dirname "$1")
    cfg=$(basename "$1")

    # Execute the test in the directory of the .cfg file
    (cd "$dir" && ../streamer "$cfg" > run.log|| { echo "FAILED $1"; return; })

    # Original output should have this name
    log_a="${1/.cfg/_rtest_orig.log}"
    # A file test_X.cfg should produce a log file output/test_X_rtest.log
    log_b="$dir/output/${cfg/.cfg/_rtest.log}"

    # Compare log files
    if "$2"/tools/compare_logs.py "$log_a" "$log_b" ; then
        echo "PASSED $1"
    else
        echo "FAILED $1"
    fi
}

if (( "$#" == 1 )); then
    # Run a single test
    run_test "$1" "$top_dir"
else
    # Run all tests
    declare -a test_dirs=("programs/standard_2d/tests")

    for dir in "${test_dirs[@]}"; do
        # Compile in parent folder and make sure 'output' exists
        (cd "$dir" && make -j -C .. --silent && mkdir -p output)

        # Loop over the .cfg files
        for cfg in "$dir"/*.cfg; do
            run_test "$cfg" "$top_dir"
        done
    done
fi
