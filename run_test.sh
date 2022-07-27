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
    dir=$(dirname "$1")
    cfg=$(basename "$1")

    # Execute the test in the directory of the .cfg file
    (cd "$dir" && ../streamer "$cfg" > run.log ||
             { cat run.log; return 1; })

    if (($? != 0)); then
        echo "FAILED $1"
        return
    fi

    log_name="${cfg/.cfg/_rtest.log}"
    log_a="$dir/$log_name"
    log_b="$dir/output/$log_name"

    # Compare log files
    if "$2"/tools/compare_logs.py "$log_a" "$log_b"; then
        echo "PASSED $1"
    else
        echo "FAILED $1"
    fi
}

# Use array to store test results
declare -a test_results

if (( "$#" == 1 )); then
    # Run a single test
    out=$(run_test "$1" "$top_dir")
    echo "$out"
    results+=("$out")
else
    # Run all tests
    declare -a test_dirs=("programs/standard_1d/tests"
                          "programs/standard_2d/tests"
                          "programs/standard_3d/tests")

    for dir in "${test_dirs[@]}"; do
        # Compile in parent folder and make sure 'output' exists
        (cd "$dir" && make -j -C .. --silent && mkdir -p output)

        # Loop over the .cfg files
        for cfg in "$dir"/*.cfg; do
            out=$(run_test "$cfg" "$top_dir")
            echo "$out"
            results+=("$out")
        done
    done
fi

# Check if any of the results contained "FAILED"
if [[ "${results[@]}" =~ "FAILED" ]]; then
    exit 1
fi
