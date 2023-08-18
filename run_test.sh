#!/bin/bash

# Usage:
# ./run_test.sh (runs all tests)
# ./run_test.sh path_to/my_test.cfg (run one or more tests)
#
# The parent directory of a cfg file should contain the "streamer" executable

set -eu                         # Exit on errors
top_dir=$(dirname "$0")         # Top directory

run_test() {
    # Arguments: config file, top directory
    dir=$(dirname "$1")
    cfg=$(basename "$1")

    # Execute the test in the directory of the .cfg file
    start=`date +%s.%N`
    (cd "$dir" && ../streamer "$cfg" > run.log ||
             { cat run.log; return 1; })
    end=`date +%s.%N`
    runtime_sec=$(printf "%.2f" $(echo "$end - $start" | bc -l))

    if (($? != 0)); then
        echo "FAILED $1"
        return
    fi

    log_name="${cfg/.cfg/_rtest.log}"
    log_a="$dir/$log_name"
    log_b="$dir/output/$log_name"

    # Compare log files
    if "$2"/tools/compare_logs.py "$log_a" "$log_b"; then
        echo "PASSED $1 (${runtime_sec})"
    else
        echo "FAILED $1 (${runtime_sec})"
    fi
}

# Use array to store test results
declare -a test_results

if (( "$#" > 0 )); then
    # Run given tests (note that this does not compile the code)
    for cfg in "$@"; do
        out=$(run_test "$cfg" "$top_dir")
        echo "$out"
        results+=("$out")
    done
else
    # Run all tests
    declare -a test_dirs=(
        "programs/standard_1d/tests"
        "programs/standard_2d/tests"
        "programs/standard_3d/tests"
        "programs/dielectric_2d/tests"
    )

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
