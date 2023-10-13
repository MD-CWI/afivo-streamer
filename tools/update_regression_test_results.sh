#!/bin/bash

# This script should be executed from the afivo-streamer top directory. It will
# update existing regression tests with the results generated in the respective
# "output" folders.
#
# Of course, this should only be done when you are confident that the resulting
# change is an improvement over the old behavior. Furthermore, any commit that
# changes results should have RESULTS_CHANGE in the title.

# Author: Jannis Teunissen

# Locate all existing regression tests
test_files=$(find . -iwholename "*tests/test_*rtest.log")

for regression_test in $test_files; do
    # Convert tests/test_name.log to tests/output/test_name.log
    new_result="${regression_test/tests\//tests\/output\/}"
    mv "$new_result" "$regression_test"
done

