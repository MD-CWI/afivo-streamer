#!/bin/bash

export PATH="$PATH:/opt/sw/intel/bin"
export LD_LIBRARY_PATH="$PATH:/opt/sw/intel/mkl/lib/intel64:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="$PATH:/opt/sw/intel/lib/intel64:$LD_LIBRARY_PATH"

./bolsigminus.exe "$@"


