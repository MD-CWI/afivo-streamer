#!/bin/bash

export OMP_NUM_THREADS=1

parallel --jobs 4 --ungroup < commands.txt
