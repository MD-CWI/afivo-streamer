#!/bin/bash

plot_script=${HOME}/git/afivo/tools/plot_raw_data.py
raw_converter=${HOME}/git/afivo/tools/silo_to_raw

npixels=128
mkdir -p npz_${npixels}

# Note that the script converts all the silo files in the output folder
parallel -j 4 --eta $plot_script {2} -min_pixels $npixels -variable {1} \
         -silo_to_raw ${raw_converter} -save_npz npz_${npixels}/{2/.}_{1}.npz \
         ::: e phi electric_fld N2_B3 N2_C3 rhs ::: output/*.silo
