#!/bin/bash

cd ../

./streamer 3d_0.20C4F7N_0.80CO2.cfg -field_given_by="field 6.0e6" -background_density=1e11 -output%name=/export/scratch3/baohong/eco_gas/0.20C4F7N_0.80CO2/3d_0.20C4F7N_0.80CO2_neg_bd1e11_6.0e6

./streamer 3d_0.20C4F7N_0.80CO2.cfg -field_given_by="field 7.0e6" -background_density=1e11 -output%name=/export/scratch2/md/baohong/eco_gas/0.20C4F7N_0.80CO2/3d_0.20C4F7N_0.80CO2_neg_bd1e11_7.0e6