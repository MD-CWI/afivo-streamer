#!/bin/bash

cd ../

./streamer 3d_0.10C4F7N_0.90CO2.cfg -field_given_by="field 5.0e6" -background_density=1e11 -output%name=/export/scratch2/baohong/eco_gas/0.10C4F7N_0.90CO2/3d_0.10C4F7N_0.90CO2_neg_bd1e11_5.0e6

./streamer 3d_0.10C4F7N_0.90CO2.cfg -field_given_by="field 5.7e6" -background_density=1e11 -output%name=/export/scratch2/baohong/eco_gas/0.10C4F7N_0.90CO2/3d_0.10C4F7N_0.90CO2_neg_bd1e11_5.7e6