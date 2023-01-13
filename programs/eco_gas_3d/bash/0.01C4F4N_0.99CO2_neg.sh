#!/bin/bash

cd ../

./streamer 3d_0.01C4F7N_0.99CO2.cfg -field_given_by="field 3.0e6" -background_density=0e9 -output%name=/export/scratch2/baohong/eco_gas/0.01C4F7N_0.99CO2/3d_d5mm_neg_bd0_3.0e6

./streamer 3d_0.01C4F7N_0.99CO2.cfg -field_given_by="field 3.2e6" -background_density=0e9 -output%name=/export/scratch2/baohong/eco_gas/0.01C4F7N_0.99CO2/3d_d5mm_neg_bd0_3.2e6

./streamer 3d_0.01C4F7N_0.99CO2.cfg -field_given_by="field 3.4e6" -background_density=0e9 -output%name=/export/scratch2/baohong/eco_gas/0.01C4F7N_0.99CO2/3d_d5mm_neg_bd0_3.4e6
