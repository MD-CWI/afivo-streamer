#!/bin/bash

cd ../

./streamer 3d_0.01C4F7N_0.99CO2.cfg -field_given_by="field -2.8e6" -background_density=1e9 -output%name=/export/scratch2/baohong/eco_gas/0.01C4F7N_0.99CO2/3d_d5mm_pos_bd1e9_2.8e6

./streamer 3d_0.01C4F7N_0.99CO2.cfg -field_given_by="field -3.0e6" -background_density=1e9 -output%name=/export/scratch2/baohong/eco_gas/0.01C4F7N_0.99CO2/3d_d5mm_pos_bd1e9_3.0e6

./streamer 3d_0.01C4F7N_0.99CO2.cfg -field_given_by="field -3.2e6" -background_density=1e9 -output%name=/export/scratch2/baohong/eco_gas/0.01C4F7N_0.99CO2/3d_d5mm_pos_bd1e9_3.2e6



