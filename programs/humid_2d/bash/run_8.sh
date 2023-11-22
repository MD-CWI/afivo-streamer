#!/bin/bash

cd ../

./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.8e6" -field_rod_r0="0 0.95725" -field_rod_radius=0.45e-3 -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch3/baohong/humid_air/comparison/d200_starikovskiy_0.03H2O_0.8e6_EL9mm_0.45mm
./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.8e6" -field_rod_r0="0 0.9525" -field_rod_radius=0.5e-3 -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch3/baohong/humid_air/comparison/d200_starikovskiy_0.03H2O_0.8e6_EL10mm_0.5mm