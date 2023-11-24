#!/bin/bash

cd ../

#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.03H2O_1.0e6
#./streamer d200mm_0.03H2O_photoi_custom.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.03H2O_1.0e6_photoifit-0.03H2O

#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.03H2O_1.4e6
#./streamer d200mm_0.03H2O_photoi_custom.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.03H2O_1.4e6_photoifit-0.03H2O

#./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.4e6_photoifit_0.10H2O
#./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.4e6_bourdon3_0.00H2O -photoi_helmh%author="Bourdon-3"

./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.8e6" -field_rod_r0="0 0.96675" -field_rod_radius=0.35e-3 -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.03H2O_0.8e6_EL7mm_0.35mm
./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.8e6" -field_rod_r0="0 0.962" -field_rod_radius=0.4e-3 -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.03H2O_0.8e6_EL8mm_0.4mm