#!/bin/bash

cd ../

#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_0.03H2O_1.0e6
#./streamer d200mm_0.03H2O_photoi_custom.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_0.03H2O_1.0e6_photoifit-0.03H2O

#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_0.03H2O_1.4e6
#./streamer d200mm_0.03H2O_photoi_custom.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_0.03H2O_1.4e6_photoifit-0.03H2O

./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.9e6" -field_rod_r0="0 0.9525" -field_rod_radius=0.5e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_0.9e6_EL10mm_0.5mm
./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.95e6" -field_rod_r0="0 0.9525" -field_rod_radius=0.5e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_0.95e6_EL10mm_0.5mm
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -field_rod_r0="0 0.9525" -field_rod_radius=0.5e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_EL10mm_0.5mm