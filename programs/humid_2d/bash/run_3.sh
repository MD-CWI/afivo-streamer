#!/bin/bash

cd ../

#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Morgan.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_Morgan
#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Phelps.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_Phelps

#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Morgan.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.4e6_Morgan
#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Phelps.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.4e6_Phelps

#./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Phelps.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_Phelps
#./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Triniti.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_Triniti

./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.95e6" -field_rod_r0="0 0.96675" -field_rod_radius=0.35e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_0.95e6_EL7mm_0.35mm
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -field_rod_r0="0 0.96675" -field_rod_radius=0.35e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_EL7mm_0.35mm

./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.8e6" -field_rod_r0="0 0.95725" -field_rod_radius=0.45e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_0.8e6_EL9mm_0.45mm