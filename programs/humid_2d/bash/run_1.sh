#!/bin/bash

cd ../

#./streamer d200mm_dry_air.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_dry_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_dry_1.0e6
#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6

#./streamer d200mm_dry_air.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_dry_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_dry_1.4e6
#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.4e6

#./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_Kawaguchi
#./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Itikawa.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_Itikawa

./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.8e6" -field_rod_r0="0 0.962" -field_rod_radius=0.4e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_0.8e6_EL8mm_0.4mm
./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.85e6" -field_rod_r0="0 0.962" -field_rod_radius=0.4e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_0.85e6_EL8mm_0.4mm
./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.95e6" -field_rod_r0="0 0.962" -field_rod_radius=0.4e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_0.95e6_EL8mm_0.4mm