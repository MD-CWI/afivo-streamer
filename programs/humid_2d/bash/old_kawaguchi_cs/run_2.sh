#!/bin/bash

cd ../

#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Itikawa.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_Itikawa
#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Itikawa_v2.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_Itikawa-v2

#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Itikawa.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.4e6_Itikawa
#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Itikawa_v2.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.4e6_Itikawa-v2

#./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Itikawa_v2.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_Itikawa-v2
#./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Morgan.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_Morgan

./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.8e6" -field_rod_r0="0 0.96675" -field_rod_radius=0.35e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_0.8e6_EL7mm_0.35mm
./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.85e6" -field_rod_r0="0 0.96675" -field_rod_radius=0.35e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_0.85e6_EL7mm_0.35mm
./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.9e6" -field_rod_r0="0 0.96675" -field_rod_radius=0.35e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_0.9e6_EL7mm_0.35mm