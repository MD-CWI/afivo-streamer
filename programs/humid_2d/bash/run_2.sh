#!/bin/bash

cd ../

./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.9e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Morgan.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_0.9e6_0.10H2O_Bourdon_Morgan
./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.9e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_0.9e6_0.10H2O_Bourdon_Phelps
./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.9e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Triniti.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_0.9e6_0.10H2O_Bourdon_Triniti