#!/bin/bash

cd ../

./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.9e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Kawaguchi.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_0.9e6_0.10H2O_Bourdon
./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.9e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Itikawa.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_0.9e6_0.10H2O_Bourdon_Itikawa
./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.9e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Itikawa_v2.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_0.9e6_0.10H2O_Bourdon_Itikawa_v2