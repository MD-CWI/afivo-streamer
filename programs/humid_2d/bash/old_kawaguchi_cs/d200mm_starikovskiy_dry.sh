#!/bin/bash

cd ../

./streamer d200mm_dry_air.cfg -field_given_by="field -0.8e6" -input_data%file=../../transport_data/humid_air/starikovskiy_dry_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_dry_0.8e6
./streamer d200mm_dry_air.cfg -field_given_by="field -1.2e6" -input_data%file=../../transport_data/humid_air/starikovskiy_dry_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_dry_1.2e6
./streamer d200mm_dry_air.cfg -field_given_by="field -1.6e6" -input_data%file=../../transport_data/humid_air/starikovskiy_dry_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_dry_1.6e6
./streamer d200mm_dry_air.cfg -field_given_by="field -2.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_dry_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_dry_2.0e6

