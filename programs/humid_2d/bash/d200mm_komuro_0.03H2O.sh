#!/bin/bash

cd ../

./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.8e6" -input_data%file=../../transport_data/humid_air/komuro_humid0.03_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_0.03H2O_0.8e6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.2e6" -input_data%file=../../transport_data/humid_air/komuro_humid0.03_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_0.03H2O_1.2e6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.6e6" -input_data%file=../../transport_data/humid_air/komuro_humid0.03_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_0.03H2O_1.6e6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -2.0e6" -input_data%file=../../transport_data/humid_air/komuro_humid0.03_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_0.03H2O_2.0e6