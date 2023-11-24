#!/bin/bash

cd ../

./streamer d200mm_pure_H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_pure_H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_pure_H2O_1.0e6_Kawaguchi
./streamer d200mm_pure_H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_pure_H2O_Itikawa.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_pure_H2O_1.0e6_Itikawa
./streamer d200mm_pure_H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_pure_H2O_Itikawa_v2.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_pure_H2O_1.0e6_Itikawa-v2
./streamer d200mm_pure_H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_pure_H2O_Morgan.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_pure_H2O_1.0e6_Morgan
./streamer d200mm_pure_H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_pure_H2O_Phelps.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_pure_H2O_1.0e6_Phelps
./streamer d200mm_pure_H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_pure_H2O_Triniti.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_pure_H2O_1.0e6_Triniti