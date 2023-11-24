#!/bin/bash

cd ../

./streamer d200mm_pure_H2O.cfg -input_data%file=../../transport_data/humid_air/malagon_pure_H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_0.9e6_pure_H2O_Kawaguchi
./streamer d200mm_pure_H2O.cfg -input_data%file=../../transport_data/humid_air/malagon_pure_H2O_Itikawa.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_0.9e6_pure_H2O_Itikawa
./streamer d200mm_pure_H2O.cfg -input_data%file=../../transport_data/humid_air/malagon_pure_H2O_Itikawa_v2.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_0.9e6_pure_H2O_Itikawa_v2
./streamer d200mm_pure_H2O.cfg -input_data%file=../../transport_data/humid_air/malagon_pure_H2O_Morgan.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_0.9e6_pure_H2O_Morgan
./streamer d200mm_pure_H2O.cfg -input_data%file=../../transport_data/humid_air/malagon_pure_H2O_Phelps.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_0.9e6_pure_H2O_Phelps
./streamer d200mm_pure_H2O.cfg -input_data%file=../../transport_data/humid_air/malagon_pure_H2O_Triniti.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_0.9e6_pure_H2O_Triniti