#!/bin/bash

cd ../

./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.9e6" -input_data%file=../../transport_data/humid_air/komuro_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/komuro_0.9e6_0.10H2O_Aints
./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.9e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_0.9e6_0.10H2O_Aints

./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/komuro_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/komuro_1.4e6_0.10H2O_Aints