#!/bin/bash

cd ../

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.8.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Naidis_ioniza_fac0.8
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.75.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Naidis_ioniza_fac0.75