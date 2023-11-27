#!/bin/bash

cd ../

./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Kawaguchi.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.10H2O_Bourdon
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.10H2O_Aints
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.10H2O_Naidis  -photoi_helmh%coeffs="6.97075E+08  8.58360E+04  3.82738E+06" -photoi_helmh%lambdas="3.30404E+04  2.83311E+03  5.18507E+03"