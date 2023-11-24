#!/bin/bash

cd ../

#./streamer d200mm_0.03H2O_photoi_custom.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_photoifit-0.03H2O
#./streamer d200mm_0.03H2O_photoi_custom.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_photoifit-0.00H2O -photoi_helmh%coeffs="6.10470E+06 8.07758E+08 1.51443E+05" -photoi_helmh%lambdas="4.03631E+03 3.57734E+04 1.06251E+03"

#./streamer d200mm_0.03H2O_photoi_custom.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.4e6_photoifit-0.03H2O
#./streamer d200mm_0.03H2O_photoi_custom.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.4e6_photoifit-0.00H2O -photoi_helmh%coeffs="6.10470E+06 8.07758E+08 1.51443E+05" -photoi_helmh%lambdas="4.03631E+03 3.57734E+04 1.06251E+03"

#./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_photoifit_0.10H2O
#./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_bourdon3_0.00H2O -photoi_helmh%author="Bourdon-3"

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -field_rod_r0="0 0.95725" -field_rod_radius=0.45e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_EL9mm_0.45mm

./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.8e6" -field_rod_r0="0 0.9525" -field_rod_radius=0.5e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_0.8e6_EL10mm_0.5mm
./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.85e6" -field_rod_r0="0 0.9525" -field_rod_radius=0.5e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_0.85e6_EL10mm_0.5mm