#!/bin/bash

cd ../

# different transport data

./streamer d200mm_dry_air.cfg -field_given_by="field -0.8e6" -input_data%file=../../transport_data/humid_air/baohong_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_baohong_dry_0.8e6
./streamer d200mm_dry_air.cfg -field_given_by="field -1.2e6" -input_data%file=../../transport_data/humid_air/baohong_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_baohong_dry_1.2e6
./streamer d200mm_dry_air.cfg -field_given_by="field -1.6e6" -input_data%file=../../transport_data/humid_air/baohong_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_baohong_dry_1.6e6
./streamer d200mm_dry_air.cfg -field_given_by="field -2.0e6" -input_data%file=../../transport_data/humid_air/baohong_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_baohong_dry_2.0e6


./streamer d200mm_dry_air.cfg -field_given_by="field -0.8e6" -input_data%file=../../transport_data/humid_air/komuro_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_dry_0.8e6
./streamer d200mm_dry_air.cfg -field_given_by="field -1.2e6" -input_data%file=../../transport_data/humid_air/komuro_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_dry_1.2e6
./streamer d200mm_dry_air.cfg -field_given_by="field -1.6e6" -input_data%file=../../transport_data/humid_air/komuro_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_dry_1.6e6
./streamer d200mm_dry_air.cfg -field_given_by="field -2.0e6" -input_data%file=../../transport_data/humid_air/komuro_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_dry_2.0e6

./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.8e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_0.03H2O_0.8e6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.2e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_0.03H2O_1.2e6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.6e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_0.03H2O_1.6e6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -2.0e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_0.03H2O_2.0e6


./streamer d200mm_dry_air.cfg -field_given_by="field -0.8e6" -input_data%file=../../transport_data/humid_air/malagon_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_dry_0.8e6
./streamer d200mm_dry_air.cfg -field_given_by="field -1.2e6" -input_data%file=../../transport_data/humid_air/malagon_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_dry_1.2e6
./streamer d200mm_dry_air.cfg -field_given_by="field -1.6e6" -input_data%file=../../transport_data/humid_air/malagon_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_dry_1.6e6
./streamer d200mm_dry_air.cfg -field_given_by="field -2.0e6" -input_data%file=../../transport_data/humid_air/malagon_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_dry_2.0e6

./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.8e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_0.8e6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.2e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.2e6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.6e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.6e6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -2.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_2.0e6


./streamer d200mm_dry_air.cfg -field_given_by="field -0.8e6" -input_data%file=../../transport_data/humid_air/starikovskiy_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_dry_0.8e6
./streamer d200mm_dry_air.cfg -field_given_by="field -1.2e6" -input_data%file=../../transport_data/humid_air/starikovskiy_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_dry_1.2e6
./streamer d200mm_dry_air.cfg -field_given_by="field -1.6e6" -input_data%file=../../transport_data/humid_air/starikovskiy_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_dry_1.6e6
./streamer d200mm_dry_air.cfg -field_given_by="field -2.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_dry_2.0e6

./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.8e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.03H2O_0.8e6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.2e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.03H2O_1.2e6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.6e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.03H2O_1.6e6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -2.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.03H2O_2.0e6



# different cross sections

./streamer d200mm_dry_air.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_dry_1.0e6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Itikawa.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_Itikawa
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Itikawa_v2.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_Itikawa-v2
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Morgan.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_Morgan
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Phelps.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_Phelps
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Triniti.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_Triniti

./streamer d200mm_dry_air.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_dry_1.4e6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.4e6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Itikawa.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.4e6_Itikawa
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Itikawa_v2.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.4e6_Itikawa-v2
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Morgan.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.4e6_Morgan
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Phelps.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.4e6_Phelps
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Triniti.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.4e6_Triniti

./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_Kawaguchi
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Itikawa.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_Itikawa
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Itikawa_v2.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_Itikawa-v2
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Morgan.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_Morgan
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Phelps.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_Phelps
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Triniti.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_Triniti

./streamer d200mm_0.20H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.20H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.20H2O_1.0e6_Kawaguchi
./streamer d200mm_0.20H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.20H2O_Itikawa.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.20H2O_1.0e6_Itikawa
./streamer d200mm_0.20H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.20H2O_Itikawa_v2.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.20H2O_1.0e6_Itikawa-v2
./streamer d200mm_0.20H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.20H2O_Morgan.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.20H2O_1.0e6_Morgan
./streamer d200mm_0.20H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.20H2O_Phelps.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.20H2O_1.0e6_Phelps
./streamer d200mm_0.20H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.20H2O_Triniti.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.20H2O_1.0e6_Triniti



# different photoionization coefficients

./streamer d200mm_0.03H2O_photoi_custom.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_photoifit-0.03H2O
./streamer d200mm_0.03H2O_photoi_custom.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_photoifit-0.00H2O -photoi_helmh%coeffs="6.10470E+06 8.07758E+08 1.51443E+05" -photoi_helmh%lambdas="4.03631E+03 3.57734E+04 1.06251E+03"

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_0.03H2O_1.0e6
./streamer d200mm_0.03H2O_photoi_custom.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_0.03H2O_1.0e6_photoifit-0.03H2O

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.03H2O_1.0e6
./streamer d200mm_0.03H2O_photoi_custom.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.03H2O_1.0e6_photoifit-0.03H2O


./streamer d200mm_0.03H2O_photoi_custom.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.4e6_photoifit-0.03H2O
./streamer d200mm_0.03H2O_photoi_custom.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.4e6_photoifit-0.00H2O -photoi_helmh%coeffs="6.10470E+06 8.07758E+08 1.51443E+05" -photoi_helmh%lambdas="4.03631E+03 3.57734E+04 1.06251E+03"

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_0.03H2O_1.4e6
./streamer d200mm_0.03H2O_photoi_custom.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_komuro_0.03H2O_1.4e6_photoifit-0.03H2O

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.03H2O_1.4e6
./streamer d200mm_0.03H2O_photoi_custom.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.03H2O_1.4e6_photoifit-0.03H2O


./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_photoifit_0.10H2O
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_bourdon3_0.00H2O -photoi_helmh%author="Bourdon-3"

./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.4e6_photoifit_0.10H2O
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.4e6_bourdon3_0.00H2O -photoi_helmh%author="Bourdon-3"


./streamer d200mm_0.20H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.20H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.20H2O_1.0e6_photoifit_0.20H2O
./streamer d200mm_0.20H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.20H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.20H2O_1.0e6_bourdon3_0.00H2O -photoi_helmh%author="Bourdon-3"

./streamer d200mm_0.20H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.20H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.20H2O_1.4e6_photoifit_0.20H2O
./streamer d200mm_0.20H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.20H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.20H2O_1.4e6_bourdon3_0.00H2O -photoi_helmh%author="Bourdon-3"


./streamer d200mm_0.20H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.20H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.20H2O_1.0e6_photoifit_0.20H2O
./streamer d200mm_0.20H2O.cfg -field_given_by="field -1.0e6" -field_rod_r0="0 0.9525" -field_rod_radius=0.5e-3 -input_data%file=../../transport_data/humid_air/malagon_0.20H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.20H2O_1.0e6_photoifit_0.20H2O_EL10mm_0.5mm
./streamer d200mm_0.20H2O.cfg -field_given_by="field -1.0e6" -field_rod_r0="0 0.9525" -field_rod_radius=0.5e-3 -input_data%file=../../transport_data/humid_air/starikovskiy_0.20H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.20H2O_1.0e6_photoifit_0.20H2O_EL10mm_0.5mm


./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.10H2O_1.0e6_photoifit_0.10H2O
./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.85e6" -field_rod_r0="0 0.962" -field_rod_radius=0.4e-3 -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_0.85e6_photoifit_0.10H2O_EL8mm_0.4mm
./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.85e6" -field_rod_r0="0 0.962" -field_rod_radius=0.4e-3 -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.10H2O_0.85e6_photoifit_0.10H2O_EL8mm_0.4mm

./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.9e6" -field_rod_r0="0 0.962" -field_rod_radius=0.4e-3 -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_0.9e6_photoifit_0.10H2O_EL8mm_0.4mm
./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.9e6" -field_rod_r0="0 0.962" -field_rod_radius=0.4e-3 -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.10H2O_0.9e6_photoifit_0.10H2O_EL8mm_0.4mm

./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -field_rod_r0="0 0.962" -field_rod_radius=0.4e-3 -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_1.0e6_photoifit_0.10H2O_EL8mm_0.4mm
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -field_rod_r0="0 0.962" -field_rod_radius=0.4e-3 -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.10H2O_1.0e6_photoifit_0.10H2O_EL8mm_0.4mm

./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.95e6" -field_rod_r0="0 0.962" -field_rod_radius=0.4e-3 -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_0.95e6_photoifit_0.10H2O_EL8mm_0.4mm
./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.95e6" -field_rod_r0="0 0.962" -field_rod_radius=0.4e-3 -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.10H2O_0.95e6_photoifit_0.10H2O_EL8mm_0.4mm


./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.95e6" -field_rod_r0="0 0.96675" -field_rod_radius=0.35e-3 -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_0.95e6_photoifit_0.10H2O_EL7mm_0.35mm
./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.95e6" -field_rod_r0="0 0.96675" -field_rod_radius=0.35e-3 -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.10H2O_0.95e6_photoifit_0.10H2O_EL7mm_0.35mm

./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.95e6" -field_rod_r0="0 0.9715" -field_rod_radius=0.3e-3 -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_0.95e6_photoifit_0.10H2O_EL6mm_0.3mm
./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.95e6" -field_rod_r0="0 0.9715" -field_rod_radius=0.3e-3 -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.10H2O_0.95e6_photoifit_0.10H2O_EL6mm_0.3mm

./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.95e6" -field_rod_r0="0 0.97625" -field_rod_radius=0.25e-3 -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_0.95e6_photoifit_0.10H2O_EL5mm_0.25mm
./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.95e6" -field_rod_r0="0 0.97625" -field_rod_radius=0.25e-3 -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.10H2O_0.95e6_photoifit_0.10H2O_EL5mm_0.25mm

./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.95e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_0.95e6_photoifit_0.10H2O_EL4mm_0.2mm
./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.95e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_starikovskiy_0.10H2O_0.95e6_photoifit_0.10H2O_EL4mm_0.2mm
./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.9e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_0.9e6_photoifit_0.10H2O
./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.85e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_0.85e6_photoifit_0.10H2O
./streamer d200mm_0.10H2O.cfg -field_given_by="field -0.8e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.10H2O_0.8e6_photoifit_0.10H2O



# include ion motion

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_ion -input_data%mobile_ions="N2_plus O2_plus H2O_plus H_plus OH_plus O_plus H2_plus O_plus_plus O2_min O_min H_min OH_min O3_min H2O3_min H4O4_min H6O5_min N4_plus O4_plus H2O3_plus H3O_plus H5O2_plus H7O3_plus H9O4_plus" -input_data%ion_mobilities="2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4"



./streamer d200mm_dry_air.cfg -field_given_by="field -1.0e6" -gas%temperature=300 -input_data%file=../../transport_data/humid_air/malagon_dry_300K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_dry_1.0e6_300K
./streamer d200mm_dry_air.cfg -field_given_by="field -1.0e6" -gas%temperature=319.3 -input_data%file=../../transport_data/humid_air/malagon_dry_319.3K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_dry_1.0e6_319.3K
./streamer d200mm_dry_air.cfg -field_given_by="field -1.0e6" -gas%temperature=333.6 -input_data%file=../../transport_data/humid_air/malagon_dry_333.6K_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_dry_1.0e6_333.6K





