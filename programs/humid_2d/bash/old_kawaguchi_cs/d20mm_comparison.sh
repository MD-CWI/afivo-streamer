#!/bin/bash

cd ../

./streamer d20mm_comparison.cfg -field_given_by="field -1.2e6" -input_data%file=../../transport_data/humid_air/baohong_dry_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/dry_baohong_1.2e6
./streamer d20mm_comparison.cfg -field_given_by="field -1.6e6" -input_data%file=../../transport_data/humid_air/baohong_dry_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/dry_baohong_1.6e6
./streamer d20mm_comparison.cfg -field_given_by="field -2.0e6" -input_data%file=../../transport_data/humid_air/baohong_dry_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/dry_baohong_2.0e6

./streamer d20mm_comparison.cfg -field_given_by="field -1.2e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/humid0.03_komuro_1.2e6
./streamer d20mm_comparison.cfg -field_given_by="field -1.6e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/humid0.03_komuro_1.6e6
./streamer d20mm_comparison.cfg -field_given_by="field -2.0e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/humid0.03_komuro_2.0e6

./streamer d20mm_comparison.cfg -field_given_by="field -1.2e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/humid0.03_malagon_1.2e6
./streamer d20mm_comparison.cfg -field_given_by="field -1.6e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/humid0.03_malagon_1.6e6
./streamer d20mm_comparison.cfg -field_given_by="field -2.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/humid0.03_malagon_2.0e6

./streamer d20mm_comparison.cfg -field_given_by="field -1.2e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/humid0.03_starikovskiy_1.2e6
./streamer d20mm_comparison.cfg -field_given_by="field -1.6e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/humid0.03_starikovskiy_1.6e6
./streamer d20mm_comparison.cfg -field_given_by="field -2.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/humid0.03_starikovskiy_2.0e6


# Old Morgan cross section data
./streamer d20mm_comparison.cfg -field_given_by="field -1.2e6" -input_data%file=../../transport_data/humid_air/old_phelps_morgan_cs/komuro_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/morgan_cs_humid0.03_komuro_1.2e6
./streamer d20mm_comparison.cfg -field_given_by="field -1.2e6" -input_data%file=../../transport_data/humid_air/old_phelps_morgan_cs/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/morgan_cs_humid0.03_malagon_1.2e6
./streamer d20mm_comparison.cfg -field_given_by="field -1.2e6" -input_data%file=../../transport_data/humid_air/old_phelps_morgan_cs/starikovskiy_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/morgan_cs_humid0.03_starikovskiy_1.2e6



