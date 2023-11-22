#!/bin/bash

cd ../

#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Triniti.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_Triniti

#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.0e6_ion -input_data%mobile_ions="N2_plus O2_plus H2O_plus H_plus OH_plus O_plus H2_plus O2_min O_min H_min OH_min O3_min H2O3_min H4O4_min H6O5_min N4_plus O4_plus H2O3_plus H3O_plus H5O2_plus H7O3_plus" -input_data%ion_mobilities="2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4 2.2e-4"


#./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Triniti.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_1.4e6_Triniti

./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.85e6" -field_rod_r0="0 0.95725" -field_rod_radius=0.45e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_0.85e6_EL9mm_0.45mm
./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.9e6" -field_rod_r0="0 0.95725" -field_rod_radius=0.45e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_0.9e6_EL9mm_0.45mm
./streamer d200mm_0.03H2O.cfg -field_given_by="field -0.95e6" -field_rod_r0="0 0.95725" -field_rod_radius=0.45e-3 -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_air_chemistry.txt -output%name=/export/scratch2/baohong/humid_air/comparison/d200_malagon_0.03H2O_0.95e6_EL9mm_0.45mm