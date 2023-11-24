#!/bin/bash

cd ../

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_1.4e6_0.03H2O_Aints
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_1.4e6_0.10H2O_Aints



./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.5.txt -output%name=/export/scratch2/baohong/humid_air/chemistry_sensitivity/malagon_1.4e6_0.03H2O_Aints_ioniza_fac0.5
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.6.txt -output%name=/export/scratch2/baohong/humid_air/chemistry_sensitivity/malagon_1.4e6_0.03H2O_Aints_ioniza_fac0.6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.7.txt -output%name=/export/scratch2/baohong/humid_air/chemistry_sensitivity/malagon_1.4e6_0.03H2O_Aints_ioniza_fac0.7
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.75.txt -output%name=/export/scratch2/baohong/humid_air/chemistry_sensitivity/malagon_1.4e6_0.03H2O_Aints_ioniza_fac0.75
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.8.txt -output%name=/export/scratch2/baohong/humid_air/chemistry_sensitivity/malagon_1.4e6_0.03H2O_Aints_ioniza_fac0.8
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.9.txt -output%name=/export/scratch2/baohong/humid_air/chemistry_sensitivity/malagon_1.4e6_0.03H2O_Aints_ioniza_fac0.9
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/chemistry_sensitivity/malagon_1.4e6_0.03H2O_Aints_ioniza_fac1.0
