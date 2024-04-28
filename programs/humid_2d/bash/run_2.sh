#!/bin/bash

cd ../

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_Kawaguchi.txt -photoi%quenching_pressure=1.44392e-2 -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_1.0e6_0.03H2O_Aints

./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_Kawaguchi.txt -photoi%quenching_pressure=5.79662e-3 -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_1.0e6_0.10H2O_Aints

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.5.txt -photoi%quenching_pressure=1.44392e-2 -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.03H2O_Aints_ioniza_fac0.5

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.6.txt -photoi%quenching_pressure=1.44392e-2 -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.03H2O_Aints_ioniza_fac0.6