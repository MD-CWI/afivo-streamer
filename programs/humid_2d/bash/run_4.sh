#!/bin/bash

cd ../

./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_1.0e6_0.10H2O_Naidis

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.9.txt -photoi%quenching_pressure=1.44392e-2 -photoi_helmh%coeffs="8.75226E+08  4.36991E+06  9.62902E+04" -photoi_helmh%lambdas="3.25733E+04  3.59107E+03  1.21351E+03" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Aints_ioniza_fac0.9
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.8.txt -photoi%quenching_pressure=1.44392e-2 -photoi_helmh%coeffs="8.75226E+08  4.36991E+06  9.62902E+04" -photoi_helmh%lambdas="3.25733E+04  3.59107E+03  1.21351E+03" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Aints_ioniza_fac0.8