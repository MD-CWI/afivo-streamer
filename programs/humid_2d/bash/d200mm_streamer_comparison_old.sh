#!/bin/bash

cd ../

# effect of cross sections

./streamer d200mm_0.00H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.00H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.00H2O_Bourdon

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Bourdon
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Itikawa.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Bourdon_Itikawa
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Itikawa_v2.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Bourdon_Itikawa_v2
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Morgan.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Bourdon_Morgan
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Bourdon_Phelps
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Triniti.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Bourdon_Triniti

./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Kawaguchi.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.10H2O_Bourdon
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Itikawa.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.10H2O_Bourdon_Itikawa
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Itikawa_v2.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.10H2O_Bourdon_Itikawa_v2
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Morgan.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.10H2O_Bourdon_Morgan
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.10H2O_Bourdon_Phelps
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Triniti.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.10H2O_Bourdon_Triniti



# effect of photoionization

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Aints
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Naidis  -photoi_helmh%coeffs="8.62978E+04  3.89447E+06  7.57232E+08" -photoi_helmh%lambdas="1.50562E+03  3.87500E+03  3.24957E+04"
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.10H2O_Aints
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.10H2O_Naidis  -photoi_helmh%coeffs="6.97075E+08  8.58360E+04  3.82738E+06" -photoi_helmh%lambdas="3.30404E+04  2.83311E+03  5.18507E+03"

./streamer d200mm_0.00H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.00H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.00H2O_Bourdon
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.03H2O_Bourdon
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.03H2O_Aints
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.03H2O_Naidis  -photoi_helmh%coeffs="8.62978E+04  3.89447E+06  7.57232E+08" -photoi_helmh%lambdas="1.50562E+03  3.87500E+03  3.24957E+04"
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Kawaguchi.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.10H2O_Bourdon
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.10H2O_Aints
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.10H2O_Naidis  -photoi_helmh%coeffs="6.97075E+08  8.58360E+04  3.82738E+06" -photoi_helmh%lambdas="3.30404E+04  2.83311E+03  5.18507E+03"



# effect of chemistry

./streamer d200mm_0.00H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/baohong_0.00H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/baohong_1.0e6_0.00H2O_Bourdon
./streamer d200mm_0.00H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/komuro_0.00H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/komuro_1.0e6_0.00H2O_Bourdon
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/komuro_1.0e6_0.03H2O_Aints
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/komuro_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/komuro_1.0e6_0.10H2O_Aints
./streamer d200mm_0.00H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.00H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_1.0e6_0.00H2O_Bourdon
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_1.0e6_0.03H2O_Aints
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_1.0e6_0.10H2O_Aints

./streamer d200mm_0.00H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/baohong_0.00H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/baohong_1.4e6_0.00H2O_Bourdon
./streamer d200mm_0.00H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/komuro_0.00H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/komuro_1.4e6_0.00H2O_Bourdon
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/komuro_1.4e6_0.03H2O_Aints
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/komuro_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/komuro_1.4e6_0.10H2O_Aints
./streamer d200mm_0.00H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.00H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_1.4e6_0.00H2O_Bourdon
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_1.4e6_0.03H2O_Aints
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_1.4e6_0.10H2O_Aints


 
# ionization sensitivity

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.5.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.03H2O_Aints_ioniza_fac0.5
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.6.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.03H2O_Aints_ioniza_fac0.6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.7.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.03H2O_Aints_ioniza_fac0.7
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.75.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.03H2O_Aints_ioniza_fac0.75
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.8.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.03H2O_Aints_ioniza_fac0.8
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.9.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.03H2O_Aints_ioniza_fac0.9
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.4e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.4e6_0.03H2O_Aints_ioniza_fac1.0





