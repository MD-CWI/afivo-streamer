#!/bin/bash

cd ../

# effect of H2O cross sections

./streamer d200mm_0.00H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.00H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.00H2O_Bourdon

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Bourdon
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Itikawa_v3.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Bourdon_Itikawa_v3
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Morgan.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Bourdon_Morgan
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Bourdon_Phelps
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Triniti.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Bourdon_Triniti

./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Kawaguchi.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.10H2O_Bourdon
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Itikawa_v3.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.10H2O_Bourdon_Itikawa_v3
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Morgan.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.10H2O_Bourdon_Morgan
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.10H2O_Bourdon_Phelps
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Triniti.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.10H2O_Bourdon_Triniti



# effect of photoionization model

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Naidis
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi.txt -photoi%quenching_pressure=1.44392e-2 -photoi_helmh%coeffs="8.75226E+08  4.36991E+06  9.62902E+04" -photoi_helmh%lambdas="3.25733E+04  3.59107E+03  1.21351E+03" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Aints

./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.10H2O_Naidis
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.10H2O_Kawaguchi.txt -photoi%quenching_pressure=5.79662e-3 -photoi_helmh%coeffs="1.04484E+09  5.17925E+06  1.13966E+05" -photoi_helmh%lambdas="3.33108E+04  4.23899E+03  1.85941E+03" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.10H2O_Aints



# effect of chemistry

./streamer d200mm_0.00H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/baohong_0.00H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/baohong_1.0e6_0.00H2O_Bourdon

./streamer d200mm_0.00H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/komuro_0.00H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/komuro_1.0e6_0.00H2O_Bourdon
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/komuro_1.0e6_0.03H2O_Naidis
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/komuro_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/komuro_1.0e6_0.10H2O_Naidis
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/komuro_0.03H2O_Kawaguchi.txt -photoi%quenching_pressure=1.44392e-2 -photoi_helmh%coeffs="8.75226E+08  4.36991E+06  9.62902E+04" -photoi_helmh%lambdas="3.25733E+04  3.59107E+03  1.21351E+03" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/komuro_1.0e6_0.03H2O_Aints
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/komuro_0.10H2O_Kawaguchi.txt -photoi%quenching_pressure=5.79662e-3 -photoi_helmh%coeffs="1.04484E+09  5.17925E+06  1.13966E+05" -photoi_helmh%lambdas="3.33108E+04  4.23899E+03  1.85941E+03" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/komuro_1.0e6_0.10H2O_Aints

./streamer d200mm_0.00H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.00H2O_Phelps.txt -photoi_helmh%author="Bourdon-3" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_1.0e6_0.00H2O_Bourdon
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_1.0e6_0.03H2O_Naidis
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_Kawaguchi.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_1.0e6_0.10H2O_Naidis
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.03H2O_Kawaguchi.txt -photoi%quenching_pressure=1.44392e-2 -photoi_helmh%coeffs="8.75226E+08  4.36991E+06  9.62902E+04" -photoi_helmh%lambdas="3.25733E+04  3.59107E+03  1.21351E+03" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_1.0e6_0.03H2O_Aints
./streamer d200mm_0.10H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/starikovskiy_0.10H2O_Kawaguchi.txt -photoi%quenching_pressure=5.79662e-3 -photoi_helmh%coeffs="1.04484E+09  5.17925E+06  1.13966E+05" -photoi_helmh%lambdas="3.33108E+04  4.23899E+03  1.85941E+03" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/starikovskiy_1.0e6_0.10H2O_Aints

 

# ionization sensitivity

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.9.txt -photoi%quenching_pressure=1.44392e-2 -photoi_helmh%coeffs="8.75226E+08  4.36991E+06  9.62902E+04" -photoi_helmh%lambdas="3.25733E+04  3.59107E+03  1.21351E+03" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Aints_ioniza_fac0.9
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.8.txt -photoi%quenching_pressure=1.44392e-2 -photoi_helmh%coeffs="8.75226E+08  4.36991E+06  9.62902E+04" -photoi_helmh%lambdas="3.25733E+04  3.59107E+03  1.21351E+03" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Aints_ioniza_fac0.8
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.75.txt -photoi%quenching_pressure=1.44392e-2 -photoi_helmh%coeffs="8.75226E+08  4.36991E+06  9.62902E+04" -photoi_helmh%lambdas="3.25733E+04  3.59107E+03  1.21351E+03" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Aints_ioniza_fac0.75
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.7.txt -photoi%quenching_pressure=1.44392e-2 -photoi_helmh%coeffs="8.75226E+08  4.36991E+06  9.62902E+04" -photoi_helmh%lambdas="3.25733E+04  3.59107E+03  1.21351E+03" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Aints_ioniza_fac0.7
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.6.txt -photoi%quenching_pressure=1.44392e-2 -photoi_helmh%coeffs="8.75226E+08  4.36991E+06  9.62902E+04" -photoi_helmh%lambdas="3.25733E+04  3.59107E+03  1.21351E+03" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Aints_ioniza_fac0.6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.5.txt -photoi%quenching_pressure=1.44392e-2 -photoi_helmh%coeffs="8.75226E+08  4.36991E+06  9.62902E+04" -photoi_helmh%lambdas="3.25733E+04  3.59107E+03  1.21351E+03" -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Aints_ioniza_fac0.5
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.4.txt -photoi%quenching_pressure=1.44392e-2 -photoi_helmh%coeffs="8.75226E+08  4.36991E+06  9.62902E+04" -photoi_helmh%lambdas="3.25733E+04  3.59107E+03  1.21351E+03" -end_time=5e-7 -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Aints_ioniza_fac0.4
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.25.txt -photoi%quenching_pressure=1.44392e-2 -photoi_helmh%coeffs="8.75226E+08  4.36991E+06  9.62902E+04" -photoi_helmh%lambdas="3.25733E+04  3.59107E+03  1.21351E+03" -end_time=1e-6 -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Aints_ioniza_fac0.25

./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.8.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Naidis_ioniza_fac0.8
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.75.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Naidis_ioniza_fac0.75
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.6.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Naidis_ioniza_fac0.6
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.5.txt -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Naidis_ioniza_fac0.5
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.4.txt -end_time=5e-7 -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Naidis_ioniza_fac0.4
./streamer d200mm_0.03H2O.cfg -field_given_by="field -1.0e6" -input_data%file=../../transport_data/humid_air/malagon_0.03H2O_Kawaguchi_ioniza_fac0.25.txt -end_time=1e-6 -output%name=/export/scratch2/baohong/humid_air/streamer_comparison/malagon_1.0e6_0.03H2O_Naidis_ioniza_fac0.25


