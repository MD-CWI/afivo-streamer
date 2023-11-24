#!/bin/bash

# This is provided by Hema

# Coefficients in Bourdon-3
# coeffs  = [1117314.935, 28692377.5, 2748842283.0]
# lambdas = [4147.85, 10950.93, 66755.67]

# After multiplying them with pO2 and pO2^2 for pO2 = 0.2
# coeffs = [4.46925974e+04, 1.14769510e+06, 1.09953691e+08]
# lambdas = [829.57 ,  2190.186, 13351.134]

echo "Bourdon Coeffs: 4.46925974e+04, 1.14769510e+06, 1.09953691e+08"

echo "Bourdon Lambdas: 829.57 , 2190.186, 13351.134"

bourdon_fit="4.46925974e+04 829.57 1.14769510e+06 2190.186 1.09953691e+08 13351.134"
fit_range="2e-6 6e-3"
echo "-----------------Dry air, fit with Bourdon coeffs from literature"
./absorption_function.py -gases O2 -pressures 0.2 -n_modes 3 -fit_range $fit_range -fit_type log -show_curve $bourdon_fit -show_Zheleznyak

echo "-----------------Moist air, 1.5% water, fit with Naidis curve from literature"
./absorption_function.py -gases O2 H2O -pressures 0.197 0.015 -n_modes 3 -fit_range $fit_range -fit_type log -show_curve $bourdon_fit -show_Zheleznyak -show_Naidis_moist -show_Aints_moist

echo "-----------------Moist air, 3% water, fit with Naidis curve from literature"
./absorption_function.py -gases O2 H2O -pressures 0.194 0.03 -n_modes 3 -fit_range $fit_range -fit_type log -show_curve $bourdon_fit -show_Zheleznyak -show_Naidis_moist -show_Aints_moist

echo "-----------------Moist air, 10% water, fit with Naidis curve from literature"
./absorption_function.py -gases O2 H2O -pressures 0.18 0.1 -n_modes 3 -fit_range $fit_range -fit_type log -show_curve $bourdon_fit -show_Zheleznyak -show_Naidis_moist -show_Aints_moist



# This is provided by Baohong (with the new version of absoprtion_function.py script by Jannis)

# 3% H2O
python absorption_function.py -gases O2 H2O -pressures 0.194 0.03 -fit_range 2e-6 8e-3 -fit_type log -fit_what Aints -ptot_for_quenching 1 -show_Zheleznyak -show_Naidis_moist -show_Aints_moist -figure_name absorption_fuction_0.03H20_Aints.png
python absorption_function.py -gases O2 H2O -pressures 0.194 0.03 -fit_range 2e-6 8e-3 -fit_type log -fit_what Zheleznyak-H2O -ptot_for_quenching 1 -show_Zheleznyak -show_Naidis_moist -show_Aints_moist -figure_name absorption_fuction_0.03H20_Naidis.png

# 10% H2O
python absorption_function.py -gases O2 H2O -pressures 0.18 0.1 -fit_range 2e-6 8e-3 -fit_type log -fit_what Aints -ptot_for_quenching 1 -show_Zheleznyak -show_Naidis_moist -show_Aints_moist -figure_name absorption_fuction_0.10H20_Aints.png
python absorption_function.py -gases O2 H2O -pressures 0.18 0.1 -fit_range 2e-6 8e-3 -fit_type log -fit_what Zheleznyak-H2O -ptot_for_quenching 1 -show_Zheleznyak -show_Naidis_moist -show_Aints_moist -figure_name absorption_fuction_0.10H20_Naidis.png



