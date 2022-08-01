#!/bin/bash

# This script can be used to convert old chemical reactions to the new format
# A backup of the file with a .bak extension will be created

if (($# != 1)); then
   echo "Usage: $0 input_file.txt"
fi

sed -i'.bak' -e  '
s@constant@c1@g;
s@k2_func@c1@g;
s@linear@c1*(Td-c2)@g;
s@exp_v1@c1*exp(-(c2/(c3+Td))**2)@g;
s@exp_v2@c1*exp(-(Td/c2)**2)@g;
s@k1_func@c1*(300/Te)**c2@g;
s@k3_func@(c1*(kB_eV*Te+c2)**2-c3)*c4@g;
s@k4_func@c1*(Tg/300)**c2*exp(-c3/Tg)@g;
s@k5_func@c1*exp(-c2/Tg)@g;
s@exp_v1@c1*exp(-(c2/(c3+Td))**2)@g;
s@exp_v2@c1*exp(-(Td/c2)**2)@g;
s@k1_func@c1*(300/Te)**c2@g;
s@k3_func@(c1*(kB_eV*Te+c2)**2-c3)*c4@g;
s@k4_func@c1*(Tg/300)**c2*exp(-c3/Tg)@g;
s@k5_func@c1*exp(-c2/Tg)@g;
s@k6_func@c1*Tg**c2@g;
s@k7_func@c1*(Tg/c2)**c3@g;
s@k8_func@c1*(300/Tg)**c2@g;
s@k9_func@c1*exp(-c2*Tg)@g;
s@k10_func@10**(c1+c2*(Tg-300))@g;
s@k11_func@c1*(300/Tg)**c2*exp(-c3/Tg)@g;
s@k12_func@c1*Tg**c2*exp(-c3/Tg)@g;
s@k13_func@c1*exp(-(c2/(c3+Td))**c4)@g;
s@k14_func@c1*exp(-(Td/c2)**c3)@g;
s@k15_func@c1*exp(-(c2/(kb*(Tg+Td/c3)))**c4)@g;
' "$1"
