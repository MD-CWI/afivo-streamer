#!/bin/bash

cd ../

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.08 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.8_sw0.08

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.1 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.8_sw0.1

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=0.9 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr0.9_sw0.01

