#!/bin/bash

cd ../

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr0.8_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr0.7_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=0.6 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr0.6_sw0.01

