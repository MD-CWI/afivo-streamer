#!/bin/bash

cd ../

./streamer d10mm_3d_new.cfg -field_given_by="field -1.0e6" -density_ratio=1.6 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.0e6_dr1.6_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.2e6" -density_ratio=1.6 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.2e6_dr1.6_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.4e6" -density_ratio=1.6 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.4e6_dr1.6_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.8e6" -density_ratio=1.6 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.8e6_dr1.6_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -2.0e6" -density_ratio=1.6 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_2.0e6_dr1.6_sw0.01
