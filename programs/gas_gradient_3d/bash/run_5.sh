#!/bin/bash

cd ../

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.02 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.8_sw0.02

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.04 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.8_sw0.04

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.06 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.8_sw0.06