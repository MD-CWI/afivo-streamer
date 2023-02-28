#!/bin/bash

cd ../

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.4 -shock_width=0.01 -line_coeff="1.5 1 0 -4" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line1.5+x-4z_1.6e6_dr1.4_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.4 -shock_width=0.01 -line_coeff="0.5 1 0 -2" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5+x-2z_1.6e6_dr1.4_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.4 -shock_width=0.01 -line_coeff="0 1 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/linex-z_1.6e6_dr1.4_sw0.01

