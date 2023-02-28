#!/bin/bash

cd ../

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.4 -shock_width=0.01 -line_coeff="-0.5 2 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line-0.5+2x-z_1.6e6_dr1.4_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.4 -shock_width=0.01 -line_coeff="-1.5 4 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line-1.5+4x-z_1.6e6_dr1.4_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.4 -shock_width=0.01 -line_coeff="-0.5 1 0 0" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line-0.5+x_1.6e6_dr1.4_sw0.01

