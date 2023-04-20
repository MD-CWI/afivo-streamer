#!/bin/bash

cd ../

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=20 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch3/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr20_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=25 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch3/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr25_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=28 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch3/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr28_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=30 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch3/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr30_sw0.35