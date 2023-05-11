#!/bin/bash

cd ../

# from low-density to high-density region

# density changes gradually
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=2 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr2_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=3 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr3_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=4 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr4_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=6 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr6_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=8 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr8_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=10 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr10_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=20 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr20_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=30 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr30_sw0.35


