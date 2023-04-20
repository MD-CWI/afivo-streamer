#!/bin/bash

cd ../

# from low-density to high-density region

# vary gradient slope
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.3 -shock_width=0.005 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.3_sw0.005

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.35 -shock_width=0.005 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.35_sw0.005

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.65 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.65_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.7 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.7_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=2.2 -shock_width=0.02 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr2.2_sw0.02

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=2.3 -shock_width=0.02 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr2.3_sw0.02

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=2.4 -shock_width=0.02 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr2.4_sw0.02


# vary background field
./streamer d10mm_3d_new.cfg -field_given_by="field -1.2e6" -density_ratio=1.6 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.2e6_dr1.6_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.4e6" -density_ratio=1.6 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.4e6_dr1.6_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.5e6" -density_ratio=1.6 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.5e6_dr1.6_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.7e6" -density_ratio=1.6 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.7e6_dr1.6_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.8e6" -density_ratio=1.6 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.8e6_dr1.6_sw0.01

