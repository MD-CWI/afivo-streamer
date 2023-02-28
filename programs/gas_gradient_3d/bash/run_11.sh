#!/bin/bash

cd ../

./streamer d10mm_3d.cfg -field_given_by="field -1.0e6" -density_ratio=0.5 -shock_width=0.35 -line_coeff="-0.35 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/linear_pos_1.0e6_dr0.5_sw0.35_line-0.35+z

./streamer d10mm_3d.cfg -field_given_by="field -0.4e6" -density_ratio=0.2 -shock_width=0.35 -line_coeff="-0.35 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/linear_pos_0.4e6_dr0.2_sw0.35_line-0.35+z

./streamer d10mm_3d.cfg -field_given_by="field -0.2e6" -density_ratio=0.1 -shock_width=0.35 -line_coeff="-0.35 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/linear_pos_0.2e6_dr0.1_sw0.35_line-0.35+z

./streamer d10mm_3d.cfg -field_given_by="field -0.1e6" -density_ratio=0.05 -shock_width=0.35 -line_coeff="-0.35 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/linear_pos_0.1e6_dr0.05_sw0.35_line-0.35+z

./streamer d10mm_3d.cfg -field_given_by="field -0.04e6" -density_ratio=0.02 -shock_width=0.35 -line_coeff="-0.35 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/linear_pos_0.04e6_dr0.02_sw0.35_line-0.35+z

./streamer d10mm_3d.cfg -field_given_by="field -0.02e6" -density_ratio=0.01 -shock_width=0.35 -line_coeff="-0.35 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/linear_pos_0.02e6_dr0.01_sw0.35_line-0.35+z



./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=2 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch3/baohong/gas_gradients/3d_line/linear_pos_1.5e6_dr2_sw0.35_line0.35-z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=3 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/linear_pos_1.5e6_dr3_sw0.35_line0.35-z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=4 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/linear_pos_1.5e6_dr4_sw0.35_line0.35-z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=5 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/linear_pos_1.5e6_dr5_sw0.35_line0.35-z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=6 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/linear_pos_1.5e6_dr6_sw0.35_line0.35-z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=8 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/linear_pos_1.5e6_dr8_sw0.35_line0.35-z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=10 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/linear_pos_1.5e6_dr10_sw0.35_line0.35-z
