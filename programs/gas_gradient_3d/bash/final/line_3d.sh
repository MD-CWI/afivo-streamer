#!/bin/bash

cd ../

# from low-density to high-density region

# vary density ratio
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.0 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.0_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.2 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.2_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.4 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.4_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.6 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.6_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.8_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=2.0 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr2.0_sw0.01


# vary line angle
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.2 -shock_width=0.01 -line_coeff="1.5 1 0 -4" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line1.5+x-4z_1.6e6_dr1.2_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.2 -shock_width=0.01 -line_coeff="0.5 1 0 -2" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5+x-2z_1.6e6_dr1.2_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.2 -shock_width=0.01 -line_coeff="0 1 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/linex-z_1.6e6_dr1.2_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.2 -shock_width=0.01 -line_coeff="-0.5 2 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line-0.5+2x-z_1.6e6_dr1.2_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.2 -shock_width=0.01 -line_coeff="-1.5 4 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line-1.5+4x-z_1.6e6_dr1.2_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.2 -shock_width=0.01 -line_coeff="-0.5 1 0 0" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line-0.5+x_1.6e6_dr1.2_sw0.01


# vary shock width
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.015 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.8_sw0.015

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.02 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.8_sw0.02

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.03 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.8_sw0.03

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.06 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.8_sw0.06

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.1 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr1.8_sw0.1


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


# density changes gradually
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=2 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr2_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=3 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr3_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=4 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr4_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=8 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr8_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=12 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr12_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=20 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr20_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=25 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr25_sw0.35




# from high-density to low-density region

# vary density ratio
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=0.9 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr0.9_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr0.8_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr0.7_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=0.6 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5-z_1.6e6_dr0.6_sw0.01


# vary line angle
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="1.5 1 0 -4" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line1.5+x-4z_1.6e6_dr0.7_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="0.5 1 0 -2" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.5+x-2z_1.6e6_dr0.7_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="0 1 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/linex-z_1.6e6_dr0.7_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="-0.5 2 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line-0.5+2x-z_1.6e6_dr0.7_sw0.01

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="-1.5 4 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line-1.5+4x-z_1.6e6_dr0.7_sw0.01


