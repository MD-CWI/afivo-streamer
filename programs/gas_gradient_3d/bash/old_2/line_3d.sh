#!/bin/bash

cd ../


# low density in left half

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="0.5 -1 0 0" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/half_pos_1.5e6_dr0.8_sw0.01_line0.5-x




# low density in  bottom or bottom/right

# vary density ratio
./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.9 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/bottom_pos_1.5e6_dr0.9_sw0.01_line0.5-z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/bottom_pos_1.5e6_dr0.8_sw0.01_line0.5-z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/bottom_pos_1.5e6_dr0.7_sw0.01_line0.5-z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/bottom_pos_1.5e6_dr0.6_sw0.01_line0.5-z


# vary line angle
./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="0.5 1 0 -2" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/bottom_pos_1.5e6_dr0.8_sw0.01_line0.5+x-2z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="0 1 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/bottom_pos_1.5e6_dr0.8_sw0.01_linex-z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="-0.5 2 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/bottom_pos_1.5e6_dr0.8_sw0.01_line-0.5+2x-z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="-2 5 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/bottom_pos_1.5e6_dr0.8_sw0.01_line-2+5x-z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="-4.5 10 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/bottom_pos_1.5e6_dr0.8_sw0.01_line-4.5+10x-z




# low density in top or top/left

# vary density ratio
./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.9 -shock_width=0.01 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.5e6_dr0.9_sw0.01_line-0.5+z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.5e6_dr0.8_sw0.01_line-0.5+z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.5e6_dr0.7_sw0.01_line-0.5+z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -shock_width=0.01 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.5e6_dr0.6_sw0.01_line-0.5+z


# vary line angle
./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="-4.5 -1 0 10" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.5e6_dr0.8_sw0.01_line-4.5-x+10z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="-2 -1 0 5" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.5e6_dr0.8_sw0.01_line-2-x+5z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="-0.5 -1 0 2" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.5e6_dr0.8_sw0.01_line-0.5-x+2z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="0 -1 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.5e6_dr0.8_sw0.01_line-x+z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="0.5 -2 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.5e6_dr0.8_sw0.01_line0.5-2x+z


# vary background field
./streamer d10mm_3d.cfg -field_given_by="field -1.8e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.8e6_dr0.7_sw0.01_line-0.5+z

./streamer d10mm_3d.cfg -field_given_by="field -1.6e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.6e6_dr0.7_sw0.01_line-0.5+z

./streamer d10mm_3d.cfg -field_given_by="field -1.4e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.4e6_dr0.7_sw0.01_line-0.5+z

./streamer d10mm_3d.cfg -field_given_by="field -1.2e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.2e6_dr0.7_sw0.01_line-0.5+z

./streamer d10mm_3d.cfg -field_given_by="field -1.0e6" -density_ratio=0.7 -shock_width=0.01 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.0e6_dr0.7_sw0.01_line-0.5+z


# vary shock width
./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -shock_width=0.02 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.5e6_dr0.6_sw0.02_line-0.5+z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -shock_width=0.04 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.5e6_dr0.6_sw0.04_line-0.5+z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -shock_width=0.06 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.5e6_dr0.6_sw0.06_line-0.5+z

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -shock_width=0.1 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.5e6_dr0.6_sw0.1_line-0.5+z




# negative streamers

./streamer d10mm_3d.cfg -field_given_by="field 1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="-2 5 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/bottom_neg_1.5e6_dr0.8_sw0.01_line-2+5x-z

./streamer d10mm_3d.cfg -field_given_by="field 1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="-0.5 -1 0 2" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_neg_1.5e6_dr0.8_sw0.01_line-0.5-x+2z




# without source factor

./streamer d10mm_3d_nsf.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -shock_width=0.01 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_pos_1.5e6_dr0.6_sw0.01_line-0.5+z_nsf




