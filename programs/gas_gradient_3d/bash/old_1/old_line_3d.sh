#!/bin/bash

cd ../

# low density in left half

./streamer d5mm_3d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="0.5 -1 0 0" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_2.0e6_dr0.8_line0.5-x



# low density in top/left top

./streamer d5mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_1.5e6_dr0.6_line-0.5+z

./streamer d5mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -line_coeff="-0.5 0 0 1" -shock_width=0.02 -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_1.5e6_dr0.6_line-0.5+z_sw0.02

./streamer d5mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -line_coeff="-0.5 0 0 1" -shock_width=0.04 -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_1.5e6_dr0.6_line-0.5+z_sw0.04

./streamer d5mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -line_coeff="-0.5 0 0 1" -shock_width=0.06 -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_1.5e6_dr0.6_line-0.5+z_sw0.06

./streamer d5mm_3d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -line_coeff="-0.5 0 0 1" -shock_width=0.1 -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_1.5e6_dr0.6_line-0.5+z_sw0.1


./streamer d5mm_3d.cfg -field_given_by="field -2.0e6" -density_ratio=0.7 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_2.0e6_dr0.7_line-0.5+z

./streamer d5mm_3d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="-0.5 0 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_2.0e6_dr0.8_line-0.5+z

./streamer d5mm_3d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="-4.5 -1 0 10" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_2.0e6_dr0.8_line-4.5-x+10z

./streamer d5mm_3d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="-2 -1 0 5" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_2.0e6_dr0.8_line-2-x+5z

./streamer d5mm_3d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="-0.5 -1 0 2" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_2.0e6_dr0.8_line-0.5-x+2z

./streamer d5mm_3d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="0 -1 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_2.0e6_dr0.8_line-x+z

./streamer d5mm_3d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="0.5 -2 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_2.0e6_dr0.8_line0.5-2x+z



# low density in bottom/right bottom

./streamer d5mm_3d.cfg -field_given_by="field -2.0e6" -density_ratio=0.7 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_2.0e6_dr0.7_line0.5-z

./streamer d5mm_3d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="0.5 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_2.0e6_dr0.8_line0.5-z

./streamer d5mm_3d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="0.5 1 0 -2" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_2.0e6_dr0.8_line0.5+x-2z

./streamer d5mm_3d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="0 1 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_2.0e6_dr0.8_linex-z

./streamer d5mm_3d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="-0.5 2 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_2.0e6_dr0.8_line-0.5+2x-z

./streamer d5mm_3d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="-2 5 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_2.0e6_dr0.8_line-2+5x-z

./streamer d5mm_3d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="-4.5 10 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_pos_EL1mm_2.0e6_dr0.8_line-4.5+10x-z



# negative

./streamer d5mm_3d.cfg -field_given_by="field 2.0e6" -density_ratio=0.8 -line_coeff="0 -1 0 1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_neg_EL1mm_2.0e6_dr0.8_line-x+z

./streamer d5mm_3d.cfg -field_given_by="field 2.0e6" -density_ratio=0.8 -line_coeff="0 1 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/3d_neg_EL1mm_2.0e6_dr0.8_linex-z




