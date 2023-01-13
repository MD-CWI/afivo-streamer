#!/bin/bash

cd ../

# low density in left half

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="0.5 -1 0" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.8_line0.5-x



# low density in top/left top

./streamer d20mm_2d.cfg -field_given_by="field -1.2e6" -density_ratio=0.5 -line_coeff="-0.5 0 1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_1.2e6_dr0.5_line-0.5+y


./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -line_coeff="-0.5 0 1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_1.5e6_dr0.6_line-0.5+y

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -line_coeff="-0.5 0 1" -shock_width=0.015 -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_1.5e6_dr0.6_line-0.5+y_sw0.015

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -line_coeff="-0.5 0 1" -shock_width=0.02 -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_1.5e6_dr0.6_line-0.5+y_sw0.02

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -line_coeff="-0.5 0 1" -shock_width=0.05 -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_1.5e6_dr0.6_line-0.5+y_sw0.05

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -line_coeff="-0.5 0 1" -shock_width=0.1 -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_1.5e6_dr0.6_line-0.5+y_sw0.1

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -density_ratio=0.6 -line_coeff="-0.5 0 1" -shock_width=0.2 -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_1.5e6_dr0.6_line-0.5+y_sw0.2


./streamer d20mm_2d.cfg -field_given_by="field -1.8e6" -density_ratio=0.7 -line_coeff="-0.5 0 1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_1.8e6_dr0.7_line-0.5+y

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.7 -line_coeff="-0.5 0 1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.7_line-0.5+y

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.9 -line_coeff="-0.5 0 1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.9_line-0.5+y

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="-0.5 0 1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.8_line-0.5+y

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="-4.5 -1 10" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.8_line-4.5-x+10y

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="-2 -1 5" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.8_line-2-x+5y

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="-0.5 -1 2" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.8_line-0.5-x+2y

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="0 -1 1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.8_line-x+y

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="0.5 -2 1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.8_line0.5-2x+y



# low density in bottom/right bottom

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.7 -line_coeff="0.5 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.7_line0.5-y

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.9 -line_coeff="0.5 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.9_line0.5-y

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="0.5 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.8_line0.5-y

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="0.5 1 -2" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.8_line0.5+x-2y

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="0 1 -1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.8_linex-y

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="-0.5 2 -1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.8_line-0.5+2x-y

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="-2 5 -1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.8_line-2+5x-y

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -density_ratio=0.8 -line_coeff="-4.5 10 -1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/pos_EL4mm_2.0e6_dr0.8_line-4.5+10x-y



# negative

./streamer d20mm_2d.cfg -field_given_by="field 2.0e6" -density_ratio=0.8 -line_coeff="0 -1 1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/neg_EL4mm_2.0e6_dr0.8_line-x+y

./streamer d20mm_2d.cfg -field_given_by="field 2.0e6" -density_ratio=0.8 -line_coeff="0 1 -1" -output%name=/export/scratch2/baohong/gas_gradients/2d_line/neg_EL4mm_2.0e6_dr0.8_linex-y




