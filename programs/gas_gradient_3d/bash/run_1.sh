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



# from low-density to high-density region (outside hemisphere)

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=2.0 -shock_width=0.01 -sphere_center="0.5 0.5 1.0" -sphere_radius=0.5 -density_ratio_outside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.6e6_dr2.0_sw0.01_sc0.5-1.0_sr0.5

