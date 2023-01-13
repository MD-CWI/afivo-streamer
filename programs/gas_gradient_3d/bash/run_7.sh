#!/bin/bash

cd ../

./streamer d10mm_3d.cfg -field_given_by="field 1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="-2 5 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/bottom_neg_1.5e6_dr0.8_sw0.01_line-2+5x-z

./streamer d10mm_3d.cfg -field_given_by="field 1.5e6" -density_ratio=0.8 -shock_width=0.01 -line_coeff="-0.5 -1 0 2" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/top_neg_1.5e6_dr0.8_sw0.01_line-0.5-x+2z


./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.5e6_dr0.8_sw0.01_sc0.5-0.5_sr0.2

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.6 0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.5e6_dr0.8_sw0.01_sc0.6-0.5_sr0.2

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.7 0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.5e6_dr0.8_sw0.01_sc0.7-0.5_sr0.2

./streamer d10mm_3d.cfg -field_given_by="field -1.8e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.8e6_dr0.8_sw0.01_sc0.5-0.5_sr0.2

./streamer d10mm_3d.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.6e6_dr0.8_sw0.01_sc0.5-0.5_sr0.2

./streamer d10mm_3d.cfg -field_given_by="field -1.4e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.4e6_dr0.8_sw0.01_sc0.5-0.5_sr0.2

