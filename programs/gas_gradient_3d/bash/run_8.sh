#!/bin/bash

cd ../

./streamer d10mm_3d.cfg -field_given_by="field -1.2e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.2e6_dr0.8_sw0.01_sc0.5-0.5_sr0.2

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.5e6_dr0.7_sw0.01_sc0.5-0.5_sr0.2

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.5e6_dr0.6_sw0.01_sc0.5-0.5_sr0.2

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.03 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.5e6_dr0.6_sw0.03_sc0.5-0.5_sr0.2

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.05 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.5e6_dr0.6_sw0.05_sc0.5-0.5_sr0.2

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.1 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.5e6_dr0.6_sw0.1_sc0.5-0.5_sr0.2

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.2 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.5e6_dr0.6_sw0.2_sc0.5-0.5_sr0.2

./streamer d10mm_3d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.5e6_dr0.8_sw0.01_sc0.5-0.5_sr0.2



