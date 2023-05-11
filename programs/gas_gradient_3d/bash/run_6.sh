#!/bin/bash

cd ../

# from low-density to high-density region (outside hemisphere)

# vary density ratio
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.2 -shock_width=0.01 -sphere_center="0.5 0.5 1.0" -sphere_radius=0.5 -density_ratio_outside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.6e6_dr1.2_sw0.01_sc0.5-1.0_sr0.5

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.3 -shock_width=0.01 -sphere_center="0.5 0.5 1.0" -sphere_radius=0.5 -density_ratio_outside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.6e6_dr1.3_sw0.01_sc0.5-1.0_sr0.5

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.4 -shock_width=0.01 -sphere_center="0.5 0.5 1.0" -sphere_radius=0.5 -density_ratio_outside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.6e6_dr1.4_sw0.01_sc0.5-1.0_sr0.5

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.6 -shock_width=0.01 -sphere_center="0.5 0.5 1.0" -sphere_radius=0.5 -density_ratio_outside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.6e6_dr1.6_sw0.01_sc0.5-1.0_sr0.5

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.8 -shock_width=0.01 -sphere_center="0.5 0.5 1.0" -sphere_radius=0.5 -density_ratio_outside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.6e6_dr1.8_sw0.01_sc0.5-1.0_sr0.5

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=2.0 -shock_width=0.01 -sphere_center="0.5 0.5 1.0" -sphere_radius=0.5 -density_ratio_outside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.6e6_dr2.0_sw0.01_sc0.5-1.0_sr0.5

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=2.2 -shock_width=0.01 -sphere_center="0.5 0.5 1.0" -sphere_radius=0.5 -density_ratio_outside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.6e6_dr2.2_sw0.01_sc0.5-1.0_sr0.5


