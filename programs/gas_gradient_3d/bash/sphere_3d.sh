#!/bin/bash

cd ../

# from low-density to high-density region (inside sphere)

# vary density ratio
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.2 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr1.2_sw0.01_sc0.5-0.5_sr0.2

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.4 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr1.4_sw0.01_sc0.5-0.5_sr0.2

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.6 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr1.6_sw0.01_sc0.5-0.5_sr0.2

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.8 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr1.8_sw0.01_sc0.5-0.5_sr0.2


# vary sphere location
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.4 -shock_width=0.01 -sphere_center="0.6 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr1.4_sw0.01_sc0.6-0.5_sr0.2

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.4 -shock_width=0.01 -sphere_center="0.7 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr1.4_sw0.01_sc0.7-0.5_sr0.2




# from high-density to low-density region (inside sphere)

# vary density ratio
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=0.9 -shock_width=0.01 -sphere_center="0.7 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr0.9_sw0.01_sc0.7-0.5_sr0.2

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.7 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr0.8_sw0.01_sc0.7-0.5_sr0.2

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.7 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr0.7_sw0.01_sc0.7-0.5_sr0.2

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.7 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr0.6_sw0.01_sc0.7-0.5_sr0.2


# vary sphere location
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.6 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr0.7_sw0.01_sc0.6-0.5_sr0.2

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr0.7_sw0.01_sc0.5-0.5_sr0.2




# from low-density to high-density region (outside hemisphere)

# vary density ratio
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.2 -shock_width=0.01 -sphere_center="0.5 0.5 1.0" -sphere_radius=0.5 -density_ratio_outside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.6e6_dr1.2_sw0.01_sc0.5-1.0_sr0.5

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.4 -shock_width=0.01 -sphere_center="0.5 0.5 1.0" -sphere_radius=0.5 -density_ratio_outside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.6e6_dr1.4_sw0.01_sc0.5-1.0_sr0.5

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.6 -shock_width=0.01 -sphere_center="0.5 0.5 1.0" -sphere_radius=0.5 -density_ratio_outside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.6e6_dr1.6_sw0.01_sc0.5-1.0_sr0.5

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.8 -shock_width=0.01 -sphere_center="0.5 0.5 1.0" -sphere_radius=0.5 -density_ratio_outside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.6e6_dr1.8_sw0.01_sc0.5-1.0_sr0.5

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=2.0 -shock_width=0.01 -sphere_center="0.5 0.5 1.0" -sphere_radius=0.5 -density_ratio_outside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.6e6_dr2.0_sw0.01_sc0.5-1.0_sr0.5

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=2.2 -shock_width=0.01 -sphere_center="0.5 0.5 1.0" -sphere_radius=0.5 -density_ratio_outside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/out_1.6e6_dr2.2_sw0.01_sc0.5-1.0_sr0.5



