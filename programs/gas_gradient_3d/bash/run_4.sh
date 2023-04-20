#!/bin/bash

cd ../

# from low-density to high-density region

# density changes gradually
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=4 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr4_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=8 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr8_sw0.35

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=12 -shock_width=0.35 -line_coeff="0.35 0 0 -1" -output%name=/export/scratch2/baohong/gas_gradients/3d_line/line0.35-z_1.6e6_dr12_sw0.35



# from low-density to high-density region (inside sphere)

# vary density ratio
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.2 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr1.2_sw0.01_sc0.5-0.5_sr0.2

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.4 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr1.4_sw0.01_sc0.5-0.5_sr0.2

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.6 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr1.6_sw0.01_sc0.5-0.5_sr0.2

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.8 -shock_width=0.01 -sphere_center="0.5 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr1.8_sw0.01_sc0.5-0.5_sr0.2


# vary sphere location
./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.4 -shock_width=0.01 -sphere_center="0.6 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr1.4_sw0.01_sc0.6-0.5_sr0.2

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -gradient_type="sphere" -density_ratio=1.4 -shock_width=0.01 -sphere_center="0.7 0.5 0.5" -sphere_radius=0.2 -density_ratio_outside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/in_1.6e6_dr1.4_sw0.01_sc0.7-0.5_sr0.2

