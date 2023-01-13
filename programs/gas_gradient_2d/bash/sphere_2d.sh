#!/bin/bash

cd ../


# low density outside sphere

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.5 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_2.0e6_dr0.8_sw0.01_sc0.5-0.5_sr0.1_F

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.55 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_2.0e6_dr0.8_sw0.01_sc0.55-0.5_sr0.1_F

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.6 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_2.0e6_dr0.8_sw0.01_sc0.6-0.5_sr0.1_F

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_2.0e6_dr0.8_sw0.01_sc0.5-0.5_sr0.2_F

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.6 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_2.0e6_dr0.8_sw0.01_sc0.6-0.5_sr0.2_F

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.7 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_2.0e6_dr0.8_sw0.01_sc0.7-0.5_sr0.2_F



./streamer d20mm_2d.cfg -field_given_by="field -1.8e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.5 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.8e6_dr0.7_sw0.01_sc0.5-0.5_sr0.1_F

./streamer d20mm_2d.cfg -field_given_by="field -1.8e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.55 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.8e6_dr0.7_sw0.01_sc0.55-0.5_sr0.1_F

./streamer d20mm_2d.cfg -field_given_by="field -1.8e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.6 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.8e6_dr0.7_sw0.01_sc0.6-0.5_sr0.1_F

./streamer d20mm_2d.cfg -field_given_by="field -1.8e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.8e6_dr0.7_sw0.01_sc0.5-0.5_sr0.2_F

./streamer d20mm_2d.cfg -field_given_by="field -1.8e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.6 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.8e6_dr0.7_sw0.01_sc0.6-0.5_sr0.2_F

./streamer d20mm_2d.cfg -field_given_by="field -1.8e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.7 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.8e6_dr0.7_sw0.01_sc0.7-0.5_sr0.2_F



./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.5 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.01_sc0.5-0.5_sr0.1_F

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.55 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.01_sc0.55-0.5_sr0.1_F

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.6 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.01_sc0.6-0.5_sr0.1_F

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.01_sc0.5-0.5_sr0.2_F

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.6 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.01_sc0.6-0.5_sr0.2_F

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.7 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.01_sc0.7-0.5_sr0.2_F



# increase shock width

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.02 -sphere_center="0.5 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.02_sc0.5-0.5_sr0.1_F

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.04 -sphere_center="0.5 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.04_sc0.5-0.5_sr0.1_F

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.045 -sphere_center="0.5 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.045_sc0.5-0.5_sr0.1_F

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.05 -sphere_center="0.5 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.05_sc0.5-0.5_sr0.1_F

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.1 -sphere_center="0.5 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.1_sc0.5-0.5_sr0.1_F


./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.015 -sphere_center="0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.015_sc0.5-0.5_sr0.2_F

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.02 -sphere_center="0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.02_sc0.5-0.5_sr0.2_F

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.05 -sphere_center="0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.05_sc0.5-0.5_sr0.2_F

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.1 -sphere_center="0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=F -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.1_sc0.5-0.5_sr0.2_F





# low density inside sphere

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.5 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_2.0e6_dr0.8_sw0.01_sc0.5-0.5_sr0.1_T

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.55 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_2.0e6_dr0.8_sw0.01_sc0.55-0.5_sr0.1_T

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.6 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_2.0e6_dr0.8_sw0.01_sc0.6-0.5_sr0.1_T

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_2.0e6_dr0.8_sw0.01_sc0.5-0.5_sr0.2_T

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.6 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_2.0e6_dr0.8_sw0.01_sc0.6-0.5_sr0.2_T

./streamer d20mm_2d.cfg -field_given_by="field -2.0e6" -gradient_type="sphere" -density_ratio=0.8 -shock_width=0.01 -sphere_center="0.7 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_2.0e6_dr0.8_sw0.01_sc0.7-0.5_sr0.2_T



./streamer d20mm_2d.cfg -field_given_by="field -1.9e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.5 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.9e6_dr0.7_sw0.01_sc0.5-0.5_sr0.1_T

./streamer d20mm_2d.cfg -field_given_by="field -1.9e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.55 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.9e6_dr0.7_sw0.01_sc0.55-0.5_sr0.1_T

./streamer d20mm_2d.cfg -field_given_by="field -1.9e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.6 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.9e6_dr0.7_sw0.01_sc0.6-0.5_sr0.1_T

./streamer d20mm_2d.cfg -field_given_by="field -1.9e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.9e6_dr0.7_sw0.01_sc0.5-0.5_sr0.2_T

./streamer d20mm_2d.cfg -field_given_by="field -1.9e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.6 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.9e6_dr0.7_sw0.01_sc0.6-0.5_sr0.2_T

./streamer d20mm_2d.cfg -field_given_by="field -1.9e6" -gradient_type="sphere" -density_ratio=0.7 -shock_width=0.01 -sphere_center="0.7 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.9e6_dr0.7_sw0.01_sc0.7-0.5_sr0.2_T



./streamer d20mm_2d.cfg -field_given_by="field -1.9e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.5 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.9e6_dr0.6_sw0.01_sc0.5-0.5_sr0.1_T

./streamer d20mm_2d.cfg -field_given_by="field -1.9e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.55 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.9e6_dr0.6_sw0.01_sc0.55-0.5_sr0.1_T

./streamer d20mm_2d.cfg -field_given_by="field -1.9e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.6 0.5" -sphere_radius=0.1 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.9e6_dr0.6_sw0.01_sc0.6-0.5_sr0.1_T

./streamer d20mm_2d.cfg -field_given_by="field -1.9e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.5 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.9e6_dr0.6_sw0.01_sc0.5-0.5_sr0.2_T

./streamer d20mm_2d.cfg -field_given_by="field -1.9e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.6 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.9e6_dr0.6_sw0.01_sc0.6-0.5_sr0.2_T

./streamer d20mm_2d.cfg -field_given_by="field -1.9e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.7 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.9e6_dr0.6_sw0.01_sc0.7-0.5_sr0.2_T


# change shock width

./streamer d20mm_2d.cfg -field_given_by="field -1.9e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.02 -sphere_center="0.6 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.9e6_dr0.6_sw0.02_sc0.6-0.5_sr0.2_T

./streamer d20mm_2d.cfg -field_given_by="field -1.9e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.04 -sphere_center="0.6 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.9e6_dr0.6_sw0.04_sc0.6-0.5_sr0.2_T


./streamer d20mm_2d.cfg -field_given_by="field -1.9e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.02 -sphere_center="0.7 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.9e6_dr0.6_sw0.02_sc0.7-0.5_sr0.2_T

./streamer d20mm_2d.cfg -field_given_by="field -1.9e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.04 -sphere_center="0.7 0.5" -sphere_radius=0.2 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.9e6_dr0.6_sw0.04_sc0.7-0.5_sr0.2_T


# half sphere

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.5 0.8" -sphere_radius=0.2 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.01_sc0.5-0.8_sr0.2_T

./streamer d20mm_2d.cfg -field_given_by="field -1.5e6" -gradient_type="sphere" -density_ratio=0.6 -shock_width=0.01 -sphere_center="0.5 1.0" -sphere_radius=0.5 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.5e6_dr0.6_sw0.01_sc0.5-1.0_sr0.5_T


./streamer d20mm_2d.cfg -field_given_by="field -1.2e6" -gradient_type="sphere" -density_ratio=0.5 -shock_width=0.01 -sphere_center="0.5 0.8" -sphere_radius=0.2 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.2e6_dr0.5_sw0.01_sc0.5-0.8_sr0.2_T

./streamer d20mm_2d.cfg -field_given_by="field -1.2e6" -gradient_type="sphere" -density_ratio=0.5 -shock_width=0.01 -sphere_center="0.5 1.0" -sphere_radius=0.5 -density_ratio_inside_sphere=T -output%name=/export/scratch2/baohong/gas_gradients/2d_sphere/pos_EL4mm_1.2e6_dr0.5_sw0.01_sc0.5-1.0_sr0.5_T


