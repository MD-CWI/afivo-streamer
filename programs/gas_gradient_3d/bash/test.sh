#!/bin/bash

cd ../

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -dt_min=1e-15 -refine_min_dx=0 -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/test_1.6e6_dr1.8_dx0

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -dt_min=1e-15 -refine_min_dx=0 -refine_electrode_dx=1e-5 -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/test_1.6e6_dr1.8_refine_elec

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -dt_min=1e-15 -refine_min_dx=0 -derefine_dx=1e-5 -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/test_1.6e6_dr1.8_derefine

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -dt_min=1e-15 -refine_min_dx=0 -derefine_dx=1e-5 -refine_electrode_dx=1e-5 -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/test_1.6e6_dr1.8_spacing


./streamer d10mm_3d_new.cfg -field_given_by="field -1.4e6" -density_ratio=1.8 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -dt_min=1e-15 -refine_min_dx=0 -derefine_dx=1e-5 -refine_electrode_dx=1e-5 -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/test_1.4e6_dr1.8

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -dt_min=1e-15 -background_density=5e11 -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/test_1.6e6_dr1.8_bd5e11

./streamer d10mm_3d_seed.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -dt_min=1e-15 -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/test_1.6e6_dr1.8_seed

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -dt_min=1e-15 -rng_seed="1000 100000 1000 100000" -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/test_1.6e6_dr1.8_rng

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -dt_min=1e-15 -photoi%method="helmholtz" -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/test_1.6e6_dr1.8_helmholtz

./streamer d10mm_3d_new.cfg -field_given_by="field -1.6e6" -density_ratio=1.8 -shock_width=0.01 -line_coeff="0.5 0 0 -1" -dt_min=1e-15 -fixes%source_factor="flux" -output%name=/export/scratch2/baohong/gas_gradients/3d_sphere/test_1.6e6_dr1.8_flux


