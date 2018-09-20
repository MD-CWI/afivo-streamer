# Test effect of parameters for adx = 1.0

# Increase number of multigrid cycles
# ../streamer_2d cyl_pos.cfg -refine_adx=1.0 -multigrid_num_vcycles=4 -simulation_name+=_adx1.0_4vcycle
../streamer_2d cyl_pos.cfg photoi.cfg -refine_adx=1.0 -multigrid_num_vcycles=4 -simulation_name+=_adx1.0_4vcycle

# Keep refinement
# ../streamer_2d cyl_pos.cfg -refine_adx=1.0 -derefine_dx=1e-7 -simulation_name+=_adx1.0_keep_ref
../streamer_2d cyl_pos.cfg photoi.cfg -refine_adx=1.0 -derefine_dx=1e-7 -simulation_name+=_adx1.0_keep_ref

# Use half the time stap
# ../streamer_2d cyl_pos.cfg -refine_adx=1.0 -dt_safety_factor=0.45 -simulation_name+=_adx1.0_halfdt
../streamer_2d cyl_pos.cfg photoi.cfg -refine_adx=1.0 -dt_safety_factor=0.45 -simulation_name+=_adx1.0_halfdt

# Different prolongation method
# ../streamer_2d cyl_pos.cfg -refine_adx=1.0 -prolong_density=linear -simulation_name+=_adx1.0_prolong_linear
../streamer_2d cyl_pos.cfg photoi.cfg -refine_adx=1.0 -prolong_density=linear -simulation_name+=_adx1.0_prolong_linear

# For photoionization: different number of steps per update
../streamer_2d cyl_pos.cfg photoi.cfg -refine_adx=1.0 -photoi%per_steps=1 -simulation_name+=_adx1.0_1step
../streamer_2d cyl_pos.cfg photoi.cfg -refine_adx=1.0 -photoi%per_steps=4 -simulation_name+=_adx1.0_4step
../streamer_2d cyl_pos.cfg photoi.cfg -refine_adx=1.0 -photoi%per_steps=8 -simulation_name+=_adx1.0_8step

# Change refinement criterion (alpha * dx < value)
# ../streamer_2d cyl_pos.cfg -refine_adx=2.0 -simulation_name+=_adx2.0
../streamer_2d cyl_pos.cfg photoi.cfg -refine_adx=2.0 -simulation_name+=_adx2.0
# ../streamer_2d cyl_pos.cfg -refine_adx=1.0 -simulation_name+=_adx1.0
../streamer_2d cyl_pos.cfg photoi.cfg -refine_adx=1.0 -simulation_name+=_adx1.0
# ../streamer_2d cyl_pos.cfg -refine_adx=0.5 -simulation_name+=_adx0.5
../streamer_2d cyl_pos.cfg photoi.cfg -refine_adx=0.5 -simulation_name+=_adx0.5
# ../streamer_2d cyl_pos.cfg -refine_adx=0.25 -simulation_name+=_adx0.25
../streamer_2d cyl_pos.cfg photoi.cfg -refine_adx=0.25 -simulation_name+=_adx0.25

# Uniform refinement
# ../streamer_2d cyl_pos.cfg -refine_adx=1.0e10 -refine_regions_dr=6.4e-6 -refine_regions_rmin="0.0 7.0e-3" -refine_regions_rmax="0.5e-3 11.e-3" -refine_regions_tstop=1.0 -simulation_name+=_unif_6.0
../streamer_2d cyl_pos.cfg photoi.cfg -refine_adx=1.0e10 -refine_regions_dr=6.4e-6 -refine_regions_rmin="0.0 7.0e-3" -refine_regions_rmax="0.5e-3 11.e-3" -refine_regions_tstop=1.0 -simulation_name+=_unif_6.0

# ../streamer_2d cyl_pos.cfg -refine_adx=1.0e10 -refine_regions_dr=3.2e-6 -refine_regions_rmin="0.0 7.0e-3" -refine_regions_rmax="0.5e-3 11.e-3" -refine_regions_tstop=1.0 -simulation_name+=_unif_3.0
../streamer_2d cyl_pos.cfg photoi.cfg -refine_adx=1.0e10 -refine_regions_dr=3.2e-6 -refine_regions_rmin="0.0 7.0e-3" -refine_regions_rmax="0.5e-3 11.e-3" -refine_regions_tstop=1.0 -simulation_name+=_unif_3.0

# ../streamer_2d cyl_pos.cfg -refine_adx=1.0e10 -refine_regions_dr=1.6e-6 -refine_regions_rmin="0.0 7.0e-3" -refine_regions_rmax="0.5e-3 11.e-3" -refine_regions_tstop=1.0 -simulation_name+=_unif_1.5
../streamer_2d cyl_pos.cfg photoi.cfg -refine_adx=1.0e10 -refine_regions_dr=1.6e-6 -refine_regions_rmin="0.0 7.0e-3" -refine_regions_rmax="0.5e-3 11.e-3" -refine_regions_tstop=1.0 -simulation_name+=_unif_1.5
