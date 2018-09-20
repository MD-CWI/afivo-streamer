# Currently, only include these two basic tests. The rest should be run on HPC
# systems (or you have to be quite patient)

../streamer_3d 3d_pos.cfg -refine_adx=1.0 -simulation_name+=_adx1.0
../streamer_3d 3d_pos.cfg photoi.cfg -refine_adx=1.0 -simulation_name+=_adx1.0

# Test effect of parameters for adx = 1.0

# Increase number of multigrid cycles
# ../streamer_3d 3d_pos.cfg -refine_adx=1.0 -multigrid_num_vcycles=4 -simulation_name+=_adx1.0_4vcycle
# ../streamer_3d 3d_pos.cfg photoi.cfg -refine_adx=1.0 -multigrid_num_vcycles=4 -simulation_name+=_adx1.0_4vcycle

# Use half the time stap
# ../streamer_3d 3d_pos.cfg -refine_adx=1.0 -dt_safety_factor=0.45 -simulation_name+=_adx1.0_halfdt
# ../streamer_3d 3d_pos.cfg photoi.cfg -refine_adx=1.0 -dt_safety_factor=0.45 -simulation_name+=_adx1.0_halfdt

# Different prolongation method
# ../streamer_3d 3d_pos.cfg -refine_adx=1.0 -prolong_density=linear -simulation_name+=_adx1.0_prolong_linear
# ../streamer_3d 3d_pos.cfg photoi.cfg -refine_adx=1.0 -prolong_density=linear -simulation_name+=_adx1.0_prolong_linear

# For photoionization: different number of steps per update
# ../streamer_3d 3d_pos.cfg photoi.cfg -refine_adx=1.0 -photoi%per_steps=2 -simulation_name+=_adx1.0_2step
# ../streamer_3d 3d_pos.cfg photoi.cfg -refine_adx=1.0 -photoi%per_steps=4 -simulation_name+=_adx1.0_4step
# ../streamer_3d 3d_pos.cfg photoi.cfg -refine_adx=1.0 -photoi%per_steps=8 -simulation_name+=_adx1.0_8step

# Change refinement criterion (alpha * dx < value)
# ../streamer_3d 3d_pos.cfg -refine_adx=2.0 -simulation_name+=_adx2.0
# ../streamer_3d 3d_pos.cfg photoi.cfg -refine_adx=2.0 -simulation_name+=_adx2.0
# ../streamer_3d 3d_pos.cfg -refine_adx=1.0 -simulation_name+=_adx1.0
# ../streamer_3d 3d_pos.cfg photoi.cfg -refine_adx=1.0 -simulation_name+=_adx1.0
# ../streamer_3d 3d_pos.cfg -refine_adx=0.5 -simulation_name+=_adx0.5
# ../streamer_3d 3d_pos.cfg photoi.cfg -refine_adx=0.5 -simulation_name+=_adx0.5
