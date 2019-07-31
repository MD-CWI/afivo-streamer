#in the middle
# mkdir -p ./output/final/middle/left/
# ./streamer positive_2d_left.cfg -output%name="output/final/middle/left/ml_high0_low0" -seed_rel_r0="0.5 1.0" -seed_rel_r1="0.5 0.95" -head%rmax=0.75 -head%rmin=0.25 -end_time=50e-9

# mkdir -p ./output/final/middle/two/
# ./streamer positive_2d_left.cfg -output%name="output/final/middle/two/mt_high0_low0" -seed_rel_r0="0.5 1.0" -seed_rel_r1="0.5 0.95" -head%rmax=0.75 -head%rmin=0.25 -dielectric_type=left_right -end_time=33e-9
# 
# mkdir -p ./output/final/middle/gas/
# ./streamer positive_2d_left.cfg -output%name="output/final/middle/gas/mg_high0_low0" -seed_rel_r0="0.5 1.0" -seed_rel_r1="0.5 0.95" -head%rmax=0.75 -head%rmin=0.25 -dielectric_type=gas -end_time=33e-9
# 
# #5mmaway
# mkdir -p ./output/final/5away/high0_low0/
# ./streamer positive_2d_left.cfg -output%name="output/final/5away/high0_low0/5away_high0_low0" -seed_rel_r0="0.375 1.0" -seed_rel_r1="0.375 0.95" -head%rmax=0.75 -head%rmin=0.25 -end_time=33e-9

# #3mmaway
# mkdir -p ./output/final/3away/high0_low0/
# ./streamer positive_2d_left.cfg -output%name="output/final/3away/high0_low0/3away_high0_low0" -seed_rel_r0="0.325 1.0" -seed_rel_r1="0.325 0.95" -head%rmax=1 -head%rmin=0 -end_time=25e-9
# 
# #2mmaway
# mkdir -p ./output/final/2away/high0_low0/
# ./streamer positive_2d_left.cfg -output%name="output/final/2away/high0_low0/2away_high0_low0" -seed_rel_r0="0.3 1.0" -seed_rel_r1="0.3 0.95" -head%rmax=1 -head%rmin=0 -end_time=25e-9

# #3mmaway
# mkdir -p ./output/final/3away/gas/
# ./streamer positive_2d_left.cfg -output%name="output/final/3away/gas/3awaygas_high0_low0" -seed_rel_r0="0.325 1.0" -seed_rel_r1="0.325 0.95" -head%rmax=1 -head%rmin=0 -end_time=25e-9 -dielectric_type=gas

#2mmaway
mkdir -p ./output/final/2away/gas/
./streamer positive_2d_left.cfg -output%name="output/final/2away/gas/2awaygas_high0_low0" -seed_rel_r0="0.3 1.0" -seed_rel_r1="0.3 0.95" -head%rmax=1 -head%rmin=0 -end_time=25e-9 -dielectric_type=gas



