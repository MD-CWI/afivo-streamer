# mkdir -p ./output/final/nphi/high0_low0/
# ./streamer positive_2d_left.cfg -output%name="output/final/nphi/high0_low0/0.25high0_low0" -end_time=20.0e-9 -photoi%enabled_ingas=f

# mkdir -p ./output/final/nphi/high0.1_low0.01/
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.1 -dielectric%gamma_se_ph_lowenergy=0.01 -output%name="output/final/nphi/high0.1_low0.01/0.25high0.1_low0.01" -end_time=20.0e-9
# 
# mkdir -p ./output/final/nphi/high0.5_low0.01/
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.5 -dielectric%gamma_se_ph_lowenergy=0.01 -output%name="output/final/nphi/high0.5_low0.01/0.25high0.5_low0.01" -end_time=20.0e-9

mkdir -p ./output/final/nphi/high0.5_low0.5/
./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.5 -dielectric%gamma_se_ph_lowenergy=0.5 -output%name="output/final/nphi/high0.5_low0.5/0.25high0.5_low0.5" -end_time=20.0e-9 -photoi%enabled_ingas=f
# 
# mkdir -p ./output/final/nphi/high0.5_low0/
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.5 -dielectric%gamma_se_ph_lowenergy=0 -output%name="output/final/nphi/high0.5_low0/0.25high0.5_low0" -end_time=20.0e-9 -photoi%enabled_ingas=f
# 
# mkdir -p ./output/final/nphi/high0_low0.5/
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0 -dielectric%gamma_se_ph_lowenergy=0.5 -output%name="output/final/nphi/high0_low0.5/0.25high0_low0.5" -end_time=20.0e-9 -photoi%enabled_ingas=f

mkdir -p ./output/final/nphi/high0.1_low0.1/
./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.1 -dielectric%gamma_se_ph_lowenergy=0.1 -output%name="output/final/nphi/high0.1_low0.1/0.25high0.1_low0.1" -end_time=20.0e-9 -photoi%enabled_ingas=f

mkdir -p ./output/final/nphi/high0.3_low0.3/
./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.3 -dielectric%gamma_se_ph_lowenergy=0.3 -output%name="output/final/nphi/high0.3_low0.3/0.25high0.3_low0.3" -end_time=20.0e-9 -photoi%enabled_ingas=f

# mkdir -p ./output/final/loc0.25/high0_low0.01/
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0 -dielectric%gamma_se_ph_lowenergy=0.01 -output%name="output/final/loc0.25/high0_low0.01/0.25high0_low0.01" -end_time=20.0e-9
# 
# mkdir -p ./output/final/loc0.25/high0.1_low0/
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.1 -dielectric%gamma_se_ph_lowenergy=0 -output%name="output/final/loc0.25/high0.1_low0/0.25high0.1_low0" -end_time=20.0e-9


# # 0.5mm away
# 
# mkdir -p ./output/final/0.5away/high0_low0/
# ./streamer positive_2d_left.cfg -output%name="output/final/0.5away/high0_low0/0.5away_high0_low0" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"
# 
# mkdir -p ./output/final/0.5away/high0.1_low0.01/
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.1 -dielectric%gamma_se_ph_lowenergy=0.01 -output%name="output/final/0.5away/high0.1_low0.01/0.5away_high0.1_low0.01" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"
# 
# mkdir -p ./output/final/0.5away/high0.5_low0.01/
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.5 -dielectric%gamma_se_ph_lowenergy=0.01 -output%name="output/final/0.5away/high0.5_low0.01/0.5away_high0.5_low0.01" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"
# 
# mkdir -p ./output/final/0.5away/high0.5_low0.5/
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.5 -dielectric%gamma_se_ph_lowenergy=0.5 -output%name="output/final/0.5away/high0.5_low0.5/0.5away_high0.5_low0.5" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"
# 
# mkdir -p ./output/final/0.5away/high0.5_low0/
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.5 -dielectric%gamma_se_ph_lowenergy=0 -output%name="output/final/0.5away/high0.5_low0/0.5away_high0.5_low0" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95" -end_time=15.0e-9
# 
# mkdir -p ./output/final/0.5away/high0_low0.5/
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0 -dielectric%gamma_se_ph_lowenergy=0.5 -output%name="output/final/0.5away/high0_low0.5/0.5away_high0_low0.5" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95" -end_time=15.0e-9
# 
# mkdir -p ./output/final/0.5away/high0_low0.01/
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0 -dielectric%gamma_se_ph_lowenergy=0.01 -output%name="output/final/0.5away/high0_low0.01/0.5away_high0_low0.01" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95" -end_time=15.0e-9
# 
# mkdir -p ./output/final/0.5away/high0.1_low0/
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.1 -dielectric%gamma_se_ph_lowenergy=0 -output%name="output/final/0.5away/high0.1_low0/0.5away_high0.1_low0" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95" -end_time=15.0e-9
# 
# # 1mmaway
# mkdir -p ./output/final/1away/high0_low0.5/
# ./streamer positive_2d_left.cfg -output%name="output/final/1away/high0_low0.5/1away_high0_low0.5" -seed_rel_r0="0.275 1.0" -seed_rel_r1="0.275 0.95" -head%rmax=0.285 -head%npoints=1500 -dielectric%gamma_se_ph_highenergy=0 -dielectric%gamma_se_ph_lowenergy=0.5
# 
# mkdir -p ./output/final/1away/high0_low0.01/
# ./streamer positive_2d_left.cfg -output%name="output/final/1away/high0_low0.01/1away_high0_low0.01" -seed_rel_r0="0.275 1.0" -seed_rel_r1="0.275 0.95" -head%rmax=0.285 -head%npoints=1500 -dielectric%gamma_se_ph_highenergy=0 -dielectric%gamma_se_ph_lowenergy=0.01


