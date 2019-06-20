mkdir -p ./output/final/loc0.25/high0_low0/
./streamer positive_2d_left.cfg -output%name="output/final/loc0.25/high0_low0/0.25high0_low0"

mkdir -p ./output/fianl/loc0.25/high0.1_low0.01/
./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.1 -dielectric%gamma_se_ph_lowenergy=0.01 -output%name="output/final/loc0.25/high0.1_low0.01/0.25high0.1_low0.01"

mkdir -p ./output/final/loc0.25/high0.5_low0.01/
./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.5 -dielectric%gamma_se_ph_lowenergy=0.01 -output%name="output/final/loc0.25/high0.5_low0.01/0.25high0.5_low0.01"

mkdir -p ./output/final/loc0.25/high0.5_low0.5/
./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.5 -dielectric%gamma_se_ph_lowenergy=0.5 -output%name="output/final/loc0.25/high0.5_low0.5/0.25high0.5_low0.5"


#0.5mm away

mkdir -p ./output/final/0.5away/high0_low0/
./streamer positive_2d_left.cfg -output%name="output/final/0.5away/high0_low0/0.5away_high0_low0" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"

mkdir -p ./output/fianl/0.5away/high0.1_low0.01/
./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.1 -dielectric%gamma_se_ph_lowenergy=0.01 -output%name="output/final/0.5away/high0.1_low0.01/0.5away_high0.1_low0.01" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"

mkdir -p ./output/final/0.5away/high0.5_low0.01/
./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.5 -dielectric%gamma_se_ph_lowenergy=0.01 -output%name="output/final/0.5away/high0.5_low0.01/0.5away_high0.5_low0.01" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"

mkdir -p ./output/final/0.5away/high0.5_low0.5/
./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.5 -dielectric%gamma_se_ph_lowenergy=0.5 -output%name="output/final/0.5away/high0.5_low0.5/0.5away_high0.5_low0.5" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"

# #5mmaway
# mkdir -p ./output/final/5away/high0_low0/
# ./streamer positive_2d_left.cfg -output%name="output/final/5away/high0_low0/5away_high0_low0" -seed_rel_r0="0.375 1.0" -seed_rel_r1="0.375 0.95"
