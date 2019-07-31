#1mmaway head%rmax=0.28
mkdir -p ./output/final/1away/high0_low0/
./streamer positive_2d_left.cfg -output%name="output/final/1away/high0_low0/1away_high0_low0" -seed_rel_r0="0.275 1.0" -seed_rel_r1="0.275 0.95" -head%rmax=0.285 -head%npoints=1500

mkdir -p ./output/final/1away/high0.1_low0.01/
./streamer positive_2d_left.cfg -output%name="output/final/1away/high0.1_low0.01/1away_high0.1_low0.01" -seed_rel_r0="0.275 1.0" -seed_rel_r1="0.275 0.95" -head%rmax=0.285 -head%npoints=1500 -dielectric%gamma_se_ph_highenergy=0.1 -dielectric%gamma_se_ph_lowenergy=0.01

mkdir -p ./output/final/1away/high0.5_low0.01/
./streamer positive_2d_left.cfg -output%name="output/final/1away/high0.5_low0.01/1away_high0.5_low0.01" -seed_rel_r0="0.275 1.0" -seed_rel_r1="0.275 0.95" -head%rmax=0.285 -head%npoints=1500 -dielectric%gamma_se_ph_highenergy=0.5 -dielectric%gamma_se_ph_lowenergy=0.01

mkdir -p ./output/final/1away/high0.5_low0.5/
./streamer positive_2d_left.cfg -output%name="output/final/1away/high0.5_low0.5/1away_high0.5_low0.5" -seed_rel_r0="0.275 1.0" -seed_rel_r1="0.275 0.95" -head%rmax=0.285 -head%npoints=1500 -dielectric%gamma_se_ph_highenergy=0.5 -dielectric%gamma_se_ph_lowenergy=0.5

mkdir -p ./output/final/1away/high0.1_low0/
./streamer positive_2d_left.cfg -output%name="output/final/1away/high0.1_low0/1away_high0.1_low0" -seed_rel_r0="0.275 1.0" -seed_rel_r1="0.275 0.95" -head%rmax=0.285 -head%npoints=1500 -dielectric%gamma_se_ph_highenergy=0.1 -dielectric%gamma_se_ph_lowenergy=0

mkdir -p ./output/final/1away/high0.5_low0/
./streamer positive_2d_left.cfg -output%name="output/final/1away/high0.5_low0/1away_high0.5_low0" -seed_rel_r0="0.275 1.0" -seed_rel_r1="0.275 0.95" -head%rmax=0.285 -head%npoints=1500 -dielectric%gamma_se_ph_highenergy=0.5 -dielectric%gamma_se_ph_lowenergy=0





