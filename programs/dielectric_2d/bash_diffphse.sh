mkdir -p ./output/phse/nophse/
./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.0 -dielectric%gamma_se_ph_lowenergy=0.0 -output%name="output/phse/nophse/nophse"

mkdir -p ./output/phse/high0.5_low0.5/
./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.5 -dielectric%gamma_se_ph_lowenergy=0.5 -output%name="output/phse/high0.5_low0.5/0.5away_high0.5_low0.5" 

mkdir -p ./output/phse/high0.5_low0/
./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0.5 -dielectric%gamma_se_ph_lowenergy=0 -output%name="output/phse/high0.5_low0/0.5away_high0.5_low0"

mkdir -p ./output/phse/high0_low0.5/
./streamer positive_2d_left.cfg -dielectric%gamma_se_ph_highenergy=0 -dielectric%gamma_se_ph_lowenergy=0.5 -output%name="output/phse/high0_low0.5/0.5away_high0_low0.5"
