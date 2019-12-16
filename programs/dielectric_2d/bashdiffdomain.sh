#in the middle
mkdir -p ./output/final/middle/10/left
./streamer positive_2d_left.cfg -output%name="output/final/middle/10/left/10left" -seed_rel_r0="0.5 1.0" -seed_rel_r1="0.5 0.95" -head%rmax=1 -head%rmin=0 -end_time=20e-9 -domain_len="10e-3 20e-3"  -field_amplitude=-3E+06

mkdir -p ./output/final/middle/10/gas
./streamer positive_2d_left.cfg -output%name="output/final/middle/10/gas/10gas" -seed_rel_r0="0.5 1.0" -seed_rel_r1="0.5 0.95" -head%rmax=1 -head%rmin=0 -dielectric_type=gas -end_time=20e-9 -domain_len="10e-3 20e-3" -field_amplitude=-3E+06

mkdir -p ./output/final/middle/20/left
./streamer positive_2d_left.cfg -output%name="output/final/middle/20/left/20left" -seed_rel_r0="0.5 1.0" -seed_rel_r1="0.5 0.95" -head%rmax=1 -head%rmin=0 -end_time=20e-9 -domain_len="20e-3 20e-3" -field_amplitude=-3E+06

mkdir -p ./output/final/middle/20/gas
./streamer positive_2d_left.cfg -output%name="output/final/middle/20/gas/20gas" -seed_rel_r0="0.5 1.0" -seed_rel_r1="0.5 0.95" -head%rmax=1 -head%rmin=0 -dielectric_type=gas -end_time=20e-9 -domain_len="20e-3 20e-3" -field_amplitude=-3E+06
 

