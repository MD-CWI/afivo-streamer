# #0mmaway head%rmax=12mm
# mkdir -p ./output/negative/0mmaway
# ./streamer positive_2d_left.cfg -output%name="output/negative/0mmaway/n0mmaway" -seed_rel_r0="0.25 1.0" -seed_rel_r1="0.25 0.95" -head%rmax=0.3 -head%npoints=1000  -field_amplitude=2.5E+06
# 
# #1mmaway head%rmax=13mm
# mkdir -p ./output/negative/1mmaway
# ./streamer positive_2d_left.cfg -output%name="output/negative/1mmaway/n1mmaway" -seed_rel_r0="0.275 1.0" -seed_rel_r1="0.275 0.95" -head%rmax=0.325 -head%npoints=1000  -field_amplitude=2.5E+06
# 
# #2mmaway head%rmax=14mm
# mkdir -p ./output/negative/2mmaway
# ./streamer positive_2d_left.cfg -output%name="output/negative/2mmaway/n2mmaway" -seed_rel_r0="0.3 1.0" -seed_rel_r1="0.3 0.95" -head%rmax=0.35 -head%npoints=1000  -field_amplitude=2.5E+06

# #5mmaway head%rmax=18mm
# mkdir -p ./output/negative/5mmaway
# ./streamer positive_2d_left.cfg -output%name="output/negative/5mmaway/n5mmaway" -seed_rel_r0="0.375 1.0" -seed_rel_r1="0.375 0.95" -head%rmax=0.45 -head%npoints=1000  -field_amplitude=2.5E+06

mkdir -p ./output/negative/gas
./streamer positive_2d_left.cfg -output%name="output/negative/gas/ngas" -seed_rel_r0="0.5 1.0" -seed_rel_r1="0.5 0.95" -head%rmax=0.75 -head%rmin=0.25 -head%npoints=1000 -dielectric_type=gas -field_amplitude=2.5E+06

