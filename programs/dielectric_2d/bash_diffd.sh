#0.5mmaway head%rmax=12mm
mkdir -p ./output/diffd/0.5mmaway
./streamer positive_2d_left.cfg -output%name="output/diffd/0.5mmaway/0.5mmaway" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95" -head%rmax=0.3

#1mmaway head%rmax=13mm
mkdir -p ./output/diffd/1mmaway
./streamer positive_2d_left.cfg -output%name="output/diffd/1mmaway/1mmaway" -seed_rel_r0="0.275 1.0" -seed_rel_r1="0.275 0.95" -head%rmax=0.325 -end_time=20e-9

#2mmaway head%rmax=14mm
mkdir -p ./output/diffd/2mmaway
./streamer positive_2d_left.cfg -output%name="output/diffd/2mmaway/2mmaway" -seed_rel_r0="0.3 1.0" -seed_rel_r1="0.3 0.95" -head%rmax=0.35 -end_time=20e-9

#5mmaway head%rmax=18mm
mkdir -p ./output/diffd/5mmaway
./streamer positive_2d_left.cfg -output%name="output/diffd/5mmaway/5mmaway" -seed_rel_r0="0.375 1.0" -seed_rel_r1="0.375 0.95" -head%rmax=0.45 -end_time=20e-9
