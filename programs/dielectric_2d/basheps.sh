mkdir -p ./output/final/eps/2.0
./streamer positive_2d_left.cfg -output%name="output/final/eps/2.0/2.0eps" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"  -dielectric_eps=2.0

mkdir -p ./output/final/eps/3.0
./streamer positive_2d_left.cfg -output%name="output/final/eps/3.0/3.0eps" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"  -dielectric_eps=3.0  -end_time=15.0e-9

mkdir -p ./output/final/eps/5.0
./streamer positive_2d_left.cfg -output%name="output/final/eps/5.0/5.0eps" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"  -dielectric_eps=5.0 -end_time=12.0e-9
