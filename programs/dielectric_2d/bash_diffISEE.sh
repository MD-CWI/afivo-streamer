mkdir -p ./output/ionISEE/ionISEE0.0
./streamer positive_2d_left.cfg -dielectric%gamma_se_ion=0.0 -output%name="output/ionISEE/ionISEE0.0/noISEE" -seed_rel_r0="0.26 1.0"  -seed_rel_r1="0.26 0.95"

mkdir -p ./output/ionISEE/ionISEE0.5
./streamer positive_2d_left.cfg -dielectric%gamma_se_ion=0.5 -output%name="output/ionISEE/ionISEE0.5/0.5ISEE" -seed_rel_r0="0.26 1.0"  -seed_rel_r1="0.26 0.95"
