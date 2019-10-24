# mkdir -p ./output/final/standard0.26
# ./streamer positive_2d_left.cfg -output%name="output/final/standard0.26/standard0.26"

# mkdir -p ./output/final/ions/ions0.1
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ion=0.1 -output%name="output/final/ions/ions0.1/0.1sion"
# 
# mkdir -p ./output/final/ions/ions0.3
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ion=0.3 -output%name="output/final/ions/ions0.3/0.3sion"

# mkdir -p ./output/final/ions/ions0.5
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ion=0.5 -output%name="output/final/ions/ions0.5/0.5sion"
# 
# mkdir -p ./output/final/ions/ions0.0
# ./streamer positive_2d_left.cfg -dielectric%gamma_se_ion=0.0 -output%name="output/final/ions/ions0.0/0.26nose"

mkdir -p ./output/final/ions/1mmions0.5
./streamer positive_2d_left.cfg -dielectric%gamma_se_ion=0.5 -output%name="output/final/ions/1mmions0.5/1mm0.5sion" -seed_rel_r0="0.275E+00 1.0E+00"  -seed_rel_r1="0.275E+00 0.95E+00" -end_time=20e-9
