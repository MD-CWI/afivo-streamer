# mkdir -p ./output/final/field/2.1
# ./streamer positive_2d_left.cfg -output%name="output/final/field/2.1/2.1field" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"  -field_amplitude=-2.1E+06
# 
# mkdir -p ./output/final/field/2.3
# ./streamer positive_2d_left.cfg -output%name="output/final/field/2.3/2.3field" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"  -field_amplitude=-2.3E+06
# 
# mkdir -p ./output/final/field/2.8
# ./streamer positive_2d_left.cfg -output%name="output/final/field/2.8/2.8field" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"  -field_amplitude=-2.8E+06

mkdir -p ./output/negative/field/2.1
./streamer positive_2d_left.cfg -output%name="output/negative/field/2.1/n2.1field" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"  -field_amplitude=2.1E+06

mkdir -p ./output/negative/field/2.3
./streamer positive_2d_left.cfg -output%name="output/negative/field/2.3/n2.3field" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"  -field_amplitude=2.3E+06

mkdir -p ./output/negative/field/2.7
./streamer positive_2d_left.cfg -output%name="output/negative/field/2.7/n2.7field" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"  -field_amplitude=2.7E+06


mkdir -p ./output/negative/eps/2.0
./streamer positive_2d_left.cfg -output%name="output/negative/eps/2.0/n2.0eps" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"  -dielectric_eps=2.0 -field_amplitude=2.5E+06

mkdir -p ./output/negative/eps/3.0
./streamer positive_2d_left.cfg -output%name="output/negative/eps/3.0/n3.0eps" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"  -dielectric_eps=3.0 -field_amplitude=2.5E+06

mkdir -p ./output/negative/eps/5.0
./streamer positive_2d_left.cfg -output%name="output/negative/eps/5.0/n5.0eps" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95"  -dielectric_eps=5.0 -field_amplitude=2.5E+06
