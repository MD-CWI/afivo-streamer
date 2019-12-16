# #planewrite on
# 
# mkdir -p ./output/final/mion/1e-5
# ./streamer positive_2d_left.cfg -output%name="output/final/mion/1e-5/1e-5mion" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95" -input_data%ion_mobilities=1.0000E-05 -end_time=15e-9
# 
mkdir -p ./output/final/mion/1e-4
./streamer positive_2d_left.cfg -output%name="output/final/mion/1e-4/1e-4mion" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95" -input_data%ion_mobilities=1.0000E-04 -end_time=15e-9

mkdir -p ./output/final/mion/1e-3
./streamer positive_2d_left.cfg -output%name="output/final/mion/1e-3/1e-3mion" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95" -input_data%ion_mobilities=1.0000E-03 -end_time=15e-9

mkdir -p ./output/final/mion/5e-4
./streamer positive_2d_left.cfg -output%name="output/final/mion/5e-4/5e-4mion" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95" -input_data%ion_mobilities=5.0000E-04 -end_time=15e-9

mkdir -p ./output/final/mion/0
./streamer positive_2d_left.cfg -output%name="output/final/mion/0/0mion" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95" -input_data%ion_mobilities=0 -end_time=15e-9

