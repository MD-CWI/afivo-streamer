# mkdir -p ./output/final/gas/O20.05
# ./streamer positive_2d_left.cfg -output%name="output/final/gas/O20.05/0.05O2" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95" -gas%fractions="0.95 0.05"
# 
# mkdir -p ./output/final/gas/O20.4
# ./streamer positive_2d_left.cfg -output%name="output/final/gas/O20.4/0.4O2" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95" -gas%fractions="0.6 0.4"
# 
# mkdir -p ./output/final/gas/O20.8
# ./streamer positive_2d_left.cfg -output%name="output/final/gas/O20.8/0.8O2" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95" -gas%fractions="0.2 0.8"

# mkdir -p ./output/final/gas/O21.0
# ./streamer positive_2d_left.cfg -output%name="output/final/gas/O21.0/1.0O2" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95" -gas%fractions="0.0 1.0"

mkdir -p ./output/final/gas/O20.01
./streamer positive_2d_left.cfg -output%name="output/final/gas/O20.01/0.01O2" -seed_rel_r0="0.2625 1.0" -seed_rel_r1="0.2625 0.95" -gas%fractions="0.99 0.01" -end_time=15e-9
