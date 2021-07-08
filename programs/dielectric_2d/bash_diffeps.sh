mkdir -p ./output/eps/2.0
./streamer positive_2d_left.cfg -output%name="output/eps/2.0/2.0eps" -dielectric_eps=2.0

mkdir -p ./output/eps/3.0
./streamer positive_2d_left.cfg -output%name="output/eps/3.0/3.0eps" -dielectric_eps=3.0

mkdir -p ./output/eps/5.0
./streamer positive_2d_left.cfg -output%name="output/eps/5.0/5.0eps" -dielectric_eps=5.0 -end_time=12.0e-9
