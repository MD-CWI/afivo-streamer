# mkdir -p ./output/final/standard0.26
# ./streamer positive_2d_left.cfg -output%name="output/final/standard0.26/standard0.26"

mkdir -p ./output/final/ions/ions0.1
./streamer positive_2d_left.cfg -dielectric%gamma_se_ion=0.1 -output%name="output/final/ions/ions0.1/0.1sion"

mkdir -p ./output/final/ions/ions0.3
./streamer positive_2d_left.cfg -dielectric%gamma_se_ion=0.3 -output%name="output/final/ions/ions0.3/0.3sion"

mkdir -p ./output/final/ions/ions0.5
./streamer positive_2d_left.cfg -dielectric%gamma_se_ion=0.5 -output%name="output/final/ions/ions0.5/0.5sion"
