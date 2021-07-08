mkdir -p ./output/field/2.1
./streamer positive_2d_left.cfg -output%name="output/field/2.1/2.1field" -field_amplitude=-2.1E+06

mkdir -p ./output/field/2.3
./streamer positive_2d_left.cfg -output%name="output/field/2.3/2.3field" -field_amplitude=-2.3E+06

mkdir -p ./output/field/2.5
./streamer positive_2d_left.cfg -output%name="output/field/2.5/2.5field" -field_amplitude=-2.5E+06

mkdir -p ./output/field/2.8
./streamer positive_2d_left.cfg -output%name="output/field/2.8/2.8field" -field_amplitude=-2.8E+06
