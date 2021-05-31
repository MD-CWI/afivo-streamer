mkdir -p ./output/um/1e-4
./streamer positive_2d_left.cfg -output%name="output/um/1e-4/1e-4mion" -input_data%ion_mobilities=1.0000E-04 -user_point%write=t

mkdir -p ./output/um/1e-3
./streamer positive_2d_left.cfg -output%name="output/um/1e-3/1e-3mion" -input_data%ion_mobilities=1.0000E-03 -user_point%write=t

mkdir -p ./output/um/5e-4
./streamer positive_2d_left.cfg -output%name="output/um/5e-4/5e-4mion" -input_data%ion_mobilities=5.0000E-04 -user_point%write=t

mkdir -p ./output/um/0
./streamer positive_2d_left.cfg -output%name="output/um/0/0mion" -input_data%ion_mobilities=0 -user_point%write=t
