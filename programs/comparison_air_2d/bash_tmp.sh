#./streamer standard.cfg -output%name="output/field_table" -field_table="voltage_time.txt"

#./streamer standard.cfg -output%name="output/no_attach"

#./streamer standard.cfg -output%name="output/transport/designed" -input_data%file="designed_transport.txt" -input_data%bulk_transport=t -input_data%bulk_scale_reactions=t

#./streamer standard.cfg -output%name="output/infinite/infinite_65ns"

#./streamer standard.cfg -output%name="output/voltage_time_2" -field_table="voltage_time_2.txt"

./streamer standard.cfg seed_rod.cfg -use_electrode=f -output%name="output/seed_electrode/seed_rod" -use_end_streamer_length=f
