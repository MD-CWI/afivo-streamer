mkdir output/transport/

./streamer standard.cfg -output%name="output/transport/IST_Lisbon" -input_data%file="IST_Lisbon_xiaoran.txt"

./streamer standard.cfg -output%name="output/transport/biagi" -input_data%file="Biagi_xiaoran.txt"

./streamer standard.cfg -output%name="output/transport/new_biagi" -input_data%file="td_biagi_air_xiaoran.txt"

./streamer standard.cfg -output%name="output/transport/new_biagi_bulk" -input_data%file="td_biagi_air_xiaoran.txt" -input_data%bulk_transport=t -input_data%bulk_scale_reactions=t

./streamer standard.cfg -output%name="output/transport/morgan" -input_data%file="morgan_xiaoran.txt"

./streamer standard.cfg -output%name="output/transport/triniti" -input_data%file="triniti_xiaoran.txt"

./streamer standard.cfg -output%name="output/transport/designed" -input_data%file="designed_transport.txt" -input_data%bulk_transport=t -input_data%bulk_scale_reactions=t
