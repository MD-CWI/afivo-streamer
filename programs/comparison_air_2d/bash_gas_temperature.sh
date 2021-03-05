#mkdir output/gas_temperature

./streamer standard.cfg -gas%temperature=290.0 -output%name="output/gas_temperature/290K"

./streamer standard.cfg -gas%temperature=310.0 -output%name="output/gas_temperature/310K"

./streamer standard.cfg -gas%temperature=360.0 -output%name="output/gas_temperature/360K"
