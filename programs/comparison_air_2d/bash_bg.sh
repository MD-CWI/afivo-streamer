#./streamer standard.cfg -background_density=0.0 -output%name="output/no_background"
mkdir output/BG/

./streamer standard.cfg -background_density=1.0E+3 -output%name="output/BG/BG_1e3"

./streamer standard.cfg -background_density=1.0E+13 -output%name="output/BG/BG_1e13"

./streamer standard.cfg -background_density=1.0E+15 -output%name="output/BG/BG_1e15"

./streamer standard.cfg seed.cfg -background_density=0.0 -output%name="output/BG/seed"
