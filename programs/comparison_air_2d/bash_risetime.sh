mkdir output/voltage

./streamer standard.cfg -field_table="UNDEFINED" -field_amplitude=-1.50E+05 -field_rise_time=0e-9 -output%name="output/voltage/rt_0ns"

./streamer standard.cfg -field_table="UNDEFINED" -field_amplitude=-1.50E+05 -field_rise_time=20e-9 -output%name="output/voltage/rt_20ns"

./streamer standard.cfg -field_table="UNDEFINED" -field_amplitude=-1.50E+05 -field_rise_time=40e-9 -output%name="output/voltage/rt_40ns"

./streamer standard.cfg -field_table="UNDEFINED" -field_amplitude=-1.50E+05 -field_rise_time=65e-9 -output%name="output/voltage/rt_65ns"

./streamer standard.cfg -field_table="UNDEFINED" -field_rise_time=65e-9 -field_amplitude=-1.250E+05 -output%name="output/voltage/sim_12_5kV"

./streamer standard.cfg -field_table="UNDEFINED" -field_rise_time=65e-9 -field_amplitude=-1.750E+05 -output%name="output/voltage/sim_17_5kV"
