# Specifying voltage {#voltage-intro}
 - There is now the option to specify whether we provide the voltage (in Volt) or electric field (in V/m). The syntax is as follows:
 	`field_given_by= voltage 34`

	`field_given_by= field 0.25e7` (Inside the code, this field value is converted into the correspoding voltage and then used for computation)

	`field_given_by= voltage_table <file_name>.txt`

	`field_given_by=  field_table <file_name>.txt`

- For the sake of backward compatiblity (with earlier versions of the code), the field (in V/m) can still be specified by `field_amplitude`.
# Syntax for applying a voltage pulse {#voltage-pulse}
The code has the functionality to generate trapezoid-shaped pulses. There is functionality to specify the pulse rise/fall time, constant value time, number of pulses (default=1), time for one complete pulse. Using the user provided input, the code computes the time for which the voltage/field is turned off. 
- `field_pulse_width= <time in seconds>`. This is the pulse width excluding the pulse rise, and fall time. By default, we have a pulse width of infinity, i.e., the voltage/field value is constantly applied throughout the duration of the simulation.

- `field_num_pulses= <integer number of pulses>`. This is the number of pulses that we want. By default, this variable is 1.

- `field_rise_time= <linear rise and fall times of the pulse (seconds)>`. By default, this value is 0 sec. 

- `field_pulse_period= <total duration of the pulse, including the off-time in sec>`. Make sure that this value is larger than the pulse-on time (`field_pulse_width + 2*field_rise_time`), else the code will exit with an error message. 

- Please look at the .cfg file `programs/standard_2d/streamer_cyl_2pulses.cfg` for an example on how to setup and run a simulation with 2 pules.
