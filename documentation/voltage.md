# Specifying the applied voltage or field {#voltage-description}

[TOC]

# Specifying the voltage {#voltage-intro}

There is now a single option to specify whether we provide the voltage (in Volt) or electric field (in V/m). The syntax is as follows:

    field_given_by = <voltage,field,voltage_table,field_table> <value>

For example:

    # To apply a voltage of 34 Volts
    field_given_by = voltage 34

    # To apply a background field of 0.25e7 V/m, which is internally converted into a voltage
    field_given_by = field 0.25e7

    # To apply a time-dependent voltage given in file_name.txt
    field_given_by = voltage_table file_name.txt

    # To apply a time-dependent background field given in file_name.txt
    field_given_by =  field_table file_name.txt

The voltage or field table should have the following structure:

    voltage_vs_time (or field_vs_time)
    --------
    t0 v0
    t1 v1
    ... ...
    --------

For the sake of backward compatiblity (with earlier versions of the code), the field (in V/m) can still be specified by `field_amplitude`.

# Applying a voltage pulse {#voltage-pulse}

The code has the functionality to generate trapezoid-shaped pulses, by specifying the pulse rise/fall time, width, period and the number of pulses (default=1).

- `field_pulse_width = <time in seconds>`. This is the pulse width excluding the pulse rise and fall time. By default, we have a pulse width of infinity, i.e., the voltage/field value is constantly applied throughout the duration of the simulation.

- `field_num_pulses = <integer number of pulses>`. The number of pulses, by default one.

- `field_rise_time = <linear rise and fall times of the pulse (seconds)>`. By default, this value is 0 sec.

- `field_pulse_period = <total duration of the pulse, including the off-time in sec>`. Make sure that this value is larger than the "pulse on time" (`field_pulse_width + 2*field_rise_time`), else the code will exit with an error message.

Please look at the .cfg file `programs/standard_2d/streamer_cyl_2pulses.cfg` for an example on how to setup and run a simulation with 2 pulses.
