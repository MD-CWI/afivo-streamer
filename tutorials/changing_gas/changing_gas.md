# Tutorial: changing the gas and the chemistry {#tutorial-changing-gas}

[TOC]

Using a different gas requires the following steps:

1. If no suitable transport data file is already avaible, electron-neutral cross sections for the gas have to be obtained and electron transport data has to be computed, see [Generating electron transport data](documentation/transport_data.md)
2. Depending on the application, additional chemical reactions have to be obtained.
3. If applicable, photoionization parameters for the gas have to be obtained. Unfortunately, there is often a lack of information on these parameters.
4. Updating the config file

Below, we will provide a basic outline of these steps for a somewhat case of simulations in pure argon, without any additional chemistry or photoionization. **Note**: all the referenced files can be found in the folder `tutorials/changing_gas`

## Cross sections and transport data

We will skip a detailed explanation here, since this step is covered at @ref td-boltzmann. Here it is assumed that the Linux version of `Bolsig+` is used, which is called `bolsigminus`. However, the output can be obtained with the Windows version, or the online version available at [lxcat](https://lxcat.net). Coupled with the program is a set of Siglo cross sections that we will use, which are also available from [lxcat](https://lxcat.net). We will use the Bolsig input file `bolsigminus_argon_input.txt` available at `tutorials/changing_gas`.

## Converting to a format that the code can read

After executing Bolsig, a file called `argon_swarm.dat` is produced. This file can be converted using `tools/bolsig_convert.py`:

    tools/bolsig_convert.py <path>/argon_swarm.dat <path>/argon_transport_data.txt

which produces a file `argon_transport_data.txt`. Sometimes an error like the following can occur `ValueError: could not convert string '-0.2757-206' to float64`. These wrongly formatted numbers should then be corrected in the Bolsig output.

We will now add a minimal chemistry, which will only consist of an ionization reaction:

    reaction_list
    -----------------------
    # Only a single ionization reaction
    e + Ar -> e + e + Ar+,field_table,C3 Ar Ionization 15.80 eV
    -----------------------

The resulting file is available as `argon_chemistry.txt`.

## Updating the config file

Note that we will not include photoionization here, so that the following parameters should be changed:

    # Whether photoionization is enabled:
    photoi%enabled = F

    # Gas component names:
    gas%components = Ar

    # Gas component fractions:
    gas%fractions = 1.0

    # Input file with transport (and reaction) data:
    input_data%file = <path>/argon_chemistry.txt

An example of config file for axisymmetric simulations is available at `tutorial_argon.cfg`.
