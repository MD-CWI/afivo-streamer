# Chemistry {#chem-head}

[TOC]

# Including chemical reactions {#chem-intro}

Chemical reactions can be defined in the `input_data%%file`, as in the following example:

    reaction_list
    -----------------------
    # Ionization
    e + N2 -> e + e + N2+,field_table,C25 N2 Ionization 15.60 eV
    e + N2 -> e + e + N2+,field_table,C26 N2 Ionization 18.80 eV
    e + O2 -> e + e + O2+,field_table,C43 O2 Ionization 12.06 eV
    # Attachment
    e + O2 + O2 -> O2- + O2,field_table,C27 O2 Attachment
    e + O2 -> O- + O,field_table,C28 O2 Attachment
    # Detachment (Panchesnyi 2013)
    O2- + M -> e + O2 + M,c1*exp(-(c2/(c3+Td))**2),1.24e-17 179.0 8.8
    O- + N2 -> e + N2O,c1*exp(-(c2/(c3+Td))**2),1.16e-18 48.9 11.0
    # Negative ion conversion (Panchesnyi 2013)
    O- + O2 -> O2- + O,c1*exp(-(c2/(c3+Td))**2),6.96e-17 198.0 5.6
    O- + O2 + M -> O3- + M,c1*exp(-(Td/c2)**2),1.1e-42 65.0
    # Positive ion conversion (Kossyi et. al. 1992, Aleksandrov and Bazelyan 1999)
    N2+ + O2 -> O2+ + N2,c1*(300/Tg)**c2,6e-17 0.5
    N2+ + N2 + M -> N4+ + M,c1*(300/Tg)**c2,5e-41 2.0
    N4+ + O2 -> O2+ + 2N2,c1,2.5e-16
    O2+ + O2 + M -> O4+ + M,c1*(300/Tg)**c2,2.4e-42 3.0
    # Electron recombination (Kossyi et. al. 1992)
    e + N4+ -> 2N2,c1*(300/Te)**c2,2e-12 0.5
    e + O4+ -> 2O2,c1*(300/Te)**c2,1.4e-12 0.5
    -----------------------

The format of these reactions is

    reaction,rate_type,value(s) [, length unit]

where:

* `reaction` is the reaction formula, as in the example above
* `rate_type` denotes the method used to obtain the reaction rate. Several
  options are described below.
* `value(s)` are one or more values needed to obtain the reaction rate, for
  example the name of a table or the coefficients of a function
* `length unit` optionally denotes the length unit used in the reaction rate. It
  can be `m` (the default) or `cm`. This setting can be used to convert reaction
  rates given in units of `cm^3/s` or `cm^6/s` to `m^3/s` or `m^6/s`
  respectively. The time unit is always a second.

All the species that are present in the chemistry file will automatically be present in the simulation, so it is not necessary to provide a list of species.

# Chemistry syntax {#chem-syntax}

## Reaction syntax {#chem-syntax-reactions}

* `M` denotes 'any' gas molecule
* `X+` means a positively charged species `X`
* `X-` means a negatively charged species `X`
* `2X` means `X + X` (same for `3X` etc.)

Note that in the output, special characters such as `+` and `-` are converted, because only alphanumeric symbols and `_` can be include in Silo variable names.

## Reaction groups {#chem-syntax-groups}

Sometimes, there are many similar reactions. To write these more compactly, the following syntax is available, somewhat similar to ZDPlaskin:

    e + @x -> e + e + @x+,field_table,@source
    @x = N2,N2,O2
    @source = C25 N2,C26 N2,C43 O2

The symbols with an `@` will be replaced by the respective values specified in the lines below. So for this example, the reaction set would become:

    e + N2 -> e + e + N2+,field_table,C25 N2
    e + N2 -> e + e + N2+,field_table,C26 N2
    e + O2 -> e + e + O2+,field_table,C43 O2

The number of such replacement groups is flexible.

## Supported reaction rate formats {#chem-syntax-rate-formats}

### Field-dependent reactions
If the reaction depends on the reduced electric field, `rate_type` should be set to `field_table`, and `value` should be set to the name of the table in `input_data%file`. For example

    e + N2 -> e + e + N2+,field_table,C25 N2 Ionization 15.60 eV

means that the rate for this reaction will be obtained by interpolating the tabulated data following the string "C25 N2 Ionization 15.60 eV". An example of such tabulated data is given below, where the first column is the field in Townsend:

    C25 N2 Ionization 15.60 eV
    -----------------------
    1.000000000000000000e+00 0.000000000000000000e+00
    1.151000000000000023e+00 0.000000000000000000e+00
    1.326000000000000068e+00 0.000000000000000000e+00
    1.526000000000000023e+00 0.000000000000000000e+00
    1.758000000000000007e+00 0.000000000000000000e+00
    2.024000000000000021e+00 0.000000000000000000e+00
    2.330000000000000071e+00 0.000000000000000000e+00
    ...
    -----------------------

### Functional expressions

Reactions with analytic formulas can be included by giving the form of the analytic expression.
This expression should not contain spaces, and it should be one of the following:

* `c1`
* `c1*(Td-c2)`
* `c1*exp(-(c2/(c3+Td))**2)`
* `c1*exp(-(Td/c2)**2)`
* `c1*(300/Te)**c2`
* `(c1*(kB_eV*Te+c2)**2-c3)*c4`
* `c1*(Tg/300)**c2*exp(-c3/Tg)`
* `c1*exp(-c2/Tg)`
* `c1*Tg**c2`
* `c1*(Tg/c2)**c3`
* `c1*(300/Tg)**c2`
* `c1*exp(-c2*Tg)`
* `10**(c1+c2*(Tg-300))`
* `c1*(300/Tg)**c2*exp(-c3/Tg)`
* `c1*Tg**c2*exp(-c3/Tg)`
* `c1*exp(-(c2/(c3+Td))**c4)`
* `c1*exp(-(Td/c2)**c3)`
* `c1*exp(-(c2/(kb*(Tg+Td/c3)))**c4)`

For these expressions, the values specified should be `c1`, `c2`, etc. So for example

    O- + O2 + M -> O3- + M,c1*exp(-(Td/c2)**2),1.1e-42 65.0

means that the reaction rate is given by `1.1e-42 * exp(-(Td/65.0)**2)`. The above symbols have the following meanings:

symbol | meaning | unit
---|---|---
Td | The reduced electric field E/N | Townsend
Te | Electron 'temperature' (given by 2*energy/(3*kB))| K
Ti | Ion temperature | K
Tg | Gas temperature | K
kB | Boltzmann constant | J/K
kB_eV | Boltzmann constant | eV/K

# Ignoring the production of certain species {#chem-ignored}

To simplify a reaction set, it is possible to ignore the production of certain species. This can be done by adding a table of the following form:

    ignored_species
    -----------------------
    NO2
    NO
    O
    N
    -----------------------

When ignored species occur on the left-hand side of a reaction, and if they do not correspond to the background gas that is initially present, the reaction will be ignored.
When ignored species occur on the right-hand side of a reaction, their production will be ignored, but the reaction will still be included.

# Adding new types of reactions {#chem-new-reactions}

First check if the new format can be computed according to one of the existing
expressions, for example with different parameters (e.g., you could one
parameter less, or allow for a minus sign). If this is the case, you only have
to make modifications in `m_chemistry::read_reactions()`.

1. Add a new `case` statement with the new reaction format
2. Determine which of the existing rate functions you can re-use
3. Modify the coefficients according to the new format

If instead you have to add a completely new rate function, this can be done as
follows:

1. Define a new parameter, e.g. `integer, parameter :: rate_analytic_kN = ...` in `m_chemistry`
2. Add a new `case` statement with the new reaction format in `m_chemistry::read_reactions()`
3. Implement the new reaction with a new case in `m_chemistry::get_rates()`

# Chemistry output and visualization {#chem-vis}

Each simulation produces the following files:
1. `<base_name>/_rates.txt`: Contains volume and time-integrated reaction rates, for each timestep.
2. `<base_name>/_amounts.txt`: Contains the amount of each species (volume-integrated densities), for each timestep.
3. `<base_name>/_reactions.txt`: Contains the list of all the reactions used in the simulation
4. `<base_name>/_species.txt`: Contains the list of all the species used in the simulation
4. `<base_name>/_stoich_matrix.txt`: This file contains the [stoichiometric matrix](https://en.wikipedia.org/wiki/Stoichiometry#Stoichiometry_matrix) of the reaction set used in the simulation.

The data in the above files can be visualized using the `chemistry_visualize_rates.py` script in the `/tools` directory. The documentation for each of the arguments can be found by typing `python chemistry_visualize_rates.py -h`.


