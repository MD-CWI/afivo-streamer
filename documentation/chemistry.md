# Chemistry

# Including chemical reactions

Chemical reactions can be defined in the `input_data%file`, as in the following example:

    reaction_list
    -----------------------
    # Ionization
    e + N2 -> e + e + N2+,field_table,C25 N2 Ionization 15.60 eV
    e + O2 -> e + e + O2+,field_table,C43 O2 Ionization 12.06 eV
    # Attachment
    e + O2 + O2 -> O2-,field_table,C27 O2 Attachment
    e + O2 -> O- + O,field_table,C28 O2 Attachment
    # Detachment
    O2- + M -> e + O2 + M,exp_v1,1.24e-17 179.0 8.8
    O- + N2 -> e + N2O,exp_v1,1.16e-18 48.9 11.0
    # Negative ion conversion
    O- + O2 -> O2- + O,exp_v1,6.96e-17 198.0 5.6
    O- + O2 + M -> O3- + M,exp_v2,1.1e-42 65.0
    O3- + O -> O2- + O2,constant,3.2e-16
    # Positive ion conversion
    # We can assume this rates are constant
    N2+ + N2 + M -> N4+ + M,constant,5.0e-41
    N4+ + O2 -> N2 + N2 + O2+,constant,2.5e-16
    O2+ + O2 + M -> O4+ + M,constant,2.4e-42
    -----------------------

The format of these reactions is

    reaction,rate_type,value(s) [, length unit]

where:

* `reaction` is the reaction text, as in the example above
* `rate_type` denotes the method used to obtain the reaction rate. Several
  options are described below. For a list of all supported rate types, see @ref
  m_chemistry::get_rates()
* `value(s)` are one or more values needed to obtain the reaction rate, for
  example the name of a table or the coefficients of a function
* `length unit` optionally denotes the length unit used in the reaction rate. It
  can be `m` (the default) or `cm`. This setting can be used to convert reaction
  rates given in units of `cm^3/s` or `cm^6/s` to `m^3/s` or `m^6/s` respectively.

## field_table

For a table of the reaction rate versus the reduced electric field (E/N) in
Townsend. Value: the name of the table in `input_data%file`

## constant

For a constant reaction rate. Value: the reaction rate.

## `linear`

For a reaction rate linear in E/N. Values: c1, c2. The rate is given by
`c1 * (Td - c2)` where `Td` is the field in Townsend

## exp_v1

For an exponential dependence on E/N. Values: c1, c2, c3. The rate is given by `c(1) * exp(-(c(2) / (c(3) + fields))**2)`

## k1_func

Values: c1. The rate is given by `c(1) * (300 / Te)**2`, with `Te` the electron temperature.

## k2_func

For a constant reaction rate. Is the same as `constant`, but was implemented seperately by the chemistry data repo. Values: c1. The rate is given by `c(1)`

## k3_func

Values: c1, c2, c3, c4. The rate is given by `(c(1) * (kB * Te + c(2))**2 - c(3)) * c(4)`, with `kB` the Boltzmann constant and `Te` the electron temperature.

## k4_func

Values: c1, c2, c3. The rate is given by `c(1) * (T / 300)**c(2) * exp(-c(3) / T)`, with `T` the gas temperature.

## k5_func

Values: c1, c2. The rate is given by `c(1) * exp(-c(2) / T)`, with `T` the gas temperature.

## k6_func

Values: c1, c2. The rate is given by `c(1) * T**c(2)`, with `T` the gas temperature.

## k7_func

Values: c1, c2, c3. The rate is given by `c(1) * (T / c(2))**c(3)`, with `T` the gas temperature.

## k8_func

Values: c1, c2. The rate is given by `c(1) * (300 / T)**c(2)`, with `T` the gas temperature.

## k9_func

Values: c1, c2. The rate is given by `c(1) * exp(-c(2) * T)`, with `T` the gas temperature.

## k10_func

Values: c1, c2. The rate is given by `10**(c(1) + c(2) * (T - 300))`, with `T` the gas temperature.

## k11_func

Values: c1, c2, c3. The rate is given by `c(1) * (300 / T)**c(2) * exp(-c(3) / T)`, with `T` the gas temperature.

## k12_func

Values: c1, c2, c3. The rate is given by `c(1) * T**c(2) * exp(-c(3) / T)`, with `T` the gas temperature.
