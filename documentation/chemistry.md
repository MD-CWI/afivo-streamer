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

    reaction, rate_type, value

Some of the `rate_type` options are described below. For a list of all supported
rate types, see @ref m_chemistry::get_rates()

## field_table

For a table of the reaction rate versus the reduced electric field (E/N) in
Townsend. Options: the name of the table in `input_data%file`

## constant

For a constant reaction rate. Options: the reaction rate.

## `linear`

For a reaction rate linear in E/N. Options: c1, c2. The rate is given by
`c1 * (Td - c2)` where `Td` is the field in Townsend

## exp_v1

For an exponential dependence on E/N. Options: c1, c2, c3. The rate is given by `c(1) * exp(-(c(2) / (c(3) + fields))**2)`