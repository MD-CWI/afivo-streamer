# Lookup table

A [lookup table](https://en.wikipedia.org/wiki/Lookup_table) is an array that
can replace an expensive computation with a linear interpolation. A good way to
understand the usage of this module is to run:

    make
    ./usage_example

After which you can inspect the file [usage_example.f90](usage_example.f90) to
see how you can incorporate the module in your own programs.

## Requirements

A modern Fortran compiler, such as `gfortran 4.8` or newer.

## Examples

* [usage_example.f90](usage_example.f90) Shows several examples of how values
  can be looked up, a good starting point.
* [test_lookup_table_performance.f90](test_lookup_table_performance.f90)
  Compares the performance of different types of lookups.
* [test_lookup_table_index.f90](test_lookup_table_index.f90)
  Tests whether the correct indices in the table are retrieved
* [test_lookup_table_2d.f90](test_lookup_table_2d.f90) Simple example of a 2D lookup table
* [test_find_index_performance.f90](test_find_index_performance.f90) Compares the performance of linear search,
  binary search and automatically switching between the two.

## Creating a lookup table

Lookup tables are created with a call to the function `LT_create`:

    lt = LT_create(x_min, x_max, n_rows, n_cols, xspacing)

Here `n_rows` controls how many `x`-values are stored. The x-spacing of the
lookup table can be set to:

* `LT_xspacing_linear` (the default, in which the spacing is given by `dx = (x_max - x_min) / (n_rows - 1)`)
* `LT_xspacing_quadratic` for quadratic spacing
* `LT_xspacing_cubic` for cubic spacing, which is a bit slower

The argument `n_cols` specifies how many different columns you want to store,
but it also possible to add extra columns later with `LT_add_col`.

The columns of a table are set using `LT_set_col`:

    call LT_set_col(lt, col_ix, x, y)

where `x` and `y` are vectors of the same length, with the only requirement that
the `x` array is sorted from low to high. These `y` values will be linearly
interpolated and stored in the table. Alternatively, it is possible to directly specify the y-values in a column:

    call LT_set_col_data(lt, col_ix, y)

## Getting values from a table

There are multiple ways to get values from a lookup table:

* You can directly look up a single column at some `x`-coordinate, using
  `LT_get_col`
* You can get all the columns at some `x`-coordinate, using
  `LT_get_mcol`
* You can also get a location in the table with a `LT_loc_t` object, and then do
  lookups using that object. This can be used to speed up multiple lookups of
  different columns

The type `lt_loc_t` is defined as:

    type LT_loc_t
       integer  :: low_ix   !< The x-value lies between low_ix and low_ix+1
       real(dp) :: low_frac !< The distance from low_ix (up to low_ix+1), given
                            !< as a real number between 0 and 1.
    end type LT_loc_t

## Contributors

* [Jannis Teunissen](http://teunissen.net/)
* Margreet Nool (Documentation)
