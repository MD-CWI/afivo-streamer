\section sect_lt Lookup table(LT) and finding indices(FI) in lookup tables

A lookup table (see <a href="https://en.wikipedia.org/wiki/Lookup_table">Wikipedia</a>)
is an array that replaces runtime computation with a simpler array indexing operation.
The savings in terms of processing time can be significant, since retrieving
a value from memory is often faster than undergoing an "expensive" computation
or input/output operation.

\subsection sect_alone Stand alone package

The package <b><code>lookup_table</code></b> can easily be used in program
environments where the use of lookup tables makes sense.
It has been developed as a by-product of <b><code>a5_streamer</code></b>,
a package to simulate streamers. However, the program shown here is completely
separate from <b><code>a5_streamer</code></b>.

The tables may be precalculated and stored in a
<a class="el" href="structm__lookup__table_1_1lookup__table__t.html">lookup_table_t</a>
type defined as
\verbatim
  !> The lookup table type. There can be one or more columns, for which values
  !> can be looked up for a given 'x-coordinate'.
  type lookup_table_t
     private
     integer  :: n_rows !< The number of rows
     integer  :: n_cols !< The number of columns
     real(dp) :: x_min  !< The minimum lookup coordinate
     real(dp) :: dx     !< The x-spacing in the lookup coordinate
     real(dp) :: inv_dx !< The inverse x-spacing

     !The table is stored in two ways, to speed up different types of lookups.
     real(dp), allocatable :: cols_rows(:, :) !< The table in column-major order
     real(dp), allocatable :: rows_cols(:, :) !< The table in row-major order
  end type lookup_table_t
\endverbatim

Besides the type 
<a class="el" href="structm__lookup__table_1_1lookup__table__t.html">lookup_table_t</a>
we need a way to find locations in the table,
which can be used to speed up multiple lookups of different columns.
Therefore, another type
<a class="el" href="structm__lookup__table_1_1lt__loc__t.html">lt_loc_t</a>
is defined: 
\verbatim
  type LT_loc_t
     private
     integer  :: low_ix   !< The x-value lies between low_ix and low_ix+1
     real(dp) :: low_frac !< The distance from low_ix (up to low_ix+1), given
                          !< as a real number between 0 and 1.
  end type LT_loc_t
\endverbatim

The values in the table are calculated once only as part of a program's initialization phase. 
Further, we need some subroutines for finding special indices in lookup tables.
Such subroutines can be found in module
<a class="el" href="m__find__index_8f90_source.html">m_find_index.f90</a>:

\verbatim
  ! Searches within an increasing sorted list for the smallest ix such that
  ! rvalue <= list(ix)
  function FI_linear_r(list, rvalue) result(ix)
    real(dp), intent(in) :: list(:)    !< Increasing sorted list
    real(dp), intent(in) :: rvalue     !< Given real value
    integer              :: ix

  ! Searches within a monotonically increasing sorted list for the ix which is
  !  'closest' to rvalue. If rvalue is not in the list list(ix) will be the
  !  first value larger than rvalue.
  function FI_binary_search_r(list, rvalue) result(ix)
    real(dp), intent(in) :: list(:)    !< Increasing sorted list
    real(dp), intent(in) :: rvalue     !< Given real value
    integer              :: ix

  ! Searches monotonically increasing sorted list for the smallest ix such that
  ! rvalue <= list(ix).
  ! Switches between linear and binary search
  function FI_adaptive_r(list, rvalue) result(ix)
    real(dp), intent(in) :: list(:)    !< Increasing sorted list
    real(dp), intent(in) :: rvalue     !< Given real value
    integer              :: ix
\endverbatim

In order to optimize the use of the lookup tables, the following public
functions and subroutines are presented in module
<a class="el" href="namespacem__lookup__table.html">m_lookup_table</a>:
\verbatim
  ! Returns a new lookup table
  function LT_create(x_min, x_max, n_rows, n_cols) result(my_lt)
    real(dp), intent(in) :: x_min  !< Minimum x-coordinate
    real(dp), intent(in) :: x_max  !< Maximum x-coordinate
    integer, intent(in)  :: n_rows !< How many x-values to store
    integer, intent(in)  :: n_cols !< Number of variables that will be looked up
    type(lookup_table_t) :: my_lt  !< My lookup table

  ! Returns the x-coordinates of the lookup table
  function LT_get_xdata(my_lt) result(xdata)
    type(lookup_table_t), intent(in) :: my_lt    !< My lookup table
    real(dp) :: xdata(my_lt%n_rows)

  ! Linearly interpolate the (x, y) input data to the new_x coordinates
  function LT_get_xdata(my_lt) result(xdata)
    type(lookup_table_t), intent(in) :: my_lt    !< My lookup table

  ! Fill the column with index col_ix using the linearly interpolated (x, y) data
  subroutine LT_set_col(my_lt, col_ix, x, y)
    type(lookup_table_t), intent(inout) :: my_lt  !< My lookup table
    integer, intent(in)                 :: col_ix !< Column of x
    real(dp), intent(in)                :: x(:)   !< x value
    real(dp), intent(in)                :: y(:)   !< y value

  ! Add a new column by linearly interpolating the (x, y) data
  subroutine LT_add_col(my_lt, x, y)
    type(lookup_table_t), intent(inout) :: my_lt   !< My lookup table
    real(dp), intent(in)                :: x(:)    !< x data
    real(dp), intent(in)                :: y(:)    !< y data

  ! Add the (x,y) data to a given column
  subroutine LT_add_to_col(my_lt, col_ix, x, y)
    type(lookup_table_t), intent(inout) :: my_lt  !< My lookup table
    integer, intent(in)                 :: col_ix !< Column of x
    real(dp), intent(in)                :: x(:)   !< x value
    real(dp), intent(in)                :: y(:)   !< y value

  ! Get a location in the lookup table
  elemental function LT_get_loc(my_lt, x) result(my_loc)
    type(lookup_table_t), intent(in) :: my_lt  !< My lookup table
    real(dp), intent(in)             :: x      !< x value
    type(LT_loc_t)                   :: my_loc

  ! Get the values of all columns at x
  function LT_get_mcol(my_lt, x) result(col_values)
    type(lookup_table_t), intent(in) :: my_lt  !< My lookup table
    real(dp), intent(in)             :: x      !< x value


  ! Get the value of a single column at x
  elemental function LT_get_col(my_lt, col_ix, x) result(col_value)
    type(lookup_table_t), intent(in) :: my_lt  !< My lookup table
    integer, intent(in)         :: col_ix      !< Column of x
    real(dp), intent(in)        :: x           !< x value
    real(dp)                    :: col_values(my_lt%n_cols)


  ! Get the values of all columns at a location
  function LT_get_mcol_at_loc(my_lt, loc) result(col_values)
    type(lookup_table_t), intent(in) :: my_lt  !< My lookup table
    type(LT_loc_t), intent(in)       :: loc    !< Location in the lookup table
    real(dp)                    :: col_values(my_lt%n_cols)

  ! Get the value of a single column at a location
  elemental function LT_get_col_at_loc(my_lt, col_ix, loc) result(col_value)
    type(lookup_table_t), intent(in) :: my_lt  !< My lookup table
    integer, intent(in)              :: col_ix !< Column of x
    type(LT_loc_t), intent(in)       :: loc    !< Location in the lookup table
    real(dp)                         :: col_values(my_lt%n_cols)

  ! Return the number of rows
  integer function LT_get_num_rows(my_lt)
    type(lookup_table_t), intent(in) :: my_lt  !< My lookup table

  ! Return the number of columns
  integer function LT_get_num_cols(my_lt)
    type(lookup_table_t), intent(in) :: my_lt  !< My lookup table

  ! Get the x-coordinates and the columns of the lookup table
  subroutine LT_get_data(my_lt, x_data, cols_rows)
    type(lookup_table_t), intent(in) :: my_lt           !< My lookup table
    real(dp), intent(out)            :: x_data(:)       !< x-coordinates
    real(dp), intent(out)            :: cols_rows(:, :) !< Colums in the lookup 
                                                        !  table
  ! Compute by use of linear interpolation the value in the middle of
  ! domain D = [x_list(1) , x_list(size(x_list))].
  ! If x_value is left of domain  D, then the value becomes the value
  ! at the left side of D,
  ! if x_value is right of domain D, then the value becomes the value
  ! at the right side of D
  subroutine LT_lin_interp_list(x_list, y_list, x_value, y_value)
    real(dp), intent(in)  :: x_list(:)   !< List with x values
    real(dp), intent(in)  :: y_list(:)   !< List with y values
    real(dp), intent(in)  :: x_value     !< A given x value
    real(dp), intent(out) :: y_value     !< y value

  ! Write the lookup table to file (in binary, potentially unportable)
  subroutine LT_to_file(my_lt, filename)
    type(lookup_table_t), intent(in) :: my_lt    !< My lookup table
    character(len=*), intent(in)     :: filename !< Filename for lookup table

  ! Read the lookup table from file (in binary, potentially unportable)
  subroutine LT_from_file(my_lt, filename)
    type(lookup_table_t), intent(inout) :: my_lt    !< My lookup table
    character(len=*), intent(in)        :: filename !< Filename for lookup table
\endverbatim

The intended usage is as follows:

    You create your lookup table and add the data to the table

An example how to add data is shown in test program
<a class="el" href="test__lookup__table_8f90_source.html">test_lookup_table</a>.
The result of the first part of this program
is:
\verbatim
 start test_lookup_table
 Testing m_lookup_table.f90 implementation...
 
  *** Vectorized calls (one call for all the data) ***
 Iteration    table_size   max. difference
           1          20   3.3582003742245958E-002
           2          40   8.2030184573754772E-003
           3          80   2.0020795099551236E-003
           4         160   4.9426083187298353E-004
           5         320   1.2283397795254114E-004
           6         640   3.0612901792514968E-005
           7        1280   7.6413075068559877E-006
           8        2560   1.9088447282822329E-006
           9        5120   4.7703337102689147E-007
          10       10240   1.1923917775380488E-007
          11       20480   2.9817307956037098E-008
          12       40960   7.4591307575033738E-009
          13       81920   1.8682921876234104E-009
          14      163840   4.7167403316450418E-010
          15      327680   1.2888212719275316E-010
          16      655360   3.8146930059212991E-011
 You should see 2nd order convergence
 Number of lookups performed :            17000000
 Number of lookups / second  :        163011688.90
 Nanosecond per lookup       :               6.135
\endverbatim
The result of the second part of 
<a class="el" href="test__lookup__table_8f90_source.html">test_lookup_table</a> is:


\verbatim
 Testing multiple column mode (2 columns)
  *** Vectorized calls (one call for all the data) ***
 Iteration    table_size   max. difference
           1          20   3.3582003742245958E-002
           2          40   8.2030184573754772E-003
           3          80   2.0020795099551236E-003
           4         160   4.9426083187298353E-004
           5         320   1.2283397795254114E-004
           6         640   3.0612901792514968E-005
           7        1280   7.6413075068559877E-006
           8        2560   1.9088447282822329E-006
           9        5120   4.7703337102689147E-007
          10       10240   1.1923917775380488E-007
          11       20480   2.9817307956037098E-008
          12       40960   7.4591307575033738E-009
          13       81920   1.8682921876234104E-009
          14      163840   4.7167403316450418E-010
          15      327680   1.2888212719275316E-010
          16      655360   3.8146930059212991E-011
 You should see 2nd order convergence
 Number of lookups performed :            17000000
 Number of lookups / second  :         74655705.45
 Nanosecond per lookup       :              13.395
\endverbatim


\subsection sect_example More examples

If you are interested on how this package is used by  <code>a5_streamer</code>
to handle with lookup tables with transport data, see subroutine 
<a class="el" href="namespacem__streamer.html#ac543e682ffced5108a9e5f33b7c6c1ba">st_load_transport_data</a>

