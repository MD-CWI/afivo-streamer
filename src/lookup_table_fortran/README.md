# Lookup table

A [lookup table](https://en.wikipedia.org/wiki/Lookup_table)
is an array that replaces runtime computation with a simpler array indexing operation.
The savings in terms of processing time can be significant, since retrieving
a value from memory is often faster than undergoing an "expensive" computation
or input/output operation.

## The lookup table type

The package `lookup_table` can easily be used in program
environments where the use of lookup tables makes sense.
The tables may be precalculated and stored in a
`lookup_table_t`
type defined as

    !> The lookup table type. There can be one or more columns, for which values
    !> can be looked up for a given 'x-coordinate'.
    type lookup_table_t
       private
       integer  :: n_rows !< The number of rows
       integer  :: n_cols !< The number of columns
       real(dp) :: x_min  !< The minimum lookup coordinate
       real(dp) :: dx     !< The x-spacing in the lookup coordinate
       real(dp) :: inv_dx !< The inverse x-spacing

       ! The table is stored in two ways, to speed up different types of lookups.
       real(dp), allocatable :: cols_rows(:, :) !< The table in column-major order
       real(dp), allocatable :: rows_cols(:, :) !< The table in row-major order
    end type lookup_table_t

## Finding indices in a table

Besides the type `lookup_table_t`
we need a way to find locations in the table,
which can be used to speed up multiple lookups of different columns.
Therefore, another type `lt_loc_t`
is defined:

    type LT_loc_t
       private
       integer  :: low_ix   !< The x-value lies between low_ix and low_ix+1
       real(dp) :: low_frac !< The distance from low_ix (up to low_ix+1), given
                            !< as a real number between 0 and 1.
    end type LT_loc_t

The values in the table are calculated once only as part of a program's initialization phase. 
Further, we need some subroutines for finding special indices in lookup tables.
Such subroutines can be found in module
`m_find_index.f90`

    !> Linear search of sorted list for the smallest ix such that list(ix) >= val.
    !> On failure, returns size(list)+1
    function find_index_linear(list, val) result(ix)

    !> Binary search of sorted list for the smallest ix such that list(ix) >= val.
    !> On failure, returns size(list)+1
    function find_index_bsearch(list, val) result(ix)

    !> Adaptive search (combination of linear and binary search) of sorted list
    !> for the smallest ix such that list(ix) >= val. On failure, returns
    !> size(list)+1
    function find_index_adaptive(list, val) result(ix)

## Public methods

In order to optimize the use of the lookup tables, the following public
functions and subroutines are presented in module `m_lookup_table`:

    ! Public methods
    public :: LT_create           ! Create a new lookup table
    public :: LT_get_xdata        ! Get the x-values of a table
    public :: LT_get_spaced_data  ! Convert values to regularly spaced
    public :: LT_set_col          ! Set one table column
    public :: LT_add_to_col       ! Add to a table column
    public :: LT_add_col          ! Add a column
    public :: LT_get_loc          ! Get the index (row) of a value
    public :: LT_get_col          ! Interpolate one column
    public :: LT_get_mcol         ! Interpolate multiple columns
    public :: LT_get_col_at_loc   ! Get one column at location
    public :: LT_get_mcol_at_loc  ! Get multiple columns at location
    public :: LT_get_num_rows
    public :: LT_get_num_cols
    public :: LT_get_data         ! Get all the data of the table
    public :: LT_lin_interp_list  ! Linearly interpolate a list
    public :: LT_to_file          ! Store lookup table in file
    public :: LT_from_file        ! Restore lookup table from file

## Tests

An example how to add data is shown in the test programs
`test_lookup_table_performance`, `test_find_index_performance`, and
`test_find_index_simple`.

