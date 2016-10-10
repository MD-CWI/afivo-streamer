!> A Fortran 90 module for creating lookup tables. These tables can be used to
!> efficiently interpolate one or more values.
!>
!> Author: Jannis Teunissen

module m_lookup_table
  implicit none
  private

  ! The precision of the real numbers used in the tables
  integer, parameter :: dp = kind(1.0d0)

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

  !> Type to indicate a location in the lookup table, which can be used to speed
  !> up multiple lookups of different columns.
  type LT_loc_t
     private
     integer  :: low_ix   !< The x-value lies between low_ix and low_ix+1
     real(dp) :: low_frac !< The distance from low_ix (up to low_ix+1), given
                          !< as a real number between 0 and 1.
  end type LT_loc_t

  ! Public types
  public :: lookup_table_t
  public :: LT_loc_t

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

contains

  !> This function returns a new lookup table
  function LT_create(x_min, x_max, n_rows, n_cols) result(my_lt)
    real(dp), intent(in) :: x_min !< Minimum x-coordinate
    real(dp), intent(in) :: x_max !< Maximum x-coordinate
    integer, intent(in)  :: n_rows !< How many x-values to store
    integer, intent(in)  :: n_cols !< Number of variables that will be looked up
    type(lookup_table_t) :: my_lt

    if (x_max <= x_min) print *, "set_xdata: x_max should be > x_min"
    if (n_rows <= 1)    print *, "set_xdata: n_rows should be bigger than 1"

    my_lt%n_rows = n_rows
    my_lt%x_min  = x_min
    my_lt%dx     = (x_max - x_min) / (n_rows - 1)
    my_lt%inv_dx = 1 / my_lt%dx

    allocate(my_lt%cols_rows(n_cols, n_rows))
    allocate(my_lt%rows_cols(n_rows, n_cols))
    my_lt%cols_rows = 0
    my_lt%rows_cols = 0
    my_lt%n_cols    = n_cols
  end function LT_create

  !> Returns the x-coordinates of the lookup table
  function LT_get_xdata(my_lt) result(xdata)
    type(lookup_table_t), intent(in) :: my_lt
    real(dp)                    :: xdata(my_lt%n_rows)
    integer                     :: ix

    do ix = 1, my_lt%n_rows
       xdata(ix) = my_lt%x_min + (ix-1) * my_lt%dx
    end do
  end function LT_get_xdata

  !> Linearly interpolate the (x, y) input data to the new_x coordinates
  function LT_get_spaced_data(in_x, in_y, new_x) result(out_yy)
    real(dp), intent(in) :: in_x(:), in_y(:), new_x(:)
    real(dp)             :: out_yy(size(new_x))
    integer              :: ix
    do ix = 1, size(new_x)
       call LT_lin_interp_list(in_x, in_y, new_x(ix), out_yy(ix))
    end do
  end function LT_get_spaced_data

  !> Fill the column with index col_ix using the linearly interpolated (x, y)
  !> data
  subroutine LT_set_col(my_lt, col_ix, x, y)
    type(lookup_table_t), intent(inout) :: my_lt
    integer, intent(in)                 :: col_ix
    real(dp), intent(in)                :: x(:), y(:)
    my_lt%cols_rows(col_ix, :) = &
         LT_get_spaced_data(x, y, LT_get_xdata(my_lt))
    my_lt%rows_cols(:, col_ix) = my_lt%cols_rows(col_ix, :)
  end subroutine LT_set_col

  !> Add a new column by linearly interpolating the (x, y) data
  subroutine LT_add_col(my_lt, x, y)
    type(lookup_table_t), intent(inout) :: my_lt
    real(dp), intent(in)                :: x(:), y(:)
    type(lookup_table_t)                :: temp_lt

    temp_lt = my_lt
    deallocate(my_lt%cols_rows)
    deallocate(my_lt%rows_cols)
    allocate(my_lt%cols_rows(my_lt%n_cols+1, my_lt%n_rows))
    allocate(my_lt%rows_cols(my_lt%n_rows, my_lt%n_cols+1))
    my_lt%cols_rows(1:my_lt%n_cols, :) = temp_lt%cols_rows
    my_lt%rows_cols(:, 1:my_lt%n_cols) = temp_lt%rows_cols
    my_lt%n_cols                       = my_lt%n_cols + 1
    my_lt%cols_rows(my_lt%n_cols, :)   = LT_get_spaced_data(x, y, LT_get_xdata(my_lt))
    my_lt%rows_cols(:, my_lt%n_cols)   = my_lt%cols_rows(my_lt%n_cols, :)
  end subroutine LT_add_col

  !> Add the (x,y) data to a given column
  subroutine LT_add_to_col(my_lt, col_ix, x, y)
    type(lookup_table_t), intent(inout) :: my_lt
    integer, intent(in)            :: col_ix
    real(dp), intent(in)           :: x(:), y(:)
    my_lt%cols_rows(col_ix, :) = &
         my_lt%cols_rows(col_ix, :) + LT_get_spaced_data(x, y, LT_get_xdata(my_lt))
    my_lt%rows_cols(:, col_ix) = my_lt%cols_rows(col_ix, :)
  end subroutine LT_add_to_col

  !> Get a location in the lookup table
  elemental function LT_get_loc(my_lt, x) result(my_loc)
    type(lookup_table_t), intent(in) :: my_lt
    real(dp), intent(in)             :: x
    type(LT_loc_t)                   :: my_loc
    real(dp)                         :: frac

    frac            = (x - my_lt%x_min) * my_lt%inv_dx
    my_loc%low_ix   = ceiling(frac)
    my_loc%low_frac = my_loc%low_ix - frac

    ! Check bounds
    if (my_loc%low_ix < 1) then
       my_loc%low_ix   = 1
       my_loc%low_frac = 1
    else if (my_loc%low_ix >= my_lt%n_rows) then
       my_loc%low_ix   = my_lt%n_rows - 1
       my_loc%low_frac = 0
    end if
  end function LT_get_loc

  !> Get the values of all columns at x
  function LT_get_mcol(my_lt, x) result(col_values)
    type(lookup_table_t), intent(in) :: my_lt
    real(dp), intent(in)             :: x
    real(dp)                         :: col_values(my_lt%n_cols)
    type(LT_loc_t)                   :: loc

    loc        = LT_get_loc(my_lt, x)
    col_values = LT_get_mcol_at_loc(my_lt, loc)
  end function LT_get_mcol

  !> Get the value of a single column at x
  elemental function LT_get_col(my_lt, col_ix, x) result(col_value)
    type(lookup_table_t), intent(in) :: my_lt
    integer, intent(in)         :: col_ix
    real(dp), intent(in)        :: x
    real(dp)                    :: col_value
    type(LT_loc_t)              :: loc

    loc       = LT_get_loc(my_lt, x)
    col_value = LT_get_col_at_loc(my_lt, col_ix, loc)
  end function LT_get_col

  !> Get the values of all columns at a location
  function LT_get_mcol_at_loc(my_lt, loc) result(col_values)
    type(lookup_table_t), intent(in) :: my_lt
    type(LT_loc_t), intent(in)       :: loc
    real(dp)                         :: col_values(my_lt%n_cols)

    col_values = loc%low_frac * my_lt%cols_rows(:, loc%low_ix) + &
         (1-loc%low_frac) * my_lt%cols_rows(:, loc%low_ix+1)
  end function LT_get_mcol_at_loc

  !> Get the value of a single column at a location
  elemental function LT_get_col_at_loc(my_lt, col_ix, loc) result(col_value)
    type(lookup_table_t), intent(in) :: my_lt
    integer, intent(in)              :: col_ix
    type(LT_loc_t), intent(in)       :: loc
    real(dp)                         :: col_value

    col_value = loc%low_frac * my_lt%rows_cols(loc%low_ix, col_ix) + &
         (1-loc%low_frac) * my_lt%rows_cols(loc%low_ix+1, col_ix)
  end function LT_get_col_at_loc

  !> Return the number of rows
  integer function LT_get_num_rows(my_lt)
    type(lookup_table_t), intent(in) :: my_lt
    LT_get_num_rows = my_lt%n_rows
  end function LT_get_num_rows

  !> Return the number of columns
  integer function LT_get_num_cols(my_lt)
    type(lookup_table_t), intent(in) :: my_lt
    LT_get_num_cols = size(my_lt%cols_rows, 1)
  end function LT_get_num_cols

  !> Get the x-coordinates and the columns of the lookup table
  subroutine LT_get_data(my_lt, x_data, cols_rows)
    type(lookup_table_t), intent(in) :: my_lt
    real(dp), intent(out)            :: x_data(:), cols_rows(:, :)

    x_data    = LT_get_xdata(my_lt)
    cols_rows = my_lt%cols_rows
  end subroutine LT_get_data

  !> Compute by use of linear interpolation the value in the middle
  ! of a domain D = [x_list(1) , x_list(size(x_list))].
  ! If x_value is left of domain  D,
  ! then the value becomes the value at the left side of D,
  ! if x_value is right of domain D,
  ! then the value becomes the value at the rigth side of D
  subroutine LT_lin_interp_list(x_list, y_list, x_value, y_value)
    use m_find_index
    real(dp), intent(in)  :: x_list(:), y_list(:)
    real(dp), intent(in)  :: x_value
    real(dp), intent(out) :: y_value

    integer               :: ix, iMin, iMax
    real(dp)              :: temp

    iMin = 1
    iMax = size(x_list)

    if (x_value <= x_list(iMin)) then
       y_value = y_list(iMin)
    else if (x_value >= x_list(iMax)) then
       y_value = y_list(iMax)
    else
       ix = find_index_adaptive(x_list, x_value)
       temp = (x_value - x_list(ix-1)) / (x_list(ix) - x_list(ix-1))
       y_value = (1 - temp) * y_list(ix-1) + temp * y_list(ix)
    end if
  end subroutine LT_lin_interp_list

  !> Write the lookup table to file (in binary, potentially unportable)
  subroutine LT_to_file(my_lt, filename)
    type(lookup_table_t), intent(in) :: my_lt
    character(len=*), intent(in) :: filename
    integer                      :: my_unit

    open(newunit=my_unit, file=trim(filename), form='UNFORMATTED', &
         access='STREAM', status='REPLACE')
    write(my_unit) my_lt%n_rows, my_lt%n_cols
    write(my_unit) my_lt%x_min, my_lt%dx, my_lt%inv_dx
    write(my_unit) my_lt%cols_rows
    close(my_unit)
  end subroutine LT_to_file

  !> Read the lookup table from file (in binary, potentially unportable)
  subroutine LT_from_file(my_lt, filename)
    type(lookup_table_t), intent(inout) :: my_lt
    character(len=*), intent(in)    :: filename
    integer                         :: my_unit

    open(newunit=my_unit, file=trim(filename), form='UNFORMATTED', &
         access='STREAM', status='OLD')
    read(my_unit) my_lt%n_rows, my_lt%n_cols
    read(my_unit) my_lt%x_min, my_lt%dx, my_lt%inv_dx

    allocate(my_lt%cols_rows(my_lt%n_cols, my_lt%n_rows))
    allocate(my_lt%rows_cols(my_lt%n_rows, my_lt%n_cols))

    read(my_unit) my_lt%cols_rows
    my_lt%rows_cols = transpose(my_lt%cols_rows)

    close(my_unit)
  end subroutine LT_from_file

end module m_lookup_table
