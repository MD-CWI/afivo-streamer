module m_lookup_table
  implicit none
  private

  integer, parameter :: dp = kind(1.0d0)

  type LT_table_t
     private
     integer               :: n_rows
     integer               :: n_cols
     real(dp)              :: x_min, inv_dx
     real(dp), allocatable :: cols_rows(:, :)
     real(dp), allocatable :: rows_cols(:, :)
  end type LT_table_t

  type LT_loc_t
     private
     integer  :: low_ix
     real(dp) :: low_frac
  end type LT_loc_t

  ! Public types
  public :: LT_table_t
  public :: LT_loc_t

  ! Public methods
  public :: LT_create
  public :: LT_set_col
  public :: LT_add_to_col
  public :: LT_add_col
  public :: LT_get_loc
  public :: LT_get_col
  public :: LT_get_mcol
  public :: LT_get_col_at_loc
  public :: LT_get_mcol_at_loc
  public :: LT_get_num_rows
  public :: LT_get_num_cols
  public :: LT_get_data
  public :: LT_get_xdata
  public :: LT_lin_interp_list
  public :: LT_to_file
  public :: LT_from_file

contains

  function LT_create(x_min, x_max, n_rows, n_cols) result(my_lt)
    real(dp), intent(in) :: x_min, x_max
    integer, intent(in)  :: n_rows, n_cols
    type(LT_table_t)      :: my_lt

    if (x_max <= x_min) print *, "set_xdata: x_max should be > x_min"
    if (n_rows <= 1)    print *, "set_xdata: n_rows should be bigger than 1"

    my_lt%n_rows = n_rows
    my_lt%x_min  = x_min
    my_lt%inv_dx = (n_rows - 1) / (x_max - x_min)

    allocate(my_lt%cols_rows(n_cols, n_rows))
    allocate(my_lt%rows_cols(n_rows, n_cols))
    my_lt%cols_rows = 0
    my_lt%rows_cols = 0
    my_lt%n_cols    = n_cols
  end function LT_create

  subroutine LT_set_col(my_lt, col_ix, xx, yy)
    type(LT_table_t), intent(inout) :: my_lt
    integer, intent(in)            :: col_ix
    real(dp), intent(in)           :: xx(:), yy(:)
    my_lt%cols_rows(col_ix, :) = &
         get_spaced_data(xx, yy, LT_get_xdata(my_lt))
    my_lt%rows_cols(:, col_ix) = my_lt%cols_rows(col_ix, :)
  end subroutine LT_set_col

  subroutine LT_add_col(my_lt, xx, yy)
    type(LT_table_t), intent(inout) :: my_lt
    real(dp), intent(in)            :: xx(:), yy(:)
    type(LT_table_t)                :: temp_lt

    temp_lt = my_lt
    deallocate(my_lt%cols_rows)
    deallocate(my_lt%rows_cols)
    allocate(my_lt%cols_rows(my_lt%n_cols+1, my_lt%n_rows))
    allocate(my_lt%rows_cols(my_lt%n_rows, my_lt%n_cols+1))
    my_lt%cols_rows(1:my_lt%n_cols, :) = temp_lt%cols_rows
    my_lt%rows_cols(:, 1:my_lt%n_cols) = temp_lt%rows_cols
    my_lt%n_cols                       = my_lt%n_cols + 1
    my_lt%cols_rows(my_lt%n_cols, :)   = get_spaced_data(xx, yy, LT_get_xdata(my_lt))
    my_lt%rows_cols(:, my_lt%n_cols)   = my_lt%cols_rows(my_lt%n_cols, :)
  end subroutine LT_add_col

  subroutine LT_add_to_col(my_lt, col_ix, xx, yy)
    type(LT_table_t), intent(inout) :: my_lt
    integer, intent(in)            :: col_ix
    real(dp), intent(in)           :: xx(:), yy(:)
    my_lt%cols_rows(col_ix, :) = &
         my_lt%cols_rows(col_ix, :) + get_spaced_data(xx, yy, LT_get_xdata(my_lt))
    my_lt%rows_cols(:, col_ix) = my_lt%cols_rows(col_ix, :)
  end subroutine LT_add_to_col

  elemental function LT_get_loc(my_lt, x) result(my_loc)
    type(LT_table_t), intent(in) :: my_lt
    real(dp), intent(in)        :: x
    type(LT_loc_t)              :: my_loc
    my_loc = get_loc(my_lt, x)
  end function LT_get_loc

  elemental function get_loc(my_lt, x) result(my_loc)
    type(LT_table_t), intent(in) :: my_lt
    real(dp), intent(in)        :: x
    type(LT_loc_t)              :: my_loc
    real(dp)                    :: frac

    frac            = (x - my_lt%x_min) * my_lt%inv_dx
    my_loc%low_ix   = ceiling(frac)
    my_loc%low_frac = my_loc%low_ix - frac

    ! Check bounds
    if (my_loc%low_ix < 1) then
       my_loc%low_ix = 1
       my_loc%low_frac = 1
    else if (my_loc%low_ix >= my_lt%n_rows) then
       my_loc%low_ix = my_lt%n_rows - 1
       my_loc%low_frac = 0
    end if
  end function get_loc

  function LT_get_mcol(my_lt, x) result(col_values)
    type(LT_table_t), intent(in) :: my_lt
    real(dp), intent(in)        :: x
    real(dp)                    :: col_values(my_lt%n_cols)
    type(LT_loc_t)              :: loc
    loc = get_loc(my_lt, x)
    col_values = LT_get_mcol_at_loc(my_lt, loc)
  end function LT_get_mcol

  elemental function LT_get_col(my_lt, col_ix, x) result(col_value)
    type(LT_table_t), intent(in) :: my_lt
    integer, intent(in)         :: col_ix
    real(dp), intent(in)        :: x
    real(dp)                    :: col_value
    type(LT_loc_t)              :: loc
    loc       = get_loc(my_lt, x)
    col_value = LT_get_col_at_loc(my_lt, col_ix, loc)
  end function LT_get_col

  function LT_get_mcol_at_loc(my_lt, loc) result(col_values)
    type(LT_table_t), intent(in) :: my_lt
    type(LT_loc_t), intent(in)  :: loc
    real(dp)                    :: col_values(my_lt%n_cols)
    col_values = loc%low_frac * my_lt%cols_rows(:, loc%low_ix) + &
         (1-loc%low_frac) * my_lt%cols_rows(:, loc%low_ix+1)
  end function LT_get_mcol_at_loc

  elemental function LT_get_col_at_loc(my_lt, col_ix, loc) result(col_value)
    type(LT_table_t), intent(in) :: my_lt
    integer, intent(in)         :: col_ix
    type(LT_loc_t), intent(in)  :: loc
    real(dp)                    :: col_value
    col_value = loc%low_frac * my_lt%rows_cols(loc%low_ix, col_ix) + &
         (1-loc%low_frac) * my_lt%rows_cols(loc%low_ix+1, col_ix)
  end function LT_get_col_at_loc

  integer function LT_get_num_rows(my_lt)
    type(LT_table_t), intent(in) :: my_lt
    LT_get_num_rows = my_lt%n_rows
  end function LT_get_num_rows

  integer function LT_get_num_cols(my_lt)
    type(LT_table_t), intent(in) :: my_lt
    LT_get_num_cols = size(my_lt%cols_rows, 1)
  end function LT_get_num_cols

  subroutine LT_get_data(my_lt, x_data, cols_rows)
    type(LT_table_t), intent(in) :: my_lt
    real(dp), intent(out)      :: x_data(:), cols_rows(:, :)
    x_data  = LT_get_xdata(my_lt)
    cols_rows = my_lt%cols_rows
  end subroutine LT_get_data

  function LT_get_xdata(my_lt) result(xdata)
    type(LT_table_t), intent(in) :: my_lt
    real(dp)                    :: xdata(my_lt%n_rows)
    integer                     :: ix
    real(dp)                    :: dx
    dx = 1 / my_lt%inv_dx
    do ix = 1, my_lt%n_rows
       xdata(ix) = my_lt%x_min + (ix-1) * dx
    end do
  end function LT_get_xdata

  function get_spaced_data(in_xx, in_yy, new_xx) result(out_yy)
    real(dp), intent(in) :: in_xx(:), in_yy(:), new_xx(:)
    real(dp)             :: out_yy(size(new_xx))
    integer              :: ix
    do ix = 1, size(new_xx)
       call LT_lin_interp_list(in_xx, in_yy, new_xx(ix), out_yy(ix))
    end do
  end function get_spaced_data

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
       ix = FI_adaptive_r(x_list, x_value)
       temp = (x_value - x_list(ix-1)) / (x_list(ix) - x_list(ix-1))
       y_value = (1 - temp) * y_list(ix-1) + temp * y_list(ix)
    end if
  end subroutine LT_lin_interp_list

  subroutine LT_to_file(my_lt, filename)
    type(LT_table_t), intent(in) :: my_lt
    character(len=*), intent(in) :: filename
    integer                      :: my_unit

    open(newunit=my_unit, file=trim(filename), form='UNFORMATTED', &
         access='STREAM', status='REPLACE')
    write(my_unit) my_lt%n_rows, my_lt%n_cols
    write(my_unit) my_lt%x_min, my_lt%inv_dx
    write(my_unit) my_lt%cols_rows
    close(my_unit)
  end subroutine LT_to_file

  subroutine LT_from_file(my_lt, filename)
    type(LT_table_t), intent(inout) :: my_lt
    character(len=*), intent(in)    :: filename
    integer                         :: my_unit

    open(newunit=my_unit, file=trim(filename), form='UNFORMATTED', &
         access='STREAM', status='OLD')
    read(my_unit) my_lt%n_rows, my_lt%n_cols
    read(my_unit) my_lt%x_min, my_lt%inv_dx

    allocate(my_lt%cols_rows(my_lt%n_cols, my_lt%n_rows))
    allocate(my_lt%rows_cols(my_lt%n_rows, my_lt%n_cols))

    read(my_unit) my_lt%cols_rows
    my_lt%rows_cols = transpose(my_lt%cols_rows)

    close(my_unit)
  end subroutine LT_from_file

end module m_lookup_table
