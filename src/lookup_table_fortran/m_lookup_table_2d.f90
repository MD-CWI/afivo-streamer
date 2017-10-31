!> A Fortran 90 module for creating 2D lookup tables. These tables can be used
!> to efficiently interpolate one or more values.
!>
!> Author: Jannis Teunissen
module m_lookup_table_2d
  use m_lookup_table

  implicit none
  private

  ! The precision of the real numbers used in the tables
  integer, parameter :: dp = kind(1.0d0)

  !> The lookup table type. There can be one or more columns, for which values
  !> can be looked up for a given (x1, x2) coordinate
  type LT2_t
     integer  :: table_size(2) !< The size of the table
     integer  :: n_cols !< The number of columns/variables
     real(dp) :: x_min(2)  !< The minimum lookup coordinate
     real(dp) :: dx(2)     !< The x-spacing in the lookup coordinate
     real(dp) :: inv_dx(2) !< The inverse x-spacing
     real(dp), allocatable :: x1_data(:) !< List of x1 coordinates
     real(dp), allocatable :: x2_data(:) !< List of x2 coordinates

     real(dp), allocatable :: rows_cols(:, :, :)
  end type LT2_t

  !> Type to indicate a location in the lookup table, which can be used to speed
  !> up multiple lookups of different columns.
  type LT2_loc_t
     private
     !> The x-value lies between low_ix and low_ix+1
     integer  :: low_ix(2)
     !> The distance from low_ix (up to low_ix+1), given as a real number
     !> between 0 and 1.
     real(dp) :: low_frac(2)
  end type LT2_loc_t

  ! Public types
  public :: LT2_t
  public :: LT2_loc_t

  ! Public methods
  public :: LT2_create           ! Create a new lookup table
  public :: LT2_set_col          ! Set one table column
  public :: LT2_get_loc          ! Get the index (row) of a value
  public :: LT2_get_col          ! Interpolate one column
  public :: LT2_get_col_at_loc   ! Get one column at location

contains

  !> This function returns a new lookup table
  function LT2_create(x_min, x_max, table_size, n_cols) result(my_lt)
    real(dp), intent(in) :: x_min(2) !< Minimum coordinate
    real(dp), intent(in) :: x_max(2) !< Maximum coordinate
    integer, intent(in)  :: table_size(2) !< How many values to store
    integer, intent(in)  :: n_cols !< Number of variables that will be looked up
    integer :: ix
    type(LT2_t) :: my_lt

    if (any(x_max <= x_min)) stop "LT2_create error: x_max <= x_min"
    if (any(table_size <= 1)) stop "LT2_create error: table_size <= 1"

    my_lt%table_size = table_size
    my_lt%x_min  = x_min
    my_lt%dx     = (x_max - x_min) / (table_size - 1)
    my_lt%inv_dx = 1 / my_lt%dx

    allocate(my_lt%x1_data(table_size(1)))
    allocate(my_lt%x2_data(table_size(2)))

    do ix = 1, table_size(1)
       my_lt%x1_data(ix) = x_min(1) + (ix-1) * my_lt%dx(1)
    end do

    do ix = 1, table_size(2)
       my_lt%x2_data(ix) = x_min(2) + (ix-1) * my_lt%dx(2)
    end do

    allocate(my_lt%rows_cols(table_size(1), table_size(2), n_cols))
    my_lt%rows_cols = 0
    my_lt%n_cols    = n_cols
  end function LT2_create

  !> Fill the column with index col_ix using linearly interpolated data
  pure subroutine LT2_set_col(my_lt, col_ix, x1, x2, y)
    type(LT2_t), intent(inout) :: my_lt
    integer, intent(in)                 :: col_ix
    real(dp), intent(in)                :: x1(:), x2(:), y(:, :)
    real(dp), allocatable :: tmp(:, :)
    integer :: ix

    allocate(tmp(my_lt%table_size(1), size(x2)))

    ! Interpolate along first coordinate
    do ix = 1, size(x2)
       tmp(:, ix) = LT_get_spaced_data(x1, y(:, ix), my_lt%x1_data)
    end do

    ! Interpolate along second coordinate
    do ix = 1, my_lt%table_size(1)
       my_lt%rows_cols(ix, :, col_ix) = &
            LT_get_spaced_data(x2, tmp(ix, :), my_lt%x2_data)
    end do
  end subroutine LT2_set_col

  !> Get a location in the lookup table
  pure function LT2_get_loc(my_lt, x) result(my_loc)
    type(LT2_t), intent(in) :: my_lt
    real(dp), intent(in)             :: x(2)
    type(LT2_loc_t)                   :: my_loc
    real(dp)                         :: frac(2)

    frac            = (x - my_lt%x_min) * my_lt%inv_dx
    my_loc%low_ix   = ceiling(frac)
    my_loc%low_frac = my_loc%low_ix - frac

    ! Check bounds
    where (my_loc%low_ix < 1)
       my_loc%low_ix   = 1
       my_loc%low_frac = 1
    end where

    where (my_loc%low_ix >= my_lt%table_size)
       my_loc%low_ix   = my_lt%table_size - 1
       my_loc%low_frac = 0
    end where
  end function LT2_get_loc

  !> Get the value of a single column at x
  pure function LT2_get_col(my_lt, col_ix, x) result(col_value)
    type(LT2_t), intent(in) :: my_lt
    integer, intent(in)         :: col_ix
    real(dp), intent(in)        :: x(2)
    real(dp)                    :: col_value
    type(LT2_loc_t)              :: loc

    loc       = LT2_get_loc(my_lt, x)
    col_value = LT2_get_col_at_loc(my_lt, col_ix, loc)
  end function LT2_get_col

  !> Get the value of a single column at a location
  pure function LT2_get_col_at_loc(my_lt, col_ix, loc) result(col_value)
    type(LT2_t), intent(in) :: my_lt
    integer, intent(in)              :: col_ix
    type(LT2_loc_t), intent(in)       :: loc
    integer :: ix(2)
    real(dp)                         :: w(2, 2)
    real(dp)                         :: col_value

    ! Bilinear interpolation
    w(1, 1) = loc%low_frac(1) * loc%low_frac(2)
    w(2, 1) = (1 - loc%low_frac(1)) * loc%low_frac(2)
    w(1, 2) = loc%low_frac(1) * (1 - loc%low_frac(2))
    w(2, 2) = (1 - loc%low_frac(1)) * (1 - loc%low_frac(2))
    ix = loc%low_ix

    col_value = sum(w * my_lt%rows_cols(ix(1):ix(1)+1, &
         ix(2):ix(2)+1, col_ix))
  end function LT2_get_col_at_loc

end module m_lookup_table_2d
