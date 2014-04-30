module m_lookup_table
   implicit none
   private

   integer, parameter :: dp = kind(1.0d0)

   type LT_xdata_t
      private
      integer               :: n_rows
      real(dp)              :: x_min, x_max, dx, inv_dx
      real(dp), allocatable :: x_values(:)
   end type LT_xdata_t

   type LT_col_t
      private
      type(LT_xdata_t)      :: xd
      real(dp), allocatable :: y_values(:)
   end type LT_col_t

   type LT_mcol_t
      private
      type(LT_xdata_t)      :: xd
      integer               :: n_cols
      real(dp), allocatable :: y_table(:, :)
   end type LT_mcol_t

   type LT_loc_t
      private
      integer  :: low_ix
      real(dp) :: low_frac
   end type LT_loc_t

   public :: LT_col_t
   public :: LT_mcol_t
   public :: LT_loc_t

   public :: LT_create_col
   public :: LT_create_mcol
   public :: LT_set_col
   public :: LT_set_mcol
   public :: LT_add_to_col
   public :: LT_add_to_mcol
   public :: LT_add_col
   public :: LT_get_loc_col
   public :: LT_get_loc_mcol
   public :: LT_get_col
   public :: LT_get_mcol
   public :: LT_get_1col
   public :: LT_get_col_at_loc
   public :: LT_get_mcol_at_loc
   public :: LT_get_1col_at_loc
   public :: LT_get_num_rows_col
   public :: LT_get_num_rows_mcol
   public :: LT_get_num_cols
   public :: LT_get_data_col
   public :: LT_get_data_mcol
   public :: LT_lin_interp_list

contains

   function LT_create_col(x_min, x_max, n_rows) result(my_lt)
      real(dp), intent(in) :: x_min, x_max
      integer, intent(in)  :: n_rows
      type(LT_col_t)       :: my_lt

      call set_xdata(my_lt%xd, x_min, x_max, n_rows)
      allocate(my_lt%y_values(n_rows))
      my_lt%y_values = 0
   end function LT_create_col

   function LT_create_mcol(x_min, x_max, n_rows, n_cols) result(my_lt)
      real(dp), intent(in) :: x_min, x_max
      integer, intent(in)  :: n_rows, n_cols
      type(LT_mcol_t)      :: my_lt

      call set_xdata(my_lt%xd, x_min, x_max, n_rows)
      allocate(my_lt%y_table(n_cols, n_rows))
      my_lt%y_table = 0
      my_lt%n_cols  = n_cols
   end function LT_create_mcol

   subroutine set_xdata(my_xd, x_min, x_max, n_rows)
      type(LT_xdata_t), intent(out) :: my_xd
      real(dp), intent(in)          :: x_min, x_max
      integer, intent(in)           :: n_rows
      integer                       :: ix
      if (x_max <= x_min) print *, "set_xdata: x_max should be > x_min"
      if (n_rows <= 1)    print *, "set_xdata: n_rows should be bigger than 1"

      my_xd%n_rows = n_rows
      my_xd%x_min  = x_min
      my_xd%x_max  = x_max
      my_xd%dx     = (x_max - x_min) / (n_rows - 1)
      my_xd%inv_dx = (n_rows - 1) / (x_max - x_min)
      allocate(my_xd%x_values(n_rows))
      my_xd%x_values = x_min + (/ (ix, ix = 0, n_rows-1) /) * my_xd%dx
   end subroutine set_xdata

   subroutine LT_set_col(my_lt, xx, yy)
      type(LT_col_t), intent(inout) :: my_lt
      real(dp), intent(in)          :: xx(:), yy(:)
      my_lt%y_values = get_spaced_data(xx, yy, my_lt%xd%x_values)
   end subroutine LT_set_col

   subroutine LT_set_mcol(my_lt, col_ix, xx, yy)
      type(LT_mcol_t), intent(inout) :: my_lt
      integer, intent(in)            :: col_ix
      real(dp), intent(in)           :: xx(:), yy(:)
      my_lt%y_table(col_ix, :) = get_spaced_data(xx, yy, my_lt%xd%x_values)
   end subroutine LT_set_mcol

   subroutine LT_add_col(my_lt, xx, yy)
      type(LT_mcol_t), intent(inout) :: my_lt
      real(dp), intent(in)           :: xx(:), yy(:)
      type(LT_mcol_t)                :: temp_lt

      temp_lt = my_lt
      deallocate(my_lt%y_table)
      allocate(my_lt%y_table(my_lt%n_cols+1, my_lt%xd%n_rows))
      my_lt%y_table(1:my_lt%n_cols, :) = temp_lt%y_table
      my_lt%n_cols                     = my_lt%n_cols + 1
      my_lt%y_table(my_lt%n_cols, :)   = get_spaced_data(xx, yy, my_lt%xd%x_values)
   end subroutine LT_add_col

   subroutine LT_add_to_col(my_lt, xx, yy)
      type(LT_col_t), intent(inout) :: my_lt
      real(dp), intent(in)          :: xx(:), yy(:)
      my_lt%y_values = my_lt%y_values + get_spaced_data(xx, yy, my_lt%xd%x_values)
   end subroutine LT_add_to_col

   subroutine LT_add_to_mcol(my_lt, col_ix, xx, yy)
      type(LT_mcol_t), intent(inout) :: my_lt
      integer, intent(in)            :: col_ix
      real(dp), intent(in)           :: xx(:), yy(:)
      my_lt%y_table(col_ix, :) = &
           my_lt%y_table(col_ix, :) + get_spaced_data(xx, yy, my_lt%xd%x_values)
   end subroutine LT_add_to_mcol

   elemental function LT_get_loc_col(my_lt, x) result(my_loc)
      type(LT_col_t), intent(in) :: my_lt
      real(dp), intent(in)       :: x
      type(LT_loc_t)             :: my_loc
      my_loc = get_loc(my_lt%xd, x)
   end function LT_get_loc_col

   elemental function LT_get_loc_mcol(my_lt, x) result(my_loc)
      type(LT_mcol_t), intent(in) :: my_lt
      real(dp), intent(in)        :: x
      type(LT_loc_t)              :: my_loc
      my_loc = get_loc(my_lt%xd, x)
   end function LT_get_loc_mcol

   elemental function get_loc(my_xd, x) result(my_loc)
      type(LT_xdata_t), intent(in) :: my_xd
      real(dp), intent(in)         :: x
      type(LT_loc_t)               :: my_loc
      real(dp)                     :: frac

      frac            = (x - my_xd%x_min) * my_xd%inv_dx
      my_loc%low_ix   = ceiling(frac)
      my_loc%low_frac = my_loc%low_ix - frac

      ! Check bounds
      if (my_loc%low_ix < 1) then
         my_loc%low_ix = 1
         my_loc%low_frac = 1
      else if (my_loc%low_ix >= my_xd%n_rows) then
         my_loc%low_ix = my_xd%n_rows - 1
         my_loc%low_frac = 0
      end if
   end function get_loc

   elemental function LT_get_col(my_lt, x) result(y_value)
      type(LT_col_t), intent(in) :: my_lt
      real(dp), intent(in)       :: x
      real(dp)                   :: y_value
      type(LT_loc_t)             :: loc
      loc = get_loc(my_lt%xd, x)
      y_value = LT_get_col_at_loc(my_lt, loc)
   end function LT_get_col

   elemental function LT_get_col_at_loc(my_lt, loc) result(y_value)
      type(LT_col_t), intent(in) :: my_lt
      type(LT_loc_t), intent(in) :: loc
      real(dp)                   :: y_value
      y_value = loc%low_frac * my_lt%y_values(loc%low_ix) + &
           (1-loc%low_frac) * my_lt%y_values(loc%low_ix+1)
   end function LT_get_col_at_loc

   function LT_get_mcol(my_lt, x) result(col_values)
      type(LT_mcol_t), intent(in) :: my_lt
      real(dp), intent(in)        :: x
      real(dp)                    :: col_values(my_lt%n_cols)
      type(LT_loc_t)              :: loc
      loc = get_loc(my_lt%xd, x)
      col_values = LT_get_mcol_at_loc(my_lt, loc)
   end function LT_get_mcol

   elemental function LT_get_1col(my_lt, col_ix, x) result(col_value)
      type(LT_mcol_t), intent(in) :: my_lt
      integer, intent(in)         :: col_ix
      real(dp), intent(in)        :: x
      real(dp)                    :: col_value
      type(LT_loc_t)              :: loc
      loc = get_loc(my_lt%xd, x)
      col_value = LT_get_1col_at_loc(my_lt, col_ix, loc)
   end function LT_get_1col

   function LT_get_mcol_at_loc(my_lt, loc) result(col_values)
      type(LT_mcol_t), intent(in) :: my_lt
      type(LT_loc_t), intent(in)  :: loc
      real(dp)                    :: col_values(my_lt%n_cols)
      col_values = loc%low_frac * my_lt%y_table(:, loc%low_ix) + &
           (1-loc%low_frac) * my_lt%y_table(:, loc%low_ix+1)
   end function LT_get_mcol_at_loc

   elemental function LT_get_1col_at_loc(my_lt, col_ix, loc) result(col_value)
      type(LT_mcol_t), intent(in) :: my_lt
      integer, intent(in)         :: col_ix
      type(LT_loc_t), intent(in)  :: loc
      real(dp)                    :: col_value
      col_value = loc%low_frac * my_lt%y_table(col_ix, loc%low_ix) + &
           (1-loc%low_frac) * my_lt%y_table(col_ix, loc%low_ix+1)
   end function LT_get_1col_at_loc

   integer function LT_get_num_rows_col(my_lt)
      type(LT_col_t), intent(in) :: my_lt
      LT_get_num_rows_col = my_lt%xd%n_rows
   end function LT_get_num_rows_col

   integer function LT_get_num_rows_mcol(my_lt)
      type(LT_mcol_t), intent(in) :: my_lt
      LT_get_num_rows_mcol = my_lt%xd%n_rows
   end function LT_get_num_rows_mcol

   integer function LT_get_num_cols(my_lt)
      type(LT_mcol_t), intent(in) :: my_lt
      LT_get_num_cols = size(my_lt%y_table, 1)
   end function LT_get_num_cols

   subroutine LT_get_data_col(my_lt, x_data, y_data)
      type(LT_col_t), intent(in) :: my_lt
      real(dp), intent(out)      :: x_data(:), y_data(:)
      x_data = my_lt%xd%x_values
      y_data = my_lt%y_values
   end subroutine LT_get_data_col

   subroutine LT_get_data_mcol(my_lt, x_data, y_table)
      type(LT_mcol_t), intent(in) :: my_lt
      real(dp), intent(out)      :: x_data(:), y_table(:, :)
      x_data  = my_lt%xd%x_values
      y_table = my_lt%y_table
   end subroutine LT_get_data_mcol

   function get_spaced_data(in_xx, in_yy, new_xx) result(out_yy)
      real(dp), intent(in) :: in_xx(:), in_yy(:), new_xx(:)
      real(dp)             :: out_yy(size(new_xx))
      integer              :: ix
      do ix = 1, size(new_xx)
         call LT_lin_interp_list(in_xx, in_yy, new_xx(ix), out_yy(ix))
      end do
   end function get_spaced_data

   ! Searches list for the interval containing value, such that list(i) <= value < list(i+1), and returns i
   integer function get_ix_in_list(list, value)
      real(dp), intent(in) :: list(:)
      real(dp), intent(in) :: value
      integer              :: iMin, iMax, iMiddle

      iMin = 1
      iMax = size(list)

      do
         iMiddle = iMin + (iMax - iMin) / 2
         if (value < list(iMiddle)) then
            iMax = iMiddle
         else
            iMin = iMiddle
         end if
         if (iMax - iMin <= 1) exit
      end do

      get_ix_in_list = iMin
   end function get_ix_in_list

   subroutine LT_lin_interp_list(x_list, y_list, x_value, y_value)
      real(dp), intent(in)    :: x_list(:), y_list(:)
      real(dp), intent(in)    :: x_value
      real(dp), intent(inout) :: y_value

      integer                 :: ix, iMin, iMax
      real(dp)                :: temp

      iMin = lbound(x_list, 1)
      iMax = ubound(x_list, 1)

      if (x_value <= x_list(iMin)) then
         y_value = y_list(iMin)
      else if (x_value >= x_list(iMax)) then
         y_value = y_list(iMax)
      else
         ix = get_ix_in_list(x_list, x_value)
         temp = (x_value - x_list(ix)) / (x_list(ix+1) - x_list(ix))
         y_value = (1 - temp) * y_list(ix) + temp * y_list(ix+1)
      end if
   end subroutine LT_lin_interp_list

end module m_lookup_table
