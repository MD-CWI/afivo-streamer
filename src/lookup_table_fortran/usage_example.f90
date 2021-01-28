program usage_example
  use m_lookup_table

  implicit none
  integer, parameter   :: dp = kind(0.0d0)

  real(dp), parameter :: x_data(4) = [1.0_dp, 2.75_dp, 3.3_dp, 5.6_dp]
  type(LT_t)          :: my_lt

  print *, "This program demonstrates the usage of m_lookup_table"
  print *, ""

  print *, "Create a lookup table with the following options"
  print *, "   x_min    = ", x_data(1)
  print *, "   x_max    = ", x_data(4)
  print *, "   n_points = ", 20
  print *, "   n_cols   = ", 2
  my_lt = LT_create(x_data(1), x_data(4), n_points=20, n_cols=2)

  print *, "The table input is unevenly spaced, using these x-values:"
  write(*, "(4e12.4)") x_data
  print *, ""

  print *, "Setting the first column to y1 = x**2"
  call LT_set_col(my_lt, 1, x_data, x_data**2)

  print *, "Setting the second column to y2 = sin(x)"
  call LT_set_col(my_lt, 2, x_data, sin(x_data))

  print *, "Adding a third column (which causes re-allocations),"
  print *, "which is set to y3 = sqrt(x)"
  call LT_add_col(my_lt, x_data, sqrt(x_data))
  print *, ""

  call example_1(my_lt)
  call example_2(my_lt)
  call example_3(my_lt)
  call example_4(my_lt)
  call example_5(my_lt)

contains

  subroutine example_1(lt)
    type(LT_t), intent(in) :: lt

    print *, "----------------------------------------"
    print *, "Example 1: Get values from a single column"
    print *, ""
    print *, "3**2 is about:       ", LT_get_col(lt, 1, 3.0_dp)
    print *, "sin(3) is about:     ", LT_get_col(lt, 2, 3.0_dp)
    print *, "sqrt(3) is about:    ", LT_get_col(lt, 3, 3.0_dp)
    print *, "----------------------------------------"
    print *, ""
  end subroutine example_1

  subroutine example_2(lt)
    type(LT_t), intent(in) :: lt
    type(LT_loc_t)         :: loc

    print *, "----------------------------------------"
    print *, "Example 2: Get values using a loc (location) object, which"
    print *, "can be faster if you look up at the same x-value often"
    print *, ""

    loc = LT_get_loc(my_lt, 3.0_dp)
    print *, "3**2 is about:       ", LT_get_col_at_loc(lt, 1, loc)
    print *, "sin(3) is about:     ", LT_get_col_at_loc(lt, 2, loc)
    print *, "sqrt(3) is about:    ", LT_get_col_at_loc(lt, 3, loc)
    print *, "----------------------------------------"
    print *, ""
  end subroutine example_2

  subroutine example_3(lt)
    type(LT_t), intent(in) :: lt
    real(dp)               :: columns(lt%n_cols)

    print *, "----------------------------------------"
    print *, "Example 3: Get all the columns at the same time."
    print *, "This can also be done using a loc object"
    print *, ""

    columns = LT_get_mcol(my_lt, 3.0_dp)
    print *, "3**2 is about:       ", columns(1)
    print *, "sin(3) is about:     ", columns(2)
    print *, "sqrt(3) is about:    ", columns(3)
    print *, "----------------------------------------"
    print *, ""
  end subroutine example_3

  subroutine example_4(lt)
    type(LT_t), intent(in) :: lt

    print *, "----------------------------------------"
    print *, "Example 4: When x > x_max or x < x_min, a lookup"
    print *, "returns y(x_max) and y(x_min), respectively"
    print *, ""

    print *, "The table starts at 1.0, so a lookup at 0.0 gives"
    print *, "LT_get_col(lt, 1, 0.0_dp):", LT_get_col(lt, 1, 0.0_dp)
    print *, "which is the same as"
    print *, "LT_get_col(lt, 1, 1.0_dp):", LT_get_col(lt, 1, 1.0_dp)
    print *, ""
    print *, "The table ends at 5.6, so a lookup at 10.0 gives"
    print *, "LT_get_col(lt, 1, 10.0_dp):", LT_get_col(lt, 1, 10.0_dp)
    print *, "which is the same as"
    print *, "LT_get_col(lt, 1, 5.6_dp): ", LT_get_col(lt, 1, 5.6_dp)
    print *, "----------------------------------------"
    print *, ""
  end subroutine example_4

  subroutine example_5(lt)
    type(LT_t), intent(in)      :: lt
    type(LT_t)                  :: new_lt
    character(len=*), parameter :: filename = "saved_table.dat"

    print *, "----------------------------------------"
    print *, "Example 5: You can store the lookup table in a file,"
    print *, "and read it back in. The table is stored in binary,"
    print *, "which is not always portable between machines."

    print *, "Storing the lookup table in    ", filename
    call LT_to_file(lt, filename)

    print *, "Reading in a lookup table from ", filename
    call LT_from_file(new_lt, filename)

    print *, ""
    print *, "The results should be the same:"
    print *, "3**2 is about:       ", LT_get_col(new_lt, 1, 3.0_dp)
    print *, "sin(3) is about:     ", LT_get_col(new_lt, 2, 3.0_dp)
    print *, "sqrt(3) is about:    ", LT_get_col(new_lt, 3, 3.0_dp)
    print *, "----------------------------------------"
    print *, ""
  end subroutine example_5

end program usage_example
