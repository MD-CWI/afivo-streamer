program test
  use m_lookup_table
  implicit none

  integer, parameter   :: dp             = kind(0.0d0)
  integer, parameter   :: test_size      = 1000*1000
  integer, parameter   :: min_table_size = 10

  integer              :: i, cntr, table_size
  type(lookup_table_t) :: lkp_tbl
  real(dp)             :: x_values(test_size)
  real(dp)             :: y_values(test_size)
  real(dp)             :: y2_values(2, test_size)
  real(dp)             :: lkp_results(test_size)
  real(dp)             :: lkp2_results(2, test_size)
  real(dp)             :: max_diff
  real(dp)             :: total_time, time_t1, time_t2

  print *, 'start test_lookup_table'

  ! Create some testing data
  do i = 1, test_size
     x_values(i) = i / (0.1_dp * test_size)
     y_values(i) = sin(x_values(i))
     y2_values(:, i) = (/sin(x_values(i)), atan(x_values(i))/)
  end do

  print *, "Testing m_lookup_table.f90 implementation..."
  print *, ""
  print *, " *** Vectorized calls (one call for all the data) ***"
  print *, "Iteration    table_size   max. difference"

  table_size = min_table_size
  cntr = 0
  total_time = 0.0_dp

  do
     cntr = cntr + 1
     table_size = table_size * 2
     if (table_size > test_size) exit

     ! Create lookup table, last two arguments are optional (defaults 1.0 and 1)
     ! First arguments are min x, max x, num points
     lkp_tbl = LT_create(x_values(1), x_values(test_size), table_size, 1)

     ! Add the data to the table
     call LT_set_col(lkp_tbl, 1, x_values, y_values)

     call cpu_time(time_t1)
     lkp_results = LT_get_col(lkp_tbl, 1, x_values)
     call cpu_time(time_t2)
     total_time = total_time + (time_t2 - time_t1)

     max_diff = maxval(abs(lkp_results-y_values))
     print *, cntr, table_size, max_diff
  end do

  print *, "You should see 2nd order convergence"
  print *, "Number of lookups performed:", cntr * test_size
  print *, "Number of lookups / second:", &
       (cntr * test_size) / (total_time + epsilon(1.0_dp))
  print *, "Nanosecond per lookup:", &
       1e9_dp * total_time / (cntr * test_size)
  print *, ""

  print *, ""
  print *, "Testing multiple column mode (2 columns)"
  print *, " *** Vectorized calls (one call for all the data) ***"
  print *, "Iteration    table_size   max. difference"

  table_size = min_table_size
  cntr = 0
  total_time = 0.0_dp

  do
     cntr = cntr + 1
     table_size = table_size * 2
     if (table_size > test_size) exit

     ! Create lookup table, last two arguments are optional (defaults 1.0 and 1)
     ! First arguments are min x, max x, num points
     lkp_tbl = LT_create(x_values(1), x_values(test_size), table_size, 2)

     ! Add the data to the table
     call LT_set_col(lkp_tbl, 1, x_values, y2_values(1, :))
     call LT_set_col(lkp_tbl, 2, x_values, y2_values(2, :))

     call cpu_time(time_t1)
     do i = 1, test_size
        lkp2_results(:, i) = LT_get_mcol(lkp_tbl, x_values(i))
     end do
     call cpu_time(time_t2)
     total_time = total_time + (time_t2 - time_t1)

     max_diff = maxval(abs(lkp2_results-y2_values))
     print *, cntr, table_size, max_diff
  end do

  print *, "You should see 2nd order convergence"
  print *, "Number of lookups performed:", cntr * test_size
  print *, "Number of lookups / second:", &
       (cntr * test_size) / (total_time + epsilon(1.0_dp))
  print *, "Nanosecond per lookup:", &
       1e9_dp * total_time / (cntr * test_size)
  print *, ""

end program test
