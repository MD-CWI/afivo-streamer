program test
  use m_lookup_table
  implicit none

  integer, parameter :: dp          = kind(0.0d0)
  integer, parameter :: in_size(2)  = [53, 39]
  integer, parameter :: tbl_size(2) = [64, 8]

  integer     :: i, j
  type(LT2_t) :: lkp_tbl
  real(dp)    :: x1(in_size(1))
  real(dp)    :: x2(in_size(2))
  real(dp)    :: y_values(in_size(1), in_size(2))
  real(dp)    :: diff, mean_diff, xy(2)

  print *, 'Start test_lookup_table_2d'
  print *, 'Here, the input data has a different shape than the table,'
  print *, 'so linear interpolation is used to create the 2D table.'

  do i = 1, in_size(1)
     x1(i) = (i-1) * acos(-1.0_dp) / (in_size(1) - 1)
  end do

  do i = 1, in_size(2)
     x2(i) = (i-1) * acos(-1.0_dp) / (in_size(2) - 1)
  end do


  ! Create some testing data
  do j = 1, in_size(2)
     do i = 1, in_size(1)
        y_values(i, j) = sol(x1(i), x2(j))
     end do
  end do

  lkp_tbl = LT2_create([x1(1), x2(1)], &
       [x1(in_size(1)), x2(in_size(2))], tbl_size, 1)
  call LT2_set_col(lkp_tbl, 1, x1, x2, y_values)

  mean_diff = 0.0_dp

  do j = 1, in_size(2)
     do i = 1, in_size(1)
        diff = abs(y_values(i, j) - &
             LT2_get_col(lkp_tbl, 1, x1(i), x2(j)))
        mean_diff = mean_diff + diff / product(in_size)
     end do
  end do

  print *, "Mean difference at input points: ", mean_diff

  mean_diff = 0.0_dp

  do j = 1, in_size(2)
     do i = 1, in_size(1)
        call random_number(xy)
        xy = xy * acos(-1.0_dp)
        diff = abs(sol(xy(1), xy(2)) - &
             LT2_get_col(lkp_tbl, 1, xy(1), xy(2)))
        mean_diff = mean_diff + diff / product(in_size)
     end do
  end do

  print *, "Mean difference at random points: ", mean_diff

contains

  real(dp) function sol(x, y)
    real(dp), intent(in) :: x, y
    sol = 2 * x**2 + 3 * y
  end function sol

end program test
