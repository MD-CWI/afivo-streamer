program usage_example
  use m_lookup_table

  implicit none
  integer, parameter :: dp         = kind(0.0d0)
  integer, parameter :: table_size = 99
  integer, parameter :: test_size  = 100*1000

  call test_linear_data(table_size, test_size)
  call test_linear_extrapolation(table_size, test_size)
  call test_valid_location(table_size, test_size)

contains

  subroutine test_linear_data(table_size, n_samples)
    integer, intent(in)   :: table_size, n_samples
    type(LT_t)            :: my_lt
    real(dp)              :: x(2)        = [-1.0_dp, 1.0_dp]
    real(dp)              :: x_sample(2) = [-2.0_dp, 2.0_dp]
    real(dp), allocatable :: x_test(:), y(:)
    real(dp)              :: max_deviation
    real(dp), parameter   :: tolerance = 5e-14_dp

    print *, "test_linear_data"

    my_lt = LT_create(x(1), x(2), n_points=table_size, n_cols=1, &
       xspacing=LT_xspacing_linear)

    call LT_set_col(my_lt, 1, x, f_linear(x))

    allocate(x_test(n_samples))
    allocate(y(n_samples))

    call random_number(x_test)
    x_test = x_sample(1) + (x_sample(2) - x_sample(1)) * x_test

    where (x_test < x(1))
       y = f_linear(x(1))
    elsewhere (x_test > x(2))
       y = f_linear(x(2))
    elsewhere
       y = f_linear(x_test)
    end where

    max_deviation = maxval(abs(LT_get_col(my_lt, 1, x_test) - y))

    if (max_deviation > tolerance) then
       print *, "FAILED: too large deviation from solution", max_deviation
    else
       print *, "PASSED"
    end if

  end subroutine test_linear_data

  subroutine test_linear_extrapolation(table_size, n_samples)
    integer, intent(in)   :: table_size, n_samples
    type(LT_t)            :: my_lt
    real(dp)              :: x(2)        = [-1.0_dp, 1.0_dp]
    real(dp)              :: x_sample(2) = [-2.0_dp, 2.0_dp]
    real(dp), allocatable :: x_test(:), y(:)
    real(dp)              :: max_deviation
    real(dp), parameter   :: tolerance = 5e-14_dp

    print *, "test_linear_extrapolation"

    my_lt = LT_create(x(1), x(2), n_points=table_size, n_cols=1, &
       xspacing=LT_xspacing_linear, extrapolate_above=.true.)

    call LT_set_col(my_lt, 1, x, f_linear(x))

    allocate(x_test(n_samples))
    allocate(y(n_samples))

    call random_number(x_test)
    x_test = x_sample(1) + (x_sample(2) - x_sample(1)) * x_test

    where (x_test < x(1))
       y = f_linear(x(1))
    elsewhere
       y = f_linear(x_test)
    end where

    max_deviation = maxval(abs(LT_get_col(my_lt, 1, x_test) - y))

    if (max_deviation > tolerance) then
       print *, "FAILED: too large deviation from solution", max_deviation
    else
       print *, "PASSED"
    end if

  end subroutine test_linear_extrapolation

  subroutine test_valid_location(table_size, n_samples)
    integer, intent(in)         :: table_size, n_samples
    type(LT_t)                  :: my_lt
    real(dp)                    :: x(2)        = [-1.0_dp, 1.0_dp]
    real(dp)                    :: x_sample(2) = [-2.0_dp, 2.0_dp]
    real(dp), allocatable       :: x_test(:)
    type(LT_loc_t), allocatable :: locs(:)
    integer                     :: n
    logical                     :: success

    print *, "test_valid_location"

    my_lt = LT_create(x(1), x(2), n_points=table_size, n_cols=0, &
         xspacing=LT_xspacing_linear)

    allocate(x_test(n_samples))
    call random_number(x_test)
    x_test = x_sample(1) + (x_sample(2) - x_sample(1)) * x_test

    allocate(locs(n_samples))
    locs = LT_get_loc(my_lt, x_test)

    success = .true.
    if (maxval(locs(:)%low_ix) > table_size - 1) then
       print *, "FAILED: too high low_ix in location"
       success = .false.
    end if

    if (minval(locs(:)%low_ix) < 1) then
       print *, "FAILED: too low low_ix in location"
       success = .false.
    end if

    do n = 1, n_samples
       if (x_test(n) < x(1)) then
          if (locs(n)%low_ix /= 1) then
             print *, "FAILED: low_ix /= 1 for x < x_min"
             exit
          end if
       else if (x_test(n) > x(2)) then
          if (locs(n)%low_ix /= table_size-1) then
             print *, "FAILED: low_ix /= table_size-1 for x > x_max"
             exit
          end if
       else if (my_lt%x(locs(n)%low_ix) > x_test(n) .or. &
            my_lt%x(locs(n)%low_ix+1) < x_test(n)) then
          print *, "FAILED: x not between x(low_ix) and x(low_ix+1)"
          exit
       end if
    end do

    if (n /= n_samples + 1) success = .false.
    if (success) print *, "PASSED"

  end subroutine test_valid_location

  real(dp) elemental function f_linear(x)
    real(dp), intent(in) :: x
    f_linear = -1.5_dp * x + 0.5_dp
  end function f_linear

end program
