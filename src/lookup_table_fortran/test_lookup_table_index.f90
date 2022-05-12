program test_lookup_table_index
  use m_lookup_table

  implicit none
  integer, parameter    :: dp        = kind(0.0d0)
  integer, parameter    :: n_points  = 1024
  integer, parameter    :: n_samples = 2**20
  real(dp), allocatable :: x(:)
  integer               :: n, k
  type(LT_t)            :: lt(3)
  type(LT_loc_t)        :: loc

  lt(1) = LT_create(0.0_dp, 1.0_dp, n_points, n_cols=0, &
       xspacing=LT_xspacing_linear)
  lt(2) = LT_create(0.0_dp, 1.0_dp, n_points, n_cols=0, &
       xspacing=LT_xspacing_quadratic)
  lt(3) = LT_create(0.0_dp, 1.0_dp, n_points, n_cols=0, &
       xspacing=LT_xspacing_cubic)

  allocate(x(n_samples))
  call random_number(x)

  do k = 1, 3
     do n = 1, n_samples
        loc = LT_get_loc(lt(k), x(n))

        ! Check if xa <= x <= xb
        if (x(n) < lt(k)%x(loc%low_ix) .or. x(n) > lt(k)%x(loc%low_ix+1)) then
           print *, "For table with xspacing of order", k
           print *, lt(k)%x(loc%low_ix), x(n), lt(k)%x(loc%low_ix+1)
           error stop "FAIL (should have xa <= x <= xb)"
        end if
     end do
  end do

  print *, "PASS"

end program test_lookup_table_index
