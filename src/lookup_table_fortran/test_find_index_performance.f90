program performance_test
  use m_find_index
  implicit none

  integer, parameter :: dp         = kind(0.0d0)
  integer, parameter :: min_size   = 10
  integer, parameter :: max_size   = 100
  integer, parameter :: delta_size = 10
  integer, parameter :: n_finds    = 1000*1000

  integer               :: i, ix_1, ix_2, ix_3, test_size, n
  real(dp), allocatable :: list(:)
  real(dp), allocatable :: search_vals(:)
  real(dp)              :: t_start, t_end, t_linear
  real(dp)              :: t_bsearch, t_adaptive

  allocate(search_vals(n_finds))
  print *, "Test size | correct? | t_linear | t_bsearch | t_adaptive"

  do test_size = min_size, max_size, delta_size
     allocate(list(test_size))
     call random_number(search_vals)
     search_vals = search_vals * test_size

     ! Create sorted list
     call random_number(list)
     list = nint(list)
     do i = 2, test_size
        list(i) = list(i-1) + list(i)
     end do

     call cpu_time(t_start)
     do n = 1, n_finds
        ix_1 = find_index_linear(list, search_vals(n))
     end do
     call cpu_time(t_end)
     t_linear = t_end - t_start

     call cpu_time(t_start)
     do n = 1, n_finds
        ix_2 = find_index_bsearch(list, search_vals(n))
     end do
     call cpu_time(t_end)
     t_bsearch = t_end - t_start

     call cpu_time(t_start)
     do n = 1, n_finds
        ix_3 = find_index_adaptive(list, search_vals(n))
     end do
     call cpu_time(t_end)
     t_adaptive = t_end - t_start

     write(*, '(i9,l10,3e12.2)') test_size, &
          ix_1 == ix_2 .and. ix_2 == ix_3, &
          t_linear, t_bsearch, t_adaptive
     deallocate(list)
  end do

end program performance_test
