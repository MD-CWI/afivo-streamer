program simple_test
  use m_find_index

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: test_size = 1000

  integer  :: i
  integer  :: ix_linear, ix_bsearch, ix_adaptive
  real(dp) :: sorted_list(test_size)
  real(dp) :: search_value

  do i = 1, test_size
     sorted_list(i) = sqrt(real(i, dp))
  end do

  search_value = 0.5_dp * (sorted_list(1) + sorted_list(test_size))

  ix_linear   = find_index_linear(sorted_list, search_value)
  ix_bsearch  = find_index_bsearch(sorted_list, search_value)
  ix_adaptive = find_index_adaptive(sorted_list, search_value)

  if (ix_linear /= ix_bsearch .or. ix_linear /= ix_adaptive) then
     error stop "Indices are unequal"
  end if

  if (sorted_list(ix_linear) < search_value) then
     error stop "Wrong index found"
  end if

  if (ix_linear > 1) then
     if (sorted_list(ix_linear-1) >= search_value) then
        error stop "Wrong index found"
     end if
  endif

  print *, "Success! Found index", ix_linear

end program simple_test
