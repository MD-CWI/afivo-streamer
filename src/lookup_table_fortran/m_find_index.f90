!> Module m_find_index contains functions for finding an index in a sorted list
module m_find_index

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! Public methods
  public :: find_index_linear
  public :: find_index_bsearch
  public :: find_index_adaptive

contains

  !> Linear search of sorted list for the smallest ix such that list(ix) >= val.
  !> On failure, returns size(list)+1
  function find_index_linear(list, val) result(ix)
    real(dp), intent(in) :: list(:) !< Sorted list
    real(dp), intent(in) :: val     !< Value to search for
    integer              :: ix

    do ix = 1, size(list)
       if (list(ix) >= val) exit
    end do
  end function find_index_linear

  !> Binary search of sorted list for the smallest ix such that list(ix) >= val.
  !> On failure, returns size(list)+1
  function find_index_bsearch(list, val) result(ix)
    real(dp), intent(in) :: list(:) !< Sorted list
    real(dp), intent(in) :: val     !< Value to search for
    integer              :: ix, i_min, i_max, i_middle

    i_min = 1
    i_max = size(list)

    do while (i_min < i_max)
       i_middle = i_min + (i_max - i_min) / 2

       if (val <= list(i_middle)) then
          i_max = i_middle
       else
          i_min = i_middle + 1
       end if
    end do

    ix = i_min
    if (val > list(ix)) ix = size(list) + 1
  end function find_index_bsearch

  !> Adaptive search (combination of linear and binary search) of sorted list
  !> for the smallest ix such that list(ix) >= val. On failure, returns
  !> size(list)+1
  function find_index_adaptive(list, val) result(ix)
    real(dp), intent(in) :: list(:) !< Sorted list
    real(dp), intent(in) :: val     !< Value to search for
    integer              :: ix
    integer, parameter   :: binary_search_limit = 40

    if (size(list) < binary_search_limit) then
       ix = find_index_linear(list, val)
    else
       ix = find_index_bsearch(list, val)
    end if
  end function find_index_adaptive

end module m_find_index
