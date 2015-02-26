module m_find_index

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! Public methods
  public :: FI_linear_r
  public :: FI_bsearch_r
  public :: FI_adaptive_r

contains

  ! Searches sorted list for the smallest ix such that
  ! val <= list(ix)
  ! On failure, returns size(list)+1
  function FI_linear_r(list, val) result(ix)
    real(dp), intent(in) :: list(:), val
    integer :: ix

    do ix = 1, size(list)
       if (val <= list(ix)) exit
    end do
  end function FI_linear_r

  ! Searches sorted list for the smallest ix such that
  ! val <= list(ix)
  ! On failure, returns size(list)+1
  function FI_bsearch_r(list, val) result(ix)
    real(dp), intent(in) :: list(:)
    real(dp), intent(in) :: val
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
  end function FI_bsearch_r

  ! Searches sorted list for the smallest ix such that
  ! val <= list(ix)
  ! On failure, returns size(list)+1
  ! Switches between linear and binary search
  function FI_adaptive_r(list, val) result(ix)
    real(dp), intent(in) :: list(:)
    real(dp), intent(in) :: val
    integer              :: ix
    integer, parameter   :: bsearch_limit = 40
    
    if (size(list) < bsearch_limit) then
       ix = FI_linear_r(list, val)
    else
       ix = FI_bsearch_r(list, val)
    end if
  end function FI_adaptive_r
  
end module m_find_index