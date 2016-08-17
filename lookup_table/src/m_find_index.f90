!> Module that contains subroutines for finding special indices in
!  monotonically increasing sorted lists
module m_find_index

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! Public methods
  public :: FI_linear_r
  public :: FI_binary_search_r
  public :: FI_adaptive_r

contains

  !> Searches monotonically increasing sorted list for the smallest ix
  !  such that rvalue <= list(ix)
  !  On failure, returns size(list)+1
  function FI_linear_r(list, rvalue) result(ix)
    real(dp), intent(in) :: list(:)    !< Increasing sorted list
    real(dp), intent(in) :: rvalue     !< Given real value
    integer :: ix

    do ix = 1, size(list)
       if (rvalue <= list(ix)) exit
    end do
  end function FI_linear_r

  !> Searches monotonically increasing sorted list for the ix which is
  !  'closest' to rvalue. If rvalue is not in the list list(ix) will be the
  !  first value larger than rvalue.
  !  On failure, returns size(list)+1
  function FI_binary_search_r(list, rvalue) result(ix)
    real(dp), intent(in) :: list(:)    !< Increasing sorted list
    real(dp), intent(in) :: rvalue     !< Given real value
    integer              :: ix, i_min, i_max, i_middle

    i_min = 1
    i_max = size(list)

    do while (i_min < i_max)
       i_middle = i_min + (i_max - i_min) / 2
       
       if (rvalue <= list(i_middle)) then
          i_max = i_middle
       else
          i_min = i_middle + 1
       end if
    end do

    ix = i_min
    if (rvalue > list(ix)) ix = size(list) + 1
  end function FI_binary_search_r

  !> Searches monotonically increasing sorted list for the smallest ix
  !  such that rvalue <= list(ix)
  !  Switches between linear and binary search
  !  On failure, returns size(list)+1
  function FI_adaptive_r(list, rvalue) result(ix)
    real(dp), intent(in) :: list(:)    !< Increasing sorted list
    real(dp), intent(in) :: rvalue     !< Given real value
    integer              :: ix
    integer, parameter   :: binary_search_limit = 40
    
    if (size(list) < binary_search_limit) then
       ix = FI_linear_r(list, rvalue)
    else
       ix = FI_binary_search_r(list, rvalue)
    end if
  end function FI_adaptive_r
  
end module m_find_index
