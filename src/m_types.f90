!> Module with basic types
module m_types

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  !> Undefined string
  character(len=*), parameter :: undefined_str = "UNDEFINED"

  !> Undefined number
  real(dp), parameter :: undefined_real = -1e100_dp

  !> Default length of strings
  integer, parameter :: string_len = 200

  !> Default length of names
  integer, parameter :: name_len = 20

  !> Default length of component names
  integer, parameter :: comp_len = 20

end module m_types
