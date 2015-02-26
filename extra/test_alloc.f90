program test_alloc
  implicit none

  integer :: i, j
  integer, parameter :: n = 1000*1000, m = 10
  integer, parameter :: dp = kind(0.0d0)

  ! type test_t
  !    real(dp), allocatable :: x1(:, :, :, :)
  !    real(dp), allocatable :: x2(:, :, :, :)
  !    real(dp), allocatable :: x3(:, :, :, :)
  !    real(dp), allocatable :: x4(:, :, :, :)
  ! end type test_t

  type test_t
     real(dp), dimension(m) :: x1, x2, x3, x4
  end type test_t

  type(test_t), allocatable :: my_array(:)

  allocate(my_array(n))

  ! do i = 1, n
  !    allocate(my_array(i)%x1(m, 1, 1, 1))
  !    allocate(my_array(i)%x2(m, 1, 1, 1))
  !    allocate(my_array(i)%x3(m, 1, 1, 1))
  !    allocate(my_array(i)%x4(m, 1, 1, 1))
  ! end do

  do i = 1, n
     my_array(i)%x1 = 0
     my_array(i)%x2 = 1
     my_array(i)%x3 = 2
     my_array(i)%x4 = 3
  end do
  do i = 1, n
     my_array(i)%x1 = sin(my_array(i)%x2)
     my_array(i)%x2 = cos(my_array(i)%x1)
     my_array(i)%x3 = tan(my_array(i)%x2)
     my_array(i)%x4 = atan(my_array(i)%x3)
  end do

  print *, size(my_array), my_array(1)%x1(1)
end program test_alloc