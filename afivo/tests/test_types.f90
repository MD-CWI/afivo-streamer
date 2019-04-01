! Test whether the different types can be included in the same program
program test_types
  use m_af_types

  implicit none

  type(af_t) :: tree_2d

  ! This should give a warning (not initialized)
  call af_print_info(tree_2d)

  print *, "Compilation: OK"
end program test_types
