! Test whether the different types can be included in the same program
program test_types_2d_3d
  use m_a2_types
  use m_a3_types

  type(a2_t) :: tree_2d
  type(a3_t) :: tree_3d

  ! This should give a warning (not initialized)
  call a2_print_info(tree_2d)
  call a3_print_info(tree_3d)

  print *, "Compilation: OK"
end program test_types_2d_3d
