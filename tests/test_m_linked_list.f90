program test_m_linked_list
   use m_linked_list

   type(LL_int_head_t) :: my_list
   integer :: ix, my_value
   integer, parameter :: NN = 10
   logical :: success

   print *, "Lists do not need to be initialized"
   print *, "Initial size:", LL_get_size(my_list)

   print *, "Adding ", NN, " values to list"
   do ix = 1, NN
      call LL_add(my_list, ix)
   end do

   print *, "Popping values until the list is empty"

   do
      call LL_pop(my_list, my_value, success)
      if (.not. success) exit
      print *, my_value
   end do

   call LL_clear(my_list)

end program test_m_linked_list