module m_linked_list

   implicit none
   private

   integer, parameter :: dp = kind(0.0d0)

   type LL_int_t
      private
      integer                 :: x
      type(LL_int_t), pointer :: next
   end type LL_int_t

   type LL_int_head_t
      private
      type(LL_int_t), pointer :: ptr => null()
      integer                 :: n_values = 0
   end type LL_int_head_t

   ! Public types
   public :: LL_int_head_t

   ! Public methods
   public :: LL_clear
   public :: LL_pop
   public :: LL_add
   public :: LL_get_size

contains

   subroutine LL_clear(head)
      type(LL_int_head_t), intent(inout) :: head
      type(LL_int_t), pointer            :: this_element
      integer                            :: ix

      do ix = 1, head%n_values
         this_element => head%ptr
         head%ptr     => this_element%next
         deallocate(this_element)
      end do

      head%n_values = 0
   end subroutine LL_clear

   subroutine LL_add(head, x)
      type(LL_int_head_t)     :: head
      integer, intent(in)     :: x
      type(LL_int_t), pointer :: new_element

      allocate(new_element)
      new_element%next => head%ptr
      new_element%x = x
      head%ptr => new_element
      head%n_values = head%n_values + 1
   end subroutine LL_add

   subroutine LL_pop(head, popped_value, success)
      type(LL_int_head_t), intent(inout) :: head
      integer, intent(inout)             :: popped_value
      logical, intent(out)               :: success
      type(LL_int_t), pointer            :: element_to_pop

      success = head%n_values > 0

      if (success) then
         element_to_pop => head%ptr
         head%ptr       => element_to_pop%next
         popped_value   = element_to_pop%x
         head%n_values  = head%n_values - 1
         deallocate(element_to_pop)
      end if
   end subroutine LL_pop

   function LL_get_size(head) result(list_size)
      type(LL_int_head_t), intent(in) :: head
      integer                         :: list_size

      list_size = head%n_values
   end function LL_get_size

end module m_linked_list