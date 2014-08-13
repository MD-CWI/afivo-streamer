module m_afivo

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: a5_size = 8

  type, abstract :: bnd_t
   contains
     procedure, deferred :: fill
  end type bnd_t

  type box_t
     integer               :: i_min(3)
     type(box_t), pointer  :: parent
     type(box_t), pointer  :: children(4)
     type(bnd_t)           :: bounds(4)
     real(dp), allocatable :: vars(:, :, :)
  end type box_t

  type box_ptr
     type(box_t), pointer  :: ptr
  end type box_ptr

  type box_array_t
     type(box_ptr), allocatable :: boxes(:)
  end type box_array_t

contains

  subroutine a5_create_box(bx, i_min, bounds, n_vars)
    type(box_t), intent(out) :: bx
    integer, intent(in)      :: i_min(3), bounds(4)

    bx%i_min  = i_min
    bx%bounds = bounds
    nullify(bx%parent)
    nullify(bx%children)
    allocate(bx%vars(0:a5_size+1, 0:a5_size+1, n_vars))
    bx%vars   = 0
  end subroutine a5_create_box

  subroutine a5_destroy_box(bx)
    type(box_t), intent(inout) :: bx
    if (associated(bx%children)) stop "still have children"
    deallocate(vars)
  end subroutine a5_destroy_box

  subroutine a5_get_child_array(bxa, bxa_children)
    type(box_array_t), intent(in) :: bxa
    type(box_array_t), intent(out) :: bxa_children

    ...
  end subroutine a5_get_child_array

  subroutine a5_prolong()

  end subroutine a5_prolong

  subroutine a5_restrict()
    
  end subroutine a5_restrict

end module m_afivo