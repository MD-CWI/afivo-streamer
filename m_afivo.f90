module m_afivo

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: a5_size = 8

  type, abstract :: bnd_t
   contains
     procedure, deferred :: fill
  end type bnd_t

  type box2_t
     integer               :: lvl
     integer               :: i_min(2)
     integer               :: parent
     integer               :: children(4)
     integer               :: neighbors(4)
     real(dp), allocatable :: vars(:, :, :)
  end type box2_t

  type box2_store
     ! For now a simple array
     integer :: n_boxes
     type(box2_t), allocatable :: bxs(:)
  end type box2_store

contains

  elemental function a5_get_box2(box_store, box_id) result(box)
    integer, intent(in)            :: box_id
    type(box2_store_t), intent(in) :: box_store
    type(box2_t)                   :: box
    box = box_store%bxs(box_id)
  end function a5_get_box

  subroutine a5_set_box2(box_store, box_id, box)
    integer, intent(in)               :: box_id
    type(box2_store_t), intent(inout) :: box_store
    type(box2_t), intent(in)          :: box
    box_store%bxs(box_id) = box
  end subroutine a5_set_box2

  subroutine a5_add_box2(box_store, box)
    type(box2_store_t), intent(inout) :: box_store
    type(box2_t), intent(in)          :: box
    integer :: box_id
    ! TODO: lock store
    box_id = get_free_box_id(box_store)
    box_store%bxs(box_id) = box
  end subroutine a5_add_box2

  subroutine a5_create_box(box, i_min, bounds, n_vars)
    type(box2_t), intent(out) :: box
    integer, intent(in)      :: i_min(3), bounds(4)

    box%i_min  = i_min
    box%bounds = bounds
    nullify(box%parent)
    nullify(box%children)
    allocate(box%vars(0:a5_size+1, 0:a5_size+1, n_vars))
    box%vars   = 0
  end subroutine a5_create_box

  subroutine a5_destroy_box(box)
    type(box2_t), intent(inout) :: box
    if (associated(box%children)) stop "still have children"
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