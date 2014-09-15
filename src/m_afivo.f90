module m_afivo

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  integer, parameter :: a5_no_neighbor = 0

  integer, parameter :: a2_ch_rev(4, 2) = (/ (/2, 1, 4, 3/), (/3, 4, 1, 2/) /)
  logical, parameter :: a2_ch_low(4, 2) = (/ &
       (/.true., .false., .true., .false./), &
       (/.true., .true., .false., .false./), /)
  logical, parameter :: a2_nb_low(4) = (/.true., .false., .true., .false./)
  integer, parameter :: a2_nb_rev(4) = (/2, 1, 4, 3/)
  integer, parameter :: a2_nb_dim(4) = (/1, 1, 2, 2/)

  integer, parameter :: nb_lx = 1
  integer, parameter :: nb_hx = 2
  integer, parameter :: nb_ly = 3
  integer, parameter :: nb_hy = 4

  ! Each box contains a tag, for which the following bits are set:
  integer, parameter :: a5_bit_in_use = 1
  integer, parameter :: a5_bit_fresh = 2

  type box2_t
     integer               :: lvl
     integer               :: tag
     integer               :: ix(2)
     integer               :: parent
     integer               :: children(2,2)
     integer               :: neighbors(2,2)
     real(dp), allocatable :: cc(:, :, :)
     real(dp), allocatable :: fx(:, :, :)
     real(dp), allocatable :: fy(:, :, :)
  end type box2_t

  type level2_t
     integer, allocatable        :: ids(:)
  end type level2_t

  type a2_t
     real(dp)                  :: r_min(2)
     real(dp)                  :: dr(2, a5_max_levels)
     type(level2_t)            :: levels(a5_max_levels)
     integer                   :: n_levels
     integer                   :: box_size
     integer                   :: n_cc
     integer                   :: n_fx
     integer                   :: n_fy
     integer                   :: n_boxes
     type(box2_t), allocatable :: boxes(:)
  end type a2_t

contains

  subroutine a2_init(self, dr, n_boxes_max, box_size, n_cc, n_fx, n_fy)
    type(a2_t)        :: self
    real(dp), intent(in) :: dr(2)
    integer, intent(in)  :: n_boxes_max
    integer, intent(in)  :: box_size, n_cc, n_fx, n_fy
    integer              :: lvl, ix

    do lvl = 1, a5_max_levels
       self%dr(:, lvl) = dr_coarse * 0.5_dp**(lvl-1)
    end do

    self%n_levels = 0
    self%box_size = box_size
    self%n_cc     = n_cc
    self%n_fx     = n_fx
    self%n_fy     = n_fy
    self%n_boxes  = 0
    self%ref_func => ref_func
    allocate(self%boxes(n_boxes_max))
  end subroutine a2_init

  subroutine a5_tidy_storage(self, new_size)
    type(a2_t) :: ! TODO
  end subroutine a5_tidy_storage

  subroutine a2_set_neighbors(self)
    type(a2_t), intent(inout) :: self

    do lvl = 1, self%n_levels
       n_boxes = size(self%levels(lvl)%box_ids)
       do i = 1, n_boxes
          ! For each "fresh" box, find possible neighbors
          id = self%levels(lvl)%box_ids(i)
          if (btest(self%boxes(id)%tag, a5_bit_fresh)) then
             call set_nbs_2d(self, id)
          end if
       end do
    end do
  end subroutine a2_set_neighbors

  subroutine set_nbs_2d(self, id)
    type(a2_t), intent(inout) :: self
    integer, intent(in)          :: id
    integer                      :: d, i

    do nb = 1, 4                 ! Dimension
       if (self%boxes(id)%neighbors(nb) == a5_no_neighbor) then
          nb_id = find_nb_2d(self%boxes, id, nb)

          if (nb_id /= -1) then
             self%boxes(id)%neighbors(nb) = nb_id
             self%boxes(nb_id)%neighbors(a2_nb_rev(nb)) = id
          end if
       end if
    end do
  end subroutine set_nbs_2d

  ! Get neighbor nb of id, through its parent
  function find_nb_2d(boxes, id, nb) result(nb_id)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)        :: id, nb
    integer                    :: p_id, c_ix

    p_id = boxes(id)%parent
    c_ix = boxes(id)%

    ! Check if neighbor is in same direction as ix is (low/high). If so,
    ! use neighbor of parent
    if (a2_ch_low(c_ix, d) .eqv. a2_nb_low(nb)) &
         p_id = boxes(p_id)%neighbors(nb)

    ! The child ix of the neighbor is swapped in direction d
    d = a2_nb_dim(nb)
    nb_id = boxes(p_id)%child(a2_ch_rev(c_ix, d))
  end function find_nb_2d

  ! After boxes have been added to level 1, this stores them as a "level"
  subroutine a2_set_base(self, ix_list, nb_list)
    use mrgrnk
    type(a2_t), intent(inout) :: self
    integer, intent(in)       :: ix_list(:, :), nb_list(:, :, :)
    integer                   :: n_boxes
    integer, allocatable      :: ixs(:)
    type(a5_box_2d)              :: new_box

    self%n_levels = 1
    n_boxes       = size(ix_list, 2)
    allocate(self%levels(1)%ids(n_boxes))

    do id = 1, n_boxes
       new_box%ix        = ix_list(:, ix)
       new_box%lvl       = 1
       new_box%parent    = 0
       new_box%children  = 0
       new_box%neighbors = nb_list(:,:, ix)

       ! Set "in_use" and "fresh" tag
       new_box%tag       = 0
       new_box%tag       = ibset(new_box%tag, a5_bit_in_use)
       new_box%tag       = ibset(new_box%tag, a5_bit_fresh)

       ! Add box to storage
       self%n_boxes   = id
       self%boxes(id) = new_box
       self%levels(1)%box_ids(id) = id
       call alloc_storage_2d(self, id)
    end do
  end subroutine a2_set_base

  subroutine a2_set_levels(self)
    type(a2_t), intent(inout) :: self

    ! Lvl 1 is always set
    do lvl = 1, self%n_levels-1
       call set_child_lvl(self%levels(lvl)%ids, &
            self%levels(lvl+1)%ids, self%boxes)
    end do
  end subroutine a2_get_child_array

  ! Refine function returns integer:
  ! 0  -> derefine
  ! 1  -> refine
  ! All other values do nothing
  subroutine a2_refine(self, refine_func)
    type(a2_t), intent(inout) :: self
    procedure(a2_to_int_f) :: refine_func

    ! Set refinement flags for all boxes
    do lvl = 1, self%n_levels
       do i = 1, size(self%levels(lvl)%ids)
          id = self%levels(lvl)%ids(i)
          ref = refine_func(self%boxes(id))
          select case (ref)
          case (0)
             self%boxes(id)%tag = ibset(self%boxes(id)%tag, a5_derefine)
          case (1)
             self%boxes(id)%tag = ibset(self%boxes(id)%tag, a5_refine)
          end select
       end do
    end do

    ! Propagate refinement flags to parents
    do lvl = 1, self%n_levels
       do i = 1, size(self%levels(lvl)%ids)
          id = self%levels(lvl)%ids(i)
          
  end subroutine a2_refine

  subroutine a2_derefine(boxes, id)
    type(box2_t), intent(inout) :: boxes
    integer, intent(in) :: id

    ! Set tag to not in use
    boxes(id)%tag = ibclr(boxes(id)%tag, a5_in_use)

    ! Remove this box at the neighbors
    do nb = 1, 4
       nb_id = boxes(id)%neighbors(nb)
       if (nb_id =/ a5_no_neighbor) then
          nb_rev = a2_nb_rev(nb)
          boxes(nb_id)%neighbors(nb_rev) = a5_no_neighbor
       end if
    end do

    ! Remove this box from parent
  end subroutine a2_derefine
  
  subroutine set_child_lvl(p_ids, c_ids, boxes)
    integer, intent(in) :: p_ids(:)
    integer, allocatable, intent(inout) :: c_ids(:)
    type(box2_t), intent(in) :: boxes(:)

    ! Count 4 times the number of refined parent blocks
    n_children = 4 * count(btest(boxes(p_ids)%tag, a5_refined))
    allocate(c_ids(n_children))

    ic = 0
    do i = 1, size(p_ids)
       ip = p_ids(i)
       if (btest(boxes(ip)%tag, a5_refined)) then
          c_ids(ic+1:ic+4) = reshape(boxes(ip)%children, ...)
          ic = ic + 4
       end if
    end do
  end subroutine set_child_lvl

  subroutine a2_prolong()

  end subroutine a2_prolong

  subroutine a2_restrict()

  end subroutine a2_restrict

end module m_afivo