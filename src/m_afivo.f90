module m_afivo

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  integer, parameter :: a5_no_neighbor = 0

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
     type(morton_t), allocatable :: mortons(:)
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
    allocate(self%mortons(n_boxes_max))
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

    do d = 1, 2                 ! Dimension
       do i = 1, 2              ! Lower / higher
          if (self%boxes(id)%neighbors(d, i) == a5_no_neighbor) then
          nb_id = get_nb_id_2d(self, id, d, i)

          if (nb_id /= -1) then
             self%boxes(id)%neighbors(d, i) = nb_id
             j = reverse_ix(i)
             self%boxes(nb_id)%neighbors(d, j) = id
          end if
       end if
    end do
  end subroutine set_nbs_2d

  ! Swap 1 -> 2, 2 -> 1
  elemental function reverse_ix(ix)
    integer, intent(in) :: ix
    reverse_ix = ieor(nb_num, 3)
  end function reverse_ix

  ! Odd to 1, even to 2
  elemental function get_child_ix(ix)
    integer, intent(in) :: ix
    get_child_ix = 2 - iand(ix, 1)
  end function get_child_ix

  ! Get neighbor_id of id, in direction d, i = 1 means lower neighbor, i = 2
  ! means higher neighbor.
  function get_nb_id_2d(boxes, id, d, i) result(nb_id)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)        :: id, d, i
    integer                    :: p_id, c_ix(2)

    p_id = boxes(ix)%parent
    c_ix = get_child_ix(boxes(id)%ix)

    ! Check if neighbor is in same direction as current id is (1 = lower ix, 2 =
    ! higher ix). If so, use neighbor of parent
    if (btest(c_ix(d), 0) .eqv. btest(i, 0)) &
         p_id = boxes(p_id)%neighbors(d, i)

    ! The child ix of the neighbor is swapped in direction d
    c_ix(d) = reverse(c_ix(d))
    nb_id = boxes(p_id)%child(c_ix(1), c_ix(2))
  end function get_nb_id_2d

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
    allocate(self%levels(1)%mortons(n_boxes))

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
       call alloc_storage_2d(self, id)

       ! Set morton number for box
       call morton_from_ix2(ix, self%mortons(id))
    end do

    call mrgrnk(self%mortons(1:n_boxes), ixs)
    self%levels(1)%box_ids(:) = ixs(:)
  end subroutine a2_set_base

  subroutine a2_set_levels(self)
    type(a2_t), intent(inout) :: self

    ! Lvl 1 is always set
    do lvl = 1, self%n_levels-1
       call set_child_lvl(self%levels(lvl), self%levels(lvl+1), self%boxes)
       n_children = count(btest(self%boxes(lvl_ixs)%tag, a5_refined))
       n_parent = size(self%levels(lvl)%ids)
       do ip = 1, n_parent
          id = self%levels(lvl)%ids(ip)
          
       end do
    end do
  end subroutine a2_get_child_array

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

    ! TODO: sort immediately?
    allocate(c_lvl%mortons(n_children))
  end subroutine set_child_lvl

  subroutine a2_prolong()

  end subroutine a2_prolong

  subroutine a2_restrict()

  end subroutine a2_restrict

  function find_box_2d(lvl_2d, morton) result(id)
    type(level2_t), intent(in) :: lvl_2d
    type(morton_t), intent(in)   :: morton
    integer                      :: id, ix

    ix = morton_bsearch(lvl_2d%mortons, morton)
    if (ix /= -1) then
       id = lvl_2d%ids(ix)
    else
       id = -1
    end if
  end function find_box_2d
  
end module m_afivo