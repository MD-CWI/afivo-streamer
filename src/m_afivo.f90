module m_afivo

  implicit none
  public

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: a5_max_levels = 20

  ! Tags that can be set for a block
  integer, parameter :: a5_do_ref = 1
  integer, parameter :: a5_rm_ref = 2
  integer, parameter :: a5_kp_ref = 3
  integer, parameter :: a5_rm_children = 4
  integer, parameter :: a5_no_box = 0

  ! Corners
  integer, parameter :: cn_lxly = 1
  integer, parameter :: cn_hxly = 2
  integer, parameter :: cn_lxhy = 3
  integer, parameter :: cn_hxhy = 4
  integer, parameter :: a2_cn_nbs(2, 4) = reshape([1,3,3,2,4,1,2,4], [2,4])

  integer, parameter :: nb_lx = 1
  integer, parameter :: nb_hx = 2
  integer, parameter :: nb_ly = 3
  integer, parameter :: nb_hy = 4

  ! Children (same as corners)
  integer, parameter :: ch_lxly = 1
  integer, parameter :: ch_hxly = 2
  integer, parameter :: ch_lxhy = 3
  integer, parameter :: ch_hxhy = 4
  integer, parameter :: a2_ch_dix(2, 4) = reshape([0,0,1,0,0,1,1,1], [2,4])
  integer, parameter :: a2_ch_rev(4, 2) = reshape([2,1,4,3,3,4,1,2], [4,2])
  logical, parameter :: a2_ch_low(4, 2) = reshape([ .true., .false., .true., &
       .false., .true., .true., .false., .false.], [4,2])

  ! Neighbor topology information (todo: add documentation)
  logical, parameter :: a2_nb_low(4) = [.true., .false., .true., .false.]
  integer, parameter :: a2_nb_rev(4) = [2, 1, 4, 3]
  integer, parameter :: a2_nb_dim(4) = [1, 1, 2, 2]

  ! Each box contains a tag, for which the following bits are set:
  integer, parameter :: a5_bit_in_use = 1
  integer, parameter :: a5_bit_new_children = 2

  type box2_cfg_t
     integer  :: n_cell
     integer  :: n_node
     integer  :: n_cc
     integer  :: n_flux
     real(dp) :: dr(2, a5_max_levels)
     real(dp) :: dbr(2, a5_max_levels)
     real(dp) :: r_min(2)
  end type box2_cfg_t

  type box2_t
     integer                   :: lvl
     integer                   :: tag
     integer                   :: ix(2)
     integer                   :: parent
     integer                   :: children(4)
     integer                   :: neighbors(4)
     type(box2_cfg_t), pointer :: cfg
     real(dp), allocatable     :: cc(:, :, :)
     real(dp), allocatable     :: fx(:, :, :)
     real(dp), allocatable     :: fy(:, :, :)
  end type box2_t

  type level2_t
     integer, allocatable :: ids(:)
  end type level2_t

  type a2_t
     type(level2_t)            :: levels(a5_max_levels)
     type(box2_cfg_t), pointer :: cfg => null()
     integer                   :: n_boxes
     type(box2_t), allocatable :: boxes(:)
  end type a2_t

  abstract interface
     integer function a2_to_int_f(box)
       import
       type(box2_t), intent(in) :: box
     end function a2_to_int_f

     subroutine a2_subr(box)
       import
       type(box2_t), intent(inout) :: box
     end subroutine a2_subr

     subroutine a2_subr_boxes(boxes, id)
       import
       type(box2_t), intent(inout) :: boxes(:)
       integer, intent(in)         :: id
     end subroutine a2_subr_boxes

     subroutine a2_subr_box_int(boxes, id, i)
       import
       type(box2_t), intent(inout) :: boxes(:)
       integer, intent(in)         :: id, i
     end subroutine a2_subr_box_int
  end interface

contains

  subroutine a2_init(tree, n_boxes_max, n_cell, n_cc, n_flux, dr, r_min)
    type(a2_t), intent(out) :: tree
    integer, intent(in)     :: n_boxes_max
    integer, intent(in)     :: n_cell, n_cc, n_flux
    real(dp), intent(in)    :: dr(2), r_min(2)
    integer                 :: lvl

    if (n_cell < 2)       stop "a2_init: n_cell should be >= 2"
    if (btest(n_cell, 0)) stop "a2_init: n_cell should be even"
    if (n_cc <= 0)        stop "a2_init: n_cc should be > 0"
    if (n_boxes_max <= 0) stop "a2_init: n_boxes_max should be > 0"

    allocate(tree%cfg)
    tree%cfg%n_cell = n_cell
    tree%cfg%n_node = n_cell + 1
    tree%cfg%n_cc   = n_cc
    tree%cfg%n_flux = n_flux
    tree%cfg%r_min  = r_min
    tree%n_boxes    = 0

    allocate(tree%boxes(n_boxes_max))
    tree%boxes(:)%tag = 0

    do lvl = 1, a5_max_levels
       allocate(tree%levels(lvl)%ids(0))
       tree%cfg%dr(:, lvl) = dr * 0.5_dp**(lvl-1)
       tree%cfg%dbr(:, lvl) = dr * n_cell * 0.5_dp**(lvl-1)
    end do
  end subroutine a2_init

  subroutine a2_destroy(tree)
    type(a2_t), intent(inout) :: tree
    integer                   :: lvl
    deallocate(tree%cfg)
    deallocate(tree%boxes)
    do lvl = 1, a5_max_levels
       deallocate(tree%levels(lvl)%ids)
    end do
    tree%n_boxes = 0
  end subroutine a2_destroy

  integer function a2_n_levels(tree)
    type(a2_t), intent(in) :: tree
    integer                :: lvl
    do lvl = 1, a5_max_levels-1
       if (size(tree%levels(lvl)%ids) == 0) exit
    end do
    a2_n_levels = lvl
  end function a2_n_levels

  subroutine a2_loop_box(tree, my_procedure)
    type(a2_t), intent(inout) :: tree
    procedure(a2_subr)        :: my_procedure
    integer                   :: lvl, i, id

    do lvl = 1, a2_n_levels(tree)
       do i = 1, size(tree%levels(lvl)%ids)
          id = tree%levels(lvl)%ids(i)
          call my_procedure(tree%boxes(id))
       end do
    end do
  end subroutine a2_loop_box

  subroutine a2_loop_boxes(tree, my_procedure)
    type(a2_t), intent(inout) :: tree
    procedure(a2_subr_boxes)  :: my_procedure
    integer                   :: lvl, i, id

    do lvl = 1, a2_n_levels(tree)
       do i = 1, size(tree%levels(lvl)%ids)
          id = tree%levels(lvl)%ids(i)
          call my_procedure(tree%boxes, id)
       end do
    end do
  end subroutine a2_loop_boxes

  subroutine a2_tidy_up(tree, frac_max, frac_goal, n_clean_min, reorder)
    use m_morton
    type(a2_t), intent(inout)      :: tree
    real(dp), intent(in)           :: frac_max, frac_goal
    integer, intent(in)            :: n_clean_min
    logical, intent(in)            :: reorder
    real(dp)                       :: frac_in_use
    logical :: only_reorder
    integer                        :: n, lvl, id, old_size, new_size, n_clean
    integer                        :: n_boxes, n_used, n_stored, n_used_lvl
    integer, allocatable           :: ixs_sort(:), ixs_map(:)
    type(box2_t), allocatable      :: boxes_cpy(:)
    integer(morton_k), allocatable :: mortons(:)

    if (frac_goal > frac_max) stop "a2_tidy_up: need frac_goal < frac_max"
    if (frac_max >= 1.0_dp)   stop "a2_tidy_up: need frac_max < 1"
    if (n_clean_min < 1)      stop "a2_tidy_up: need n_clean_min > 0"

    n_boxes     = tree%n_boxes
    n_used      = count(btest(tree%boxes(1:n_boxes)%tag, a5_bit_in_use))
    old_size    = size(tree%boxes)
    frac_in_use = n_used / real(old_size, dp)
    n_clean     = nint((frac_goal - frac_in_use) * old_size)

    if (frac_in_use > frac_max .or. (frac_in_use < frac_goal .and. &
         n_clean > n_clean_min)) then
       new_size = max(1, nint(n_used/frac_goal))
    else
       new_size = old_size
    end if

    if (new_size /= old_size .or. reorder) then
       only_reorder = (new_size == old_size)
       print *, "a2_tidy_up:", old_size, new_size, only_reorder

       if (only_reorder) then
          allocate(boxes_cpy(n_used))  ! Need just enough space
       else
          allocate(boxes_cpy(new_size))
       end if

       allocate(ixs_map(0:n_boxes))
       ixs_map(0)       = 0
       n_stored         = 0

       do lvl = 1, a2_n_levels(tree)
          n_used_lvl = size(tree%levels(lvl)%ids)
          allocate(mortons(n_used_lvl))
          allocate(ixs_sort(n_used_lvl))

          do n = 1, n_used_lvl
             id = tree%levels(lvl)%ids(n)
             ! Note the -1, since our indices start at 1
             mortons(n) = morton_from_ix2(tree%boxes(id)%ix-1)
          end do

          call morton_rank(mortons, ixs_sort)
          tree%levels(lvl)%ids = tree%levels(lvl)%ids(ixs_sort)

          do n = 1, n_used_lvl
             id = tree%levels(lvl)%ids(n)
             boxes_cpy(n_stored + n) = tree%boxes(id)
             ixs_map(tree%levels(lvl)%ids(n)) = n_stored + n
          end do

          tree%levels(lvl)%ids = [(n_stored+n, n=1,n_used_lvl)]
          n_stored = n_stored + n_used_lvl
          deallocate(mortons)
          deallocate(ixs_sort)
       end do

       ! Update id's to new indices
       do n = 1, n_used
          boxes_cpy(n)%parent = ixs_map(boxes_cpy(n)%parent)
          boxes_cpy(n)%children = ixs_map(boxes_cpy(n)%children)
          boxes_cpy(n)%neighbors = ixs_map(boxes_cpy(n)%neighbors)
       end do

       if (only_reorder) then
          tree%boxes(1:n_used) = boxes_cpy
          tree%boxes(n_used+1:n_boxes)%tag = 0
          do n = n_used+1, n_boxes
             deallocate(tree%boxes(n)%cc)
             deallocate(tree%boxes(n)%fx)
             deallocate(tree%boxes(n)%fy)
          end do
          tree%n_boxes = n_used
       else
          call move_alloc(boxes_cpy, tree%boxes)
       end if

       tree%n_boxes = n_used
    end if
  end subroutine a2_tidy_up

  subroutine a2_set_new_child_neighbors(boxes, ids)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:)
    integer                     :: i, id, i_c

    do i = 1, size(ids)
       ! For each newly refined box, find children neighbors
       id = ids(i)
       if (btest(boxes(id)%tag, a5_bit_new_children)) then
          do i_c = 1, 4
             call set_nbs_2d(boxes, boxes(id)%children(i_c))
          end do
       end if
    end do
  end subroutine a2_set_new_child_neighbors

  subroutine set_nbs_2d(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)          :: id
    integer                      :: nb, nb_id

    do nb = 1, 4                 ! Dimension
       if (boxes(id)%neighbors(nb) == a5_no_box) then
          nb_id = find_nb_2d(boxes, id, nb)
          if (nb_id /= a5_no_box) then
             boxes(id)%neighbors(nb) = nb_id
             boxes(nb_id)%neighbors(a2_nb_rev(nb)) = id
          end if
       end if
    end do
  end subroutine set_nbs_2d

  integer function a2_ix_to_cix(ix)
    integer, intent(in) :: ix(2)
    ! Second index odd: -2, first ix odd: -1
    a2_ix_to_cix = 4 - 2 * iand(ix(2), 1) - iand(ix(1), 1)
  end function a2_ix_to_cix

  ! Get neighbor nb of id, through its parent
  function find_nb_2d(boxes, id, nb) result(nb_id)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id, nb
    integer                  :: nb_id, p_id, c_ix, d, old_pid

    p_id = boxes(id)%parent
    old_pid = p_id
    c_ix = a2_ix_to_cix(boxes(id)%ix)
    d = a2_nb_dim(nb)

    ! Check if neighbor is in same direction as ix is (low/high). If so,
    ! use neighbor of parent
    if (a2_ch_low(c_ix, d) .eqv. a2_nb_low(nb)) &
         p_id = boxes(p_id)%neighbors(nb)

    ! The child ix of the neighbor is reversed in direction d
    nb_id = boxes(p_id)%children(a2_ch_rev(c_ix, d))
  end function find_nb_2d

  ! After boxes have been added to level 1, this stores them as a "level"
  subroutine a2_set_base(tree, ix_list, nb_list)
    type(a2_t), intent(inout) :: tree
    integer, intent(in)       :: ix_list(:, :), nb_list(:, :)
    integer                   :: n_boxes, i, id

    if (any(ix_list < 1)) stop "a2_set_base: need ix_list > 0"

    n_boxes = size(ix_list, 2)
    deallocate(tree%levels(1)%ids)
    allocate(tree%levels(1)%ids(n_boxes))
    call a2_get_free_ids(tree, tree%levels(1)%ids)

    do i = 1, n_boxes
       id                       = tree%levels(1)%ids(i)
       tree%boxes(id)%ix        = ix_list(:, i)
       tree%boxes(id)%lvl       = 1
       tree%boxes(id)%parent    = 0
       tree%boxes(id)%children  = 0
       tree%boxes(id)%neighbors = nb_list(:, i)
       tree%boxes(id)%cfg       => tree%cfg
       tree%boxes(id)%tag       = ibset(0, a5_bit_in_use)
    end do
  end subroutine a2_set_base

  ! On input, tree should be balanced. On output, tree is still balanced, and
  ! its refinement is updated (with at most one level per call).
  ! Sets the following bits: a5_bit_uninitialized, a5_bit_new_children
  subroutine a2_adjust_refinement(tree, ref_func)
    type(a2_t), intent(inout) :: tree
    procedure(a2_to_int_f)    :: ref_func
    integer                   :: lvl, id, i, c_ids(4), n_boxes_start
    integer :: n_boxes_req
    integer, allocatable      :: ref_vals(:)
    type(box2_t), allocatable :: boxes_cpy(:)

    n_boxes_start = tree%n_boxes
    allocate(ref_vals(n_boxes_start))
    call a2_set_ref_vals(tree%levels, tree%boxes, ref_vals, ref_func)

    ! Check whether there is enough free space, otherwise extend the list
    n_boxes_req = n_boxes_start + 4 * count(ref_vals == a5_do_ref)
    if (n_boxes_req > size(tree%boxes)) then
       print *, "Resizing box storage for refinement", n_boxes_req
       allocate(boxes_cpy(n_boxes_req))
       boxes_cpy(1:n_boxes_start)      = tree%boxes(1:n_boxes_start)
       boxes_cpy(n_boxes_start+1:)%tag = 0 ! Set to not in use
       call move_alloc(boxes_cpy, tree%boxes)
    end if

    do lvl = 1, a2_n_levels(tree)
       do i = 1, size(tree%levels(lvl)%ids)
          id = tree%levels(lvl)%ids(i)

          if (id > n_boxes_start) then
             cycle              ! This is a newly added box
          else if (ref_vals(id) == a5_do_ref) then
             ! Add children. First need to get 4 free id's
             call a2_get_free_ids(tree, c_ids)
             call a2_add_children(tree%boxes, id, c_ids)
          else if (ref_vals(id) == a5_rm_children) then
             ! Remove children
             call a2_remove_children(tree%boxes, id)
          end if
       end do

       ! Set the next level
       if (lvl < a5_max_levels) then
          deallocate(tree%levels(lvl+1)%ids)
          call set_child_lvl(tree%levels(lvl)%ids, &
               tree%levels(lvl+1)%ids, tree%boxes)
          ! Set connectivity
          call a2_set_new_child_neighbors(tree%boxes, tree%levels(lvl)%ids)
       end if
    end do
  end subroutine a2_adjust_refinement

  subroutine a2_get_free_ids(tree, ids)
    type(a2_t), intent(inout) :: tree
    integer, intent(out)      :: ids(:)
    integer                   :: i, n_ids, id, nc

    n_ids = size(ids)
    ! Critical part
    ids(1)       = tree%n_boxes + 1
    tree%n_boxes = tree%n_boxes + n_ids
    ! End critical
    do i = 2, n_ids
       ids(i) = ids(1) + i - 1
    end do

    ! Allocate storage
    nc = tree%cfg%n_cell
    do i = 1, n_ids
       id = ids(i)
       allocate(tree%boxes(id)%cc(0:nc+1, 0:nc+1, tree%cfg%n_cc))
       allocate(tree%boxes(id)%fx(nc+1, nc, tree%cfg%n_flux))
       allocate(tree%boxes(id)%fy(nc, nc+1, tree%cfg%n_flux))
    end do
  end subroutine a2_get_free_ids

  subroutine a2_set_ref_vals(levels, boxes, ref_vals, ref_func)
    type(level2_t), intent(in)  :: levels(:)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(inout)      :: ref_vals(:)
    procedure(a2_to_int_f)      :: ref_func

    integer                     :: n_levels
    integer                     :: lvl, i, id, c_ids(4)
    integer                     :: nb, p_id, nb_id, p_nb_id

    n_levels = size(levels)
    ref_vals = a5_kp_ref        ! Used indices are overwritten

    ! Set refinement flags for all boxes using ref_func
    do lvl = 1, n_levels
       do i = 1, size(levels(lvl)%ids)
          id           = levels(lvl)%ids(i)
          ref_vals(id) = ref_func(boxes(id))
       end do
    end do

    ! Cannot derefine lvl 1
    do i = 1, size(levels(1)%ids)
       id = levels(1)%ids(i)
       if (ref_vals(id) == a5_rm_ref) ref_vals(id) = a5_kp_ref
    end do

    ! Cannot refine beyond max level
    do i = 1, size(levels(a5_max_levels)%ids)
       id = levels(a5_max_levels)%ids(i)
       if (ref_vals(id) == a5_do_ref) ref_vals(id) = a5_kp_ref
    end do

    ! Ensure 2-1 balance.
    do lvl = n_levels, 2, -1
       do i = 1, size(levels(lvl)%ids)
          id = levels(lvl)%ids(i)
          ! We only need check leaf boxes
          if (a2_has_children(boxes(id))) cycle

          if (ref_vals(id) == a5_do_ref) then
             ! Ensure we will have the necessary neighbors
             do nb = 1, 4
                nb_id = boxes(id)%neighbors(nb)
                if (nb_id == a5_no_box) then
                   ! Mark the parent containing neighbor for refinement
                   p_id = boxes(id)%parent
                   p_nb_id = boxes(p_id)%neighbors(nb)
                   ref_vals(p_nb_id) = a5_do_ref
                end if
             end do
          else if (ref_vals(id) == a5_rm_ref) then
             ! Ensure we do not remove a required neighbor
             do nb = 1, 4
                nb_id = boxes(id)%neighbors(nb)
                if (nb_id /= a5_no_box) then
                   if (a2_has_children(boxes(nb_id)) .or. &
                        ref_vals(nb_id) == a5_do_ref) then
                      ref_vals(id) = a5_kp_ref
                   end if
                end if
             end do
          end if
       end do
    end do

    ! Make the (de)refinement flags consistent for blocks with children.
    do lvl = n_levels-1, 1, -1
       do i = 1, size(levels(lvl)%ids)
          id = levels(lvl)%ids(i)
          if (.not. a2_has_children(boxes(id))) cycle

          ! Can only remove children if they are all marked for
          ! derefinement, and the box ittree not for refinement.
          c_ids = boxes(id)%children
          if (all(ref_vals(c_ids) == a5_rm_ref) .and. &
               ref_vals(id) /= a5_do_ref) then
             ! The children are removed
             ref_vals(id) = a5_rm_children
          else
             ! The children cannot be removed
             ref_vals(id) = a5_kp_ref
             where (ref_vals(c_ids) == a5_rm_ref)
                ref_vals(c_ids) = a5_kp_ref
             end where
          end if
       end do
    end do

  end subroutine a2_set_ref_vals

  subroutine a2_remove_children(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    integer                     :: ic, c_id, nb_id, nb_rev, nb

    do ic = 1, 4
       c_id            = boxes(id)%children(ic)
       boxes(c_id)%tag = ibclr(boxes(c_id)%tag, a5_bit_in_use)

       do nb = 1, 4             ! Remove from neighbors
          nb_id = boxes(c_id)%neighbors(nb)
          if (nb_id /= a5_no_box) then
             nb_rev = a2_nb_rev(nb)
             boxes(nb_id)%neighbors(nb_rev) = a5_no_box
          end if
       end do
    end do

    boxes(id)%children = a5_no_box
  end subroutine a2_remove_children

  subroutine a2_add_children(boxes, id, c_ids)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    integer, intent(in)         :: c_ids(4)
    integer                     :: i, c_id, c_ix_base(2)

    boxes(id)%children = c_ids
    c_ix_base          = 2 * boxes(id)%ix - 1
    boxes(id)%tag      = ibset(boxes(id)%tag, a5_bit_new_children)

    do i = 1, 4
       c_id                  = c_ids(i)
       boxes(c_id)%ix        = c_ix_base + a2_ch_dix(:,i)
       boxes(c_id)%lvl       = boxes(id)%lvl+1
       boxes(c_id)%parent    = id
       boxes(c_id)%children  = a5_no_box
       boxes(c_id)%neighbors = a5_no_box
       boxes(c_id)%cfg       => boxes(id)%cfg

       ! Set "in_use" tag
       boxes(c_id)%tag       = ibset(0, a5_bit_in_use)
    end do
  end subroutine a2_add_children

  subroutine set_child_lvl(p_ids, c_ids, boxes)
    integer, intent(in)                 :: p_ids(:)
    integer, allocatable, intent(inout) :: c_ids(:)
    type(box2_t), intent(in)            :: boxes(:)
    integer                             :: i, ip, i_c, n_children

    ! Count 4 times the number of refined parent blocks
    n_children = 4 * count(a2_has_children(boxes(p_ids)))
    allocate(c_ids(n_children))

    i_c = 0
    do i = 1, size(p_ids)
       ip = p_ids(i)
       if (a2_has_children(boxes(ip))) then
          c_ids(i_c+1:i_c+4) = boxes(ip)%children
          i_c = i_c + 4
       end if
    end do
  end subroutine set_child_lvl

  elemental logical function a2_has_children(box)
    type(box2_t), intent(in) :: box
    a2_has_children = (box%children(1) /= a5_no_box)
  end function a2_has_children

  pure function a2_dr(box) result(dr)
    type(box2_t), intent(in) :: box
    real(dp) :: dr(2)
    dr = box%cfg%dr(:, box%lvl)
  end function a2_dr

  pure function a2_r_min(box) result(r_min)
    type(box2_t), intent(in) :: box
    real(dp) :: r_min(2)
    r_min = box%cfg%r_min + (box%ix-1) * box%cfg%dbr(:, box%lvl)
  end function a2_r_min

  ! Location of cell center with index cc_ix
  pure function a2_r_cc(box, cc_ix) result(r)
    type(box2_t), intent(in) :: box
    integer, intent(in)      :: cc_ix(2)
    real(dp)                 :: r(2)
    r = a2_r_min(box) + (cc_ix-0.5_dp) * box%cfg%dr(:, box%lvl)
  end function a2_r_cc

  ! Location of node with index nd_ix
  pure function a2_r_node(box, nd_ix) result(r)
    type(box2_t), intent(in) :: box
    integer, intent(in)      :: nd_ix(2)
    real(dp)                 :: r(2)
    r = a2_r_min(box) + (nd_ix-1) * box%cfg%dr(:, box%lvl)
  end function a2_r_node

  elemental integer function l2i(my_bool)
    logical, intent(in) :: my_bool
    if (my_bool) then
       l2i = 1
    else
       l2i = 0
    end if
  end function l2i

  ! Injection to children.
  subroutine a2_prolong0_from(boxes, id, v_ixs, fill_gc)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, v_ixs(:)
    logical, intent(in)         :: fill_gc
    integer                     :: nc, i_c, c_id, ix_offset(2)
    integer                     :: i, j, i_c1, j_c1, iv, v

    nc = boxes(id)%cfg%n_cell
    do i_c = 1, 4
       c_id = boxes(id)%children(i_c)
       ix_offset = a2_ch_dix(:, i_c) * ishft(nc, -1) ! Offset child

       do v = 1, size(v_ixs)
          iv = v_ixs(v)
          do j = 1-l2i(fill_gc), nc+l2i(fill_gc)
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             do i = 1-l2i(fill_gc), nc+l2i(fill_gc)
                i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
                boxes(c_id)%cc(i, j, iv) = boxes(id)%cc(i_c1, j_c1, iv)
             end do
          end do
       end do
    end do
  end subroutine a2_prolong0_from

  ! Bilinear prolongation to children. Uses ghost cells and corners.
  subroutine a2_prolong1_from(boxes, id, v_ixs, fill_gc)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, v_ixs(:)
    logical, intent(in)         :: fill_gc
    real(dp), parameter         :: f1=1/16.0_dp, f3=3/16.0_dp, f9=9/16.0_dp
    integer                     :: v, iv, nc, i_c, c_id, ix_offset(2)
    integer                     :: i, j, i_c1, i_c2, j_c1, j_c2

    nc = boxes(id)%cfg%n_cell
    do i_c = 1, 4
       c_id = boxes(id)%children(i_c)

       ! Offset of child w.r.t. parent
       ix_offset = a2_ch_dix(:, i_c) * ishft(nc, -1)

       ! In these loops, we calculate the closest coarse index (_c1), and the
       ! one-but-closest (_c2). The fine cell lies in between.
       do v = 1, size(v_ixs)
          iv = v_ixs(v)
          do j = 1-l2i(fill_gc), nc+l2i(fill_gc)
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
             do i = 1-l2i(fill_gc), nc+l2i(fill_gc)
                i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
                i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1

                boxes(c_id)%cc(i, j, iv) = &
                     f9 * boxes(id)%cc(i_c1, j_c1, iv) + &
                     f3 * boxes(id)%cc(i_c2, j_c1, iv) + &
                     f3 * boxes(id)%cc(i_c1, j_c2, iv) + &
                     f1 * boxes(id)%cc(i_c2, j_c2, iv)
             end do
          end do
       end do
    end do
  end subroutine a2_prolong1_from

  ! Partial prolongation from parent using injection.
  subroutine a2_prolong0_to(boxes, id, i_ixs, j_ixs, v_ixs)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, i_ixs(:), j_ixs(:), v_ixs(:)
    integer                     :: v, iv, nc, p_id, ix_offset(2)
    integer                     :: ii, jj, i, j, i_c1, j_c1

    nc        = boxes(id)%cfg%n_cell
    p_id      = boxes(id)%parent
    ! Offset of child w.r.t. parent
    ix_offset = (boxes(id)%ix - 2*boxes(p_id)%ix + 1) * ishft(nc, -1)

    do v = 1, size(v_ixs)
       iv = v_ixs(v)
       do jj = 1, size(j_ixs)
          j = j_ixs(jj)
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          do ii = 1, size(i_ixs)
             i = i_ixs(ii)
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             boxes(id)%cc(i, j, iv) = boxes(p_id)%cc(i_c1, j_c1, iv)
          end do
       end do
    end do
  end subroutine a2_prolong0_to

  ! Partial bilinear prolongation from parent.
  subroutine a2_prolong1_to(boxes, id, i_ixs, j_ixs, v_ixs)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, i_ixs(:), j_ixs(:), v_ixs(:)
    real(dp), parameter         :: f1=1/16.0_dp, f3=3/16.0_dp, f9=9/16.0_dp
    integer                     :: v, iv, nc, p_id, ix_offset(2)
    integer                     :: ii, jj, i, j, i_c1, i_c2, j_c1, j_c2

    nc        = boxes(id)%cfg%n_cell
    p_id      = boxes(id)%parent
    ! Offset of child w.r.t. parent
    ix_offset = (boxes(id)%ix - 2 * boxes(p_id)%ix + 1) * ishft(nc, -1)

    ! In these loops, we calculate the closest coarse index (i_c1, j_c1), and
    ! the one-but-closest (i_c2, j_c2). The fine cell lies in between.
    do v = 1, size(v_ixs)
       iv = v_ixs(v)
       do jj = 1, size(j_ixs)
          j = j_ixs(jj)
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          do ii = 1, size(i_ixs)
             i = i_ixs(ii)
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             boxes(id)%cc(i, j, iv) = &
                  f9 * boxes(p_id)%cc(i_c1, j_c1, iv) + &
                  f3 * boxes(p_id)%cc(i_c2, j_c1, iv) + &
                  f3 * boxes(p_id)%cc(i_c1, j_c2, iv) + &
                  f1 * boxes(p_id)%cc(i_c2, j_c2, iv)
          end do
       end do
    end do
  end subroutine a2_prolong1_to

  ! Conservative restriction. 4 fine cells to one parent cell
  subroutine a2_restrict_to(boxes, id, v_ixs)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, v_ixs(:)
    integer                     :: nc, i_c, c_id, ix_offset(2)
    integer                     :: v, iv, i, j, i_c1, j_c1

    nc = boxes(id)%cfg%n_cell
    do i_c = 1, 4
       c_id = boxes(id)%children(i_c)
       ! Offset of child w.r.t. parent
       ix_offset = a2_ch_dix(:, i_c) * ishft(nc, -1)

       do v = 1, size(v_ixs)
          iv = v_ixs(v)
          do j = 1, nc, 2
             j_c1 = ix_offset(2) + ishft(j+1, -1)  ! (j+1)/2
             do i = 1, nc, 2
                i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2

                boxes(id)%cc(i_c1, j_c1, iv) = 0.25_dp * (&
                     boxes(c_id)%cc(i,   j,   iv) + &
                     boxes(c_id)%cc(i+1, j,   iv) + &
                     boxes(c_id)%cc(i,   j+1, iv) + &
                     boxes(c_id)%cc(i+1, j+1, iv) )
             end do
          end do
       end do
    end do
  end subroutine a2_restrict_to

  subroutine a2_gc_sides(tree, side_subr)
    type(a2_t), intent(inout)  :: tree
    procedure(a2_subr_box_int) :: side_subr
    integer                    :: lvl, i, id

    do lvl = 1, a2_n_levels(tree)
       do i = 1, size(tree%levels(lvl)%ids)
          id = tree%levels(lvl)%ids(i)
          call a2_gc_box_sides(tree%boxes, id, side_subr)
       end do
    end do
  end subroutine a2_gc_sides

  subroutine a2_sides_from_parent(boxes, id, nb)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb
    integer                     :: n, nc, nv

    nc = boxes(id)%cfg%n_cell
    nv = boxes(id)%cfg%n_cc

    select case (nb)
    case (nb_lx)
       call a2_prolong1_to(boxes, id, [0], [(n, n=1,nc)], [(n, n=1,nv)])
    case (nb_hx)
       call a2_prolong1_to(boxes, id, [nc+1], [(n, n=1,nc)], [(n, n=1,nv)])
    case (nb_ly)
       call a2_prolong1_to(boxes, id, [(n, n=1,nc)], [0], [(n, n=1,nv)])
    case (nb_hy)
       call a2_prolong1_to(boxes, id, [(n, n=1,nc)], [nc+1], [(n, n=1,nv)])
    end select
  end subroutine a2_sides_from_parent

  subroutine a2_gc_box_sides(boxes, id, side_subr)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    procedure(a2_subr_box_int)  :: side_subr
    integer                     :: nb

    do nb = 1, 4
       if (boxes(id)%neighbors(nb) > a5_no_box) then
          call a2_gc_side_from_nb(boxes, id, nb)
       else
          call side_subr(boxes, id, nb)
       end if
    end do
  end subroutine a2_gc_box_sides

  subroutine a2_gc_side_from_nb(boxes, id, nb)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb
    integer                     :: nc, nb_id

    nc    = boxes(id)%cfg%n_cell
    nb_id = boxes(id)%neighbors(nb)

    select case (nb)
    case (nb_lx)
       boxes(id)%cc(0, 1:nc, :)    = boxes(nb_id)%cc(nc, 1:nc, :)
    case (nb_hx)
       boxes(id)%cc(nc+1, 1:nc, :) = boxes(nb_id)%cc(1, 1:nc, :)
    case (nb_ly)
       boxes(id)%cc(1:nc, 0, :)    = boxes(nb_id)%cc(1:nc, nc, :)
    case (nb_hy)
       boxes(id)%cc(1:nc, nc+1, :) = boxes(nb_id)%cc(1:nc, 1, :)
    end select
  end subroutine a2_gc_side_from_nb

  subroutine a2_gc_corners(tree, corner_subr)
    type(a2_t), intent(inout)  :: tree
    procedure(a2_subr_box_int) :: corner_subr
    integer                    :: lvl, i, id
    do lvl = 1, a2_n_levels(tree)
       do i = 1, size(tree%levels(lvl)%ids)
          id = tree%levels(lvl)%ids(i)
          call a2_gc_box_corners(tree%boxes, id, corner_subr)
       end do
    end do
  end subroutine a2_gc_corners

  subroutine a2_corners_from_parent(boxes, id, cn)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, cn
    integer                     :: n, nc, nv

    nc = boxes(id)%cfg%n_cell
    nv = boxes(id)%cfg%n_cc

    select case (cn)
    case (cn_lxly)
       call a2_prolong1_to(boxes, id, [0], [0], [(n, n=1,nv)])
    case (cn_hxly)
       call a2_prolong1_to(boxes, id, [nc+1], [0], [(n, n=1,nv)])
    case (cn_lxhy)
       call a2_prolong1_to(boxes, id, [0], [nc+1], [(n, n=1,nv)])
    case (cn_hxhy)
       call a2_prolong1_to(boxes, id, [nc+1], [nc+1], [(n, n=1,nv)])
    end select
  end subroutine a2_corners_from_parent

  subroutine a2_gc_box_corners(boxes, id, corner_subr)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    procedure(a2_subr_box_int)  :: corner_subr
    integer                     :: cn, nbs(2)

    do cn = 1, 4
       nbs = a2_cn_nbs(:, cn)
       if (boxes(id)%neighbors(nbs(1)) > a5_no_box) then
          call a2_gc_corner_from_nb(boxes, id, cn, nbs(1))
       else if (boxes(id)%neighbors(nbs(2)) > a5_no_box) then
          call a2_gc_corner_from_nb(boxes, id, cn, nbs(2))
       else
          call corner_subr(boxes, id, cn)
       end if
    end do
  end subroutine a2_gc_box_corners

  subroutine a2_gc_corner_from_nb(boxes, id, cn, nb)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, cn, nb
    integer                     :: nc, nb_id

    nb_id = boxes(id)%neighbors(nb)
    nc    = boxes(id)%cfg%n_cell

    select case (cn)
    case (cn_lxly)
       if (nb == nb_lx) then
          boxes(id)%cc(0, 0, :) = boxes(nb_id)%cc(nc, 0, :)
       else                     ! nb_ly
          boxes(id)%cc(0, 0, :) = boxes(nb_id)%cc(0, nc, :)
       end if
    case (cn_hxly)
       if (nb == nb_hx) then
          boxes(id)%cc(nc+1, 0, :) = boxes(nb_id)%cc(1, 0, :)
       else                     ! nb_ly
          boxes(id)%cc(nc+1, 0, :) = boxes(nb_id)%cc(nc+1, nc, :)
       end if
    case (cn_lxhy)
       if (nb == nb_lx) then
          boxes(id)%cc(0, nc+1, :) = boxes(nb_id)%cc(nc, nc+1, :)
       else                     ! nb_hy
          boxes(id)%cc(0, nc+1, :) = boxes(nb_id)%cc(0, 1, :)
       end if
    case (cn_hxhy)
       if (nb == nb_hx) then
          boxes(id)%cc(nc+1, nc+1, :) = boxes(nb_id)%cc(1, nc+1, :)
       else                     ! nb_hy
          boxes(id)%cc(nc+1, nc+1, :) = boxes(nb_id)%cc(nc+1, 1, :)
       end if
    end select
  end subroutine a2_gc_corner_from_nb

  subroutine a2_consisent_fluxes(tree, f_ixs)
    type(a2_t), intent(inout) :: tree
    integer, intent(in) :: f_ixs(:)
    integer :: lvl, i, id, nb, nb_id

    do lvl = a2_n_levels(tree)-1, 1, -1
       do i = 1, size(tree%levels(lvl)%ids)
          id = tree%levels(lvl)%ids(i)
          if (.not. a2_has_children(tree%boxes(id))) cycle
          do nb = 1, 4
             nb_id = tree%boxes(id)%neighbors(nb)
             if (nb_id <= a5_no_box) cycle
             call a2_flux_from_children(tree%boxes, id, nb, f_ixs)
          end do
       end do
    end do
  end subroutine a2_consisent_fluxes

  ! The neighbor nb has no children and id does, so get flux from children for
  ! consisency at refinement boundary. TODO: TEST
  subroutine a2_flux_from_children(boxes, id, nb, f_ixs)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, f_ixs(:)
    integer                     :: nc, nch, c_id

    nc  = boxes(id)%cfg%n_cc
    nch = ishft(nc, -1) ! nc/2

    select case (nb)
    case (nb_lx)
       c_id = boxes(id)%children(ch_lxly)
       boxes(id)%fx(1, 1:nch, f_ixs) = 0.5_dp * ( &
            boxes(c_id)%fx(1, 1:nc:2, f_ixs) + &
            boxes(c_id)%fx(1, 2:nc:2, f_ixs))
       c_id = boxes(id)%children(ch_lxhy)
       boxes(id)%fx(1, nch+1:, f_ixs) = 0.5_dp * ( &
            boxes(c_id)%fx(1, 1:nc:2, f_ixs) + &
            boxes(c_id)%fx(1, 2:nc:2, f_ixs))
    case (nb_hx)
       c_id = boxes(id)%children(ch_hxly)
       boxes(id)%fx(nc, 1:nch, f_ixs) = 0.5_dp * ( &
            boxes(c_id)%fx(nc, 1:nc:2, f_ixs) + &
            boxes(c_id)%fx(nc, 2:nc:2, f_ixs))
       c_id = boxes(id)%children(ch_hxhy)
       boxes(id)%fx(nc, nch+1:, f_ixs) = 0.5_dp * ( &
            boxes(c_id)%fx(nc, 1:nc:2, f_ixs) + &
            boxes(c_id)%fx(nc, 2:nc:2, f_ixs))
    case (nb_ly)
       c_id = boxes(id)%children(ch_lxly)
       boxes(id)%fx(1:nch, 1, f_ixs) = 0.5_dp * ( &
            boxes(c_id)%fx(1:nc:2, 1, f_ixs) + &
            boxes(c_id)%fx(2:nc:2, 1, f_ixs))
       c_id = boxes(id)%children(ch_hxly)
       boxes(id)%fx(nch+1:, 1, f_ixs) = 0.5_dp * ( &
            boxes(c_id)%fx(1:nc:2, 1, f_ixs) + &
            boxes(c_id)%fx(2:nc:2, 1, f_ixs))
    case (nb_hy)
       c_id = boxes(id)%children(ch_lxhy)
       boxes(id)%fx(1:nch, nc, f_ixs) = 0.5_dp * ( &
            boxes(c_id)%fx(1:nc:2, nc, f_ixs) + &
            boxes(c_id)%fx(2:nc:2, nc, f_ixs))
       c_id = boxes(id)%children(ch_hxhy)
       boxes(id)%fx(nch+1:, nc, f_ixs) = 0.5_dp * ( &
            boxes(c_id)%fx(1:nc:2, nc, f_ixs) + &
            boxes(c_id)%fx(2:nc:2, nc, f_ixs))
    end select
  end subroutine a2_flux_from_children

  subroutine a2_write_tree(tree, filename, cc_names, n_cycle, time)
    use m_vtk
    type(a2_t), intent(in) :: tree
    character(len=*)       :: filename, cc_names(:)
    integer, intent(in)    :: n_cycle
    real(dp), intent(in)   :: time
    integer                :: lvl, bc, bn, n, n_cells, n_nodes, n_grids
    integer                :: ig, i, j, id, n_ix, c_ix
    integer                :: cell_ix, node_ix
    integer                :: nodes_per_box, cells_per_box
    real(dp), allocatable  :: coords(:), cc_vars(:,:)
    integer, allocatable   :: offsets(:), connects(:), cell_types(:)
    type(vtk_t)            :: vtkf

    bc = tree%cfg%n_cell         ! number of Box Cells
    bn = tree%cfg%n_node         ! number of Box Nodes
    nodes_per_box = bn**2
    cells_per_box = bc**2

    n_grids = 0
    do lvl = 1, a2_n_levels(tree)
       do n = 1, size(tree%levels(lvl)%ids)
          id = tree%levels(lvl)%ids(n)
          if (.not. a2_has_children(tree%boxes(id))) n_grids = n_grids + 1
       end do
    end do
    n_nodes = nodes_per_box * n_grids
    n_cells = cells_per_box * n_grids

    allocate(coords(2 * n_nodes))
    allocate(cc_vars(n_cells, tree%cfg%n_cc))
    allocate(offsets(cells_per_box * n_grids))
    allocate(cell_types(cells_per_box * n_grids))
    allocate(connects(4 * cells_per_box * n_grids))

    cell_types = 8              ! VTK pixel type

    ig = 0
    do lvl = 1, a2_n_levels(tree)
       do n = 1, size(tree%levels(lvl)%ids)
          id = tree%levels(lvl)%ids(n)
          if (a2_has_children(tree%boxes(id))) cycle

          ig = ig + 1
          cell_ix = (ig-1) * cells_per_box
          node_ix = (ig-1) * nodes_per_box

          do j = 1, bn
             do i = 1, bn
                n_ix = 2 * (node_ix + (j-1) * bn + i) - 1
                coords(n_ix:n_ix+1) = a2_r_min(tree%boxes(id)) + &
                     [i-1,j-1] * a2_dr(tree%boxes(id))
             end do
          end do

          do j = 1, bc
             do i = 1, bc
                ! In vtk, indexing starts at 0, so subtract 1
                n_ix                      = node_ix + (j-1) * bn + i - 1
                c_ix                      = cell_ix + (j-1) * bc + i
                cc_vars(c_ix, :)          = tree%boxes(id)%cc(i, j, :)
                offsets(c_ix)             = 4 * c_ix
                connects(4*c_ix-3:4*c_ix) = [n_ix, n_ix+1, n_ix+bn, n_ix+bn+1]
             end do
          end do
       end do
    end do

    call vtk_ini_xml(vtkf, trim(filename), 'UnstructuredGrid')
    call vtk_dat_xml(vtkf, "UnstructuredGrid", .true.)
    call vtk_geo_xml(vtkf, coords, n_nodes, n_cells, 2, n_cycle, time)
    call vtk_con_xml(vtkf, connects, offsets, cell_types, n_cells)
    do n = 1, tree%cfg%n_cc
       call vtk_dat_xml(vtkf, "CellData", .true.)
       call vtk_var_r8_xml(vtkf, trim(cc_names(n)), cc_vars(:, n), n_cells)
       call vtk_dat_xml(vtkf, "CellData", .false.)
    end do
    call vtk_geo_xml_close(vtkf)
    call vtk_dat_xml(vtkf, "UnstructuredGrid", .false.)
    call vtk_end_xml(vtkf)
    print *, "Written ", trim(filename), ", n_grids", n_grids
  end subroutine a2_write_tree

  ! subroutine a2_write_tree(tree, filename, cc_names, cc_units, n_cycle, time)
  !   use m_write_silo
  !   type(a2_t), intent(in)          :: tree
  !   character(len=*)                :: filename, cc_names(:), cc_units(:)
  !   integer, intent(in)             :: n_cycle
  !   real(dp), intent(in)            :: time
  !   character(len=*), parameter     :: grid_name = "gg", var_name  = "vv"
  !   character(len=*), parameter     :: amr_name  = "amr"
  !   character(len=100), allocatable :: grid_list(:), var_list(:, :)
  !   integer                         :: lvl, i, id, ig, iv, bs, n_grids, dbix

  !   bs = tree%box_cells
  !   n_grids = 0
  !   do lvl = 1, a2_n_levels(tree)
  !      n_grids = n_grids + size(tree%levels(lvl)%ids)
  !   end do

  !   allocate(grid_list(n_grids))
  !   allocate(var_list(tree%n_cc, n_grids))

  !   call SILO_create_file(filename, dbix)
  !   ig = 0

  !   do lvl = 1, a2_n_levels(tree)
  !      do i = 1, size(tree%levels(lvl)%ids)
  !         id = tree%levels(lvl)%ids(i)
  !         ig = ig + 1
  !         write(grid_list(ig), "(A,I0)") grid_name, ig
  !         call SILO_add_grid(dbix, grid_list(ig), 2, &
  !              [bs+1, bs+1], tree%boxes(id)%r_min, tree%boxes(id)%dr)
  !         print *, id, tree%boxes(id)%r_min, tree%boxes(id)%dr
  !         do iv = 1, tree%n_cc
  !            write(var_list(iv, ig), "(A,I0)") trim(cc_names(iv)) // "_", ig
  !            call SILO_add_var(dbix, var_list(iv, ig), grid_list(ig), &
  !                 pack(tree%boxes(id)%cc(1:bs, 1:bs, iv), .true.), [bs, bs], &
  !                 trim(cc_units(iv)))
  !         end do
  !      end do
  !   end do

  !   call SILO_set_mmesh_grid(dbix, amr_name, grid_list, n_cycle, time)
  !   do iv = 1, tree%n_cc
  !      call SILO_set_mmesh_var(dbix, trim(cc_names(iv)), amr_name, &
  !           var_list(iv, :), n_cycle, time)
  !   end do
  !   call SILO_close_file(dbix)
  !   print *, "Number of grids", ig
  ! end subroutine a2_write_tree

end module m_afivo