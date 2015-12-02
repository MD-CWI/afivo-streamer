!> AFiVO code for $D-dimensional simulations
!> \author Jannis Teunissen
!> \copyright GPLv3

! The following replacements take place on this code:
! 1. $D -> 2 or 3 (dimension of code)
! 2. preprocess file with cpp
! 3. cat -s (merge multiple blank lines)

module m_a$D_core
  use m_a$D_t

  implicit none
  public

  private :: set_leaves_parents
  private :: set_nbs_$Dd
  private :: find_nb_$Dd
  private :: get_free_ids
  private :: set_ref_info
  private :: consistent_ref_flags
  private :: remove_children
  private :: add_children
  private :: set_child_ids
  private :: flux_from_children

contains

  !> Get tree info
  subroutine a$D_print_info(tree)
    type(a$D_t), intent(in)        :: tree       !< The tree

    if (.not. allocated(tree%lvls)) then
       print *, "a$D_init has not been called for this tree"
    else if (.not. tree%ready) then
       print *, "a$D_set_base has not been called for this tree"
    else
       write(*, "(A,I0)") " maximum allowed level:  ", tree%lvls_max
       write(*, "(A,I0)") " current maximum level:  ", tree%max_lvl
       write(*, "(A,I0)") " max index in box list:  ", tree%max_id
       write(*, "(A,I0)") " Boxes storage size:     ", size(tree%boxes)
       write(*, "(A,I0)") " Boxes used:             ", &
            count(tree%boxes(1:tree%max_id)%in_use)
       write(*, "(A,I0)") " box size (cells):       ", tree%n_cell
       write(*, "(A,I0)") " number of cc variables: ", tree%n_var_cell
       write(*, "(A,I0)") " number of fc variables: ", tree%n_var_face
       write(*, "(A,I0)") " Type of coordinates:    ", tree%coord_t
       write(*, "(A,2E12.4)") " min. coords:        ", tree%r_base
       write(*, "(A,2E12.4)") " dx at lvl 1:        ", tree%dr_base
    end if
  end subroutine a$D_print_info

  !> Initialize a $Dd tree type.
  subroutine a$D_init(tree, n_cell, n_var_cell, n_var_face, &
       dr, r_min, lvls_max, n_boxes, coarsen_to, coord)
    type(a$D_t), intent(out)        :: tree       !< The tree to initialize
    integer, intent(in)            :: n_cell     !< Boxes have n_cell^dim cells
    integer, intent(in)            :: n_var_cell !< Number of cell-centered variables
    integer, intent(in)            :: n_var_face !< Number of face-centered variables
    real(dp), intent(in)           :: dr         !< spacing of a cell at lvl 1
    !> Lowest coordinate of box at 1,1. Default is (0, 0)
    real(dp), intent(in), optional :: r_min($D)
    !> Create additional coarse grids down to this size. Default is -1 (which
    !> means don't do this)
    integer, intent(in), optional  :: coarsen_to
    !> Maximum number of levels. Default is 30
    integer, intent(in), optional  :: lvls_max
    !> Allocate initial storage for n_boxes. Default is 100
    integer, intent(in), optional  :: n_boxes
    integer, intent(in), optional  :: coord

    integer                        :: lvls_max_a, n_boxes_a, coarsen_to_a
    real(dp)                       :: r_min_a($D)
    integer                        :: lvl, min_lvl, coord_a

    ! Set default arguments if not present
    lvls_max_a = 30;   if (present(lvls_max)) lvls_max_a = lvls_max
    n_boxes_a = 100;   if (present(n_boxes)) n_boxes_a = n_boxes
    coarsen_to_a = -1; if (present(coarsen_to)) coarsen_to_a = coarsen_to
    r_min_a = 0.0_dp;  if (present(r_min)) r_min_a = r_min
    coord_a = a5_xyz;  if (present(coord)) coord_a = coord

    if (n_cell < 2)       stop "a$D_init: n_cell should be >= 2"
    if (btest(n_cell, 0)) stop "a$D_init: n_cell should be even"
    if (n_var_cell <= 0)  stop "a$D_init: n_var_cell should be > 0"
    if (n_boxes_a <= 0)   stop "a$D_init: n_boxes should be > 0"
    if (lvls_max_a <= 0)  stop "a$D_init: lvls_max should be > 0"
#if $D == 3
    if (coord_a == a5_cyl) stop "a$D_init: cannot have 3d cyl coords"
#endif

    allocate(tree%boxes(n_boxes_a))

    if (coarsen_to_a > 0) then
       ! Determine number of lvls for subtree
       ! TODO: remove subtree in future
       min_lvl = 1 - nint(log(real(n_cell, dp)/coarsen_to_a)/log(2.0_dp))

       if (2**(1-min_lvl) * coarsen_to_a /= n_cell) &
            stop "a$D_set_base: cannot coarsen to given value"
    else
       min_lvl = 1
    end if

    ! up to lvls_max_a+1 to add dummies that are always of size zero
    allocate(tree%lvls(min_lvl:lvls_max_a+1))

    do lvl = min_lvl, lvls_max_a+1
       allocate(tree%lvls(lvl)%ids(0))
       allocate(tree%lvls(lvl)%leaves(0))
       allocate(tree%lvls(lvl)%parents(0))
    end do

    tree%n_cell          = n_cell
    tree%n_var_cell      = n_var_cell
    tree%n_var_face      = n_var_face
    tree%r_base          = r_min_a
    tree%dr_base         = dr
    tree%lvls_max        = lvls_max_a
    tree%max_id          = 0
    tree%max_lvl         = 0
    tree%coord_t         = coord_a
  end subroutine a$D_init

  !> "Destroy" the data in a tree. Since we don't use pointers, you can also
  !> just let a tree get out of scope
  subroutine a$D_destroy(tree)
    type(a$D_t), intent(inout) :: tree
    integer                   :: lvl

    if (.not. tree%ready) stop "a$D_destroy: Tree was not fully initialized"

    deallocate(tree%boxes)
    do lvl = lbound(tree%lvls, 1), tree%lvls_max
       deallocate(tree%lvls(lvl)%ids)
       deallocate(tree%lvls(lvl)%leaves)
       deallocate(tree%lvls(lvl)%parents)
    end do
    tree%max_id = 0
  end subroutine a$D_destroy

  !> Create the base level of the tree, ix_list(:, id) stores the spatial index
  !> of box(id), nb_list(:, id) stores the neighbors of box(id)
  subroutine a$D_set_base(tree, ix_list, nb_list)
    type(a$D_t), intent(inout) :: tree !< Tree for which we set the base
    integer, intent(in)       :: ix_list(:, :) !< List of spatial indices for the initial boxes
    integer, intent(inout)    :: nb_list(:, :) !< Neighbors for the initial boxes
    integer                   :: n_boxes, i, id, nb, nb_id
    integer                   :: ix($D), lvl, offset

    if (any(ix_list < 1)) stop "a$D_set_base: need all ix_list > 0"
    if (tree%max_id > 0)  stop "a$D_set_base: this tree already has boxes"
    if (.not. allocated(tree%lvls)) stop "a$D_set_base: tree not initialized"

    ! Non-periodic and non-boundary condition neighbors only have to be
    ! specified from one side
    n_boxes = size(ix_list, 2)
    do i = 1, n_boxes
       do nb = 1, a$D_num_neighbors
          nb_id = nb_list(nb, i)
          if (nb_id > a5_no_box .and. nb_id /= i) &
               nb_list(a$D_nb_rev(nb), nb_id) = i
       end do
    end do

    if (any(nb_list == a5_no_box)) stop "a$D_set_base: unresolved neighbors"

    ! Check if we have enough space, if not, increase space
    if (n_boxes > size(tree%boxes(:))) then
       call a$D_resize_box_storage(tree, n_boxes)
    end if

    ! Create coarser levels which are copies of lvl 1
    do lvl = lbound(tree%lvls, 1), 1
       deallocate(tree%lvls(lvl)%ids)
       allocate(tree%lvls(lvl)%ids(n_boxes))

       call get_free_ids(tree, tree%lvls(lvl)%ids)
       offset = tree%lvls(lvl)%ids(1) - 1

       do i = 1, n_boxes
          id                         = tree%lvls(lvl)%ids(i)
          ix                         = ix_list(:, i)
          tree%boxes(id)%lvl         = lvl
          tree%boxes(id)%ix          = ix
          tree%boxes(id)%dr          = tree%dr_base * 0.5_dp**(lvl-1)
          tree%boxes(id)%r_min       = tree%r_base + &
               (ix - 1) * tree%dr_base * tree%n_cell
          tree%boxes(id)%n_cell      = tree%n_cell / (2**(1-lvl))
          tree%boxes(id)%coord_t     = tree%coord_t

          tree%boxes(id)%parent      = a5_no_box
          tree%boxes(id)%children(:) = a5_no_box ! Gets overwritten, see below

          ! Connectivity is the same for all lvls
          where (nb_list(:, i) > a5_no_box)
             tree%boxes(id)%neighbors = nb_list(:, i) + offset
          elsewhere
             tree%boxes(id)%neighbors = nb_list(:, i)
          end where

          call init_box(tree%boxes(id), tree%boxes(id)%n_cell, &
               tree%n_var_cell, tree%n_var_face)
       end do

       if (lvl == 1) then
          deallocate(tree%lvls(lvl)%leaves)
          allocate(tree%lvls(lvl)%leaves(n_boxes))
          tree%lvls(lvl)%leaves = tree%lvls(lvl)%ids
       else
          deallocate(tree%lvls(lvl)%parents)
          allocate(tree%lvls(lvl)%parents(n_boxes))
          tree%lvls(lvl)%parents = tree%lvls(lvl)%ids
       end if

       if (lvl > lbound(tree%lvls, 1)) then
          tree%boxes(tree%lvls(lvl-1)%ids)%children(1) = &
               tree%lvls(lvl)%ids
          tree%boxes(tree%lvls(lvl)%ids)%parent = &
               tree%lvls(lvl-1)%ids
       end if
    end do

    tree%max_lvl = 1
    tree%ready = .true.

  end subroutine a$D_set_base

  !> Reorder and resize the list of boxes. If the argument reorder is true,
  !> reorder the boxes but do not resize.
  subroutine a$D_tidy_up(tree, max_frac_used, goal_frac_used, &
       n_clean_min, reorder)
    use m_morton
    type(a$D_t), intent(inout)      :: tree          !< Tree for the list of boxes is reordered/resized
    real(dp), intent(in)           :: max_frac_used !< Maximum fraction of box-memory used
    real(dp), intent(in)           :: goal_frac_used !< If resizing, what fraction should be in use?
    integer, intent(in)            :: n_clean_min    !< Free up memory if at least this many boxes can be cleaned up
    logical, intent(in)            :: reorder   !< Do not resize the box list; only reorder it
    real(dp)                       :: frac_in_use
    integer                        :: n, lvl, id, old_size, new_size, n_clean
    integer                        :: max_id, n_used, n_stored, n_used_lvl
    integer, allocatable           :: ixs_sort(:), ixs_map(:)
    type(box$D_t), allocatable      :: boxes_cpy(:)
    integer(morton_k), allocatable :: mortons(:)

    if (.not. tree%ready) stop "Tree not ready"
    if (goal_frac_used > max_frac_used) &
         stop "a$D_tidy_up: need goal_frac_used < max_frac_used"
    if (max_frac_used > 1.0_dp) stop "a$D_tidy_up: need max_frac_used < 1"
    if (n_clean_min < 1)        stop "a$D_tidy_up: need n_clean_min > 0"

    max_id      = tree%max_id
    n_used      = count(tree%boxes(1:max_id)%in_use)
    old_size    = size(tree%boxes)
    frac_in_use = n_used / real(old_size, dp)
    n_clean     = nint((goal_frac_used - frac_in_use) * old_size)
    new_size    = old_size

    if (.not. reorder) then
       if (max_id > old_size * max_frac_used .or. &
            (frac_in_use < goal_frac_used .and. &
            n_clean > n_clean_min)) then
          new_size = max(1, nint(n_used/goal_frac_used))
       end if
    end if

    if (new_size /= old_size .or. reorder) then
       print *, "a$D_tidy_up: new size = ", new_size

       if (reorder) then
          allocate(boxes_cpy(n_used))  ! Need just enough space
       else
          allocate(boxes_cpy(new_size))
       end if

       allocate(ixs_map(0:max_id))
       ixs_map(0)       = 0
       n_stored         = 0

       do lvl = lbound(tree%lvls, 1), tree%max_lvl
          n_used_lvl = size(tree%lvls(lvl)%ids)
          allocate(mortons(n_used_lvl))
          allocate(ixs_sort(n_used_lvl))

          do n = 1, n_used_lvl
             id = tree%lvls(lvl)%ids(n)
             ! Note the -1, since our indices start at 1
             mortons(n) = morton_from_ix$D(tree%boxes(id)%ix-1)
          end do

          call morton_rank(mortons, ixs_sort)
          tree%lvls(lvl)%ids = tree%lvls(lvl)%ids(ixs_sort)

          do n = 1, n_used_lvl
             id = tree%lvls(lvl)%ids(n)
             boxes_cpy(n_stored + n) = tree%boxes(id)
             ixs_map(tree%lvls(lvl)%ids(n)) = n_stored + n
          end do

          tree%lvls(lvl)%ids = [(n_stored+n, n=1,n_used_lvl)]
          call set_leaves_parents(boxes_cpy, tree%lvls(lvl))
          n_stored = n_stored + n_used_lvl
          deallocate(mortons)
          deallocate(ixs_sort)
       end do

       ! Update id's to new indices
       do n = 1, n_used
          boxes_cpy(n)%parent = ixs_map(boxes_cpy(n)%parent)
          boxes_cpy(n)%children = ixs_map(boxes_cpy(n)%children)
          where (boxes_cpy(n)%neighbors > a5_no_box)
             boxes_cpy(n)%neighbors = ixs_map(boxes_cpy(n)%neighbors)
          end where
       end do

       if (reorder) then
          tree%boxes(1:n_used) = boxes_cpy ! Copy ordered data
          do n = n_used+1, max_id
             if (tree%boxes(n)%in_use) then
                ! Remove moved data
                call clear_box(tree%boxes(n))
             end if
          end do
       else
          deallocate(tree%boxes)
          call move_alloc(boxes_cpy, tree%boxes)
       end if

       tree%max_id = n_used
    end if

  end subroutine a$D_tidy_up

  !> Create a list of leaves and a list of parents for a level
  subroutine set_leaves_parents(boxes, level)
    type(box$D_t), intent(in)   :: boxes(:) !< List of boxes
    type(lvl_t), intent(inout) :: level !< Level type which contains the indices of boxes
    integer                    :: i, id, i_leaf, i_parent
    integer                    :: n_parents, n_leaves

    n_parents = count(a$D_has_children(boxes(level%ids)))
    n_leaves = size(level%ids) - n_parents

    if (n_parents /= size(level%parents)) then
       deallocate(level%parents)
       allocate(level%parents(n_parents))
    end if

    if (n_leaves /= size(level%leaves)) then
       deallocate(level%leaves)
       allocate(level%leaves(n_leaves))
    end if

    i_leaf   = 0
    i_parent = 0
    do i = 1, size(level%ids)
       id = level%ids(i)
       if (a$D_has_children(boxes(id))) then
          i_parent                = i_parent + 1
          level%parents(i_parent) = id
       else
          i_leaf               = i_leaf + 1
          level%leaves(i_leaf) = id
       end if
    end do
  end subroutine set_leaves_parents

  !> Mark box as active and allocate data storage for a box, for its cell- and
  !> face-centered data
  subroutine init_box(box, n_cell, n_cc, n_fc)
    type(box$D_t), intent(inout) :: box !< Box for which we allocate memory
    integer, intent(in)         :: n_cell !< Number of cells per dimension in the box
    integer, intent(in)         :: n_cc   !< Number of cell-centered variables
    integer, intent(in)         :: n_fc   !< Number of face-centered variables

    box%in_use = .true.

#if $D == 2
    allocate(box%cc(0:n_cell+1, 0:n_cell+1, n_cc))
    allocate(box%fx(n_cell+1,   n_cell,     n_fc))
    allocate(box%fy(n_cell,     n_cell+1,   n_fc))
#elif $D == 3
    allocate(box%cc(0:n_cell+1, 0:n_cell+1, 0:n_cell+1, n_cc))
    allocate(box%fx(n_cell+1,   n_cell,     n_cell,     n_fc))
    allocate(box%fy(n_cell,     n_cell+1,   n_cell,     n_fc))
    allocate(box%fz(n_cell,     n_cell,     n_cell+1,   n_fc))
#endif
  end subroutine init_box

  !> Deallocate data storage for a box and mark inactive
  subroutine clear_box(box)
    type(box$D_t), intent(inout) :: box

    box%in_use = .false.

    deallocate(box%cc)
    deallocate(box%fx)
    deallocate(box%fy)
#if $D == 3
    deallocate(box%fz)
#endif
    if (associated(box%ud)) deallocate(box%ud)
  end subroutine clear_box

  ! Set the neighbors of id (using their parent)
  subroutine set_nbs_$Dd(boxes, id)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    integer                     :: nb, nb_id

    do nb = 1, a$D_num_neighbors
       if (boxes(id)%neighbors(nb) == a5_no_box) then
          nb_id = find_nb_$Dd(boxes, id, nb)
          if (nb_id > a5_no_box) then
             boxes(id)%neighbors(nb) = nb_id
             boxes(nb_id)%neighbors(a$D_nb_rev(nb)) = id
          end if
       end if
    end do
  end subroutine set_nbs_$Dd

  !> Compute the child index for a box with spatial index ix.
  integer function a$D_ix_to_cix(ix)
    integer, intent(in) :: ix($D) !< Spatial index of the box
    ! The index can range from 1 (all ix odd) and 2**$D (all ix even)
#if $D == 2
    a$D_ix_to_cix = 4 - 2 * iand(ix(2), 1) - iand(ix(1), 1)
#elif $D == 3
    a3_ix_to_cix = 8 - 4 * iand(ix(3), 1) - 2 * iand(ix(2), 1) - iand(ix(1), 1)
#endif
  end function a$D_ix_to_cix

  !> Get the id of neighbor nb of boxes(id), through its parent
  function find_nb_$Dd(boxes, id, nb) result(nb_id)
    type(box$D_t), intent(in) :: boxes(:) !< List with all the boxes
    integer, intent(in)      :: id       !< Box whose neighbor we are looking for
    integer, intent(in)      :: nb       !< Neighbor index
    integer                  :: nb_id, p_id, c_ix, d, old_pid

    p_id    = boxes(id)%parent
    old_pid = p_id
    c_ix    = a$D_ix_to_cix(boxes(id)%ix)
    d       = a$D_nb_dim(nb)

    ! Check if neighbor is in same direction as ix is (low/high). If so,
    ! use neighbor of parent
    if (a$D_ch_low(c_ix, d) .eqv. a$D_nb_low(nb)) &
         p_id = boxes(p_id)%neighbors(nb)

    ! The child ix of the neighbor is reversed in direction d
    nb_id = boxes(p_id)%children(a$D_ch_rev(c_ix, d))
  end function find_nb_$Dd

  !> Resize box storage to new_size
  subroutine a$D_resize_box_storage(tree, new_size)
    type(a$D_t), intent(inout) :: tree    !< Tree to resize
    integer, intent(in)       :: new_size !< New size for the array boxes(:)
    type(box$D_t), allocatable :: boxes_cpy(:)

    if (.not. tree%ready) stop "Tree not ready"

    ! Store boxes in larger array boxes_cpy
    allocate(boxes_cpy(new_size))
    boxes_cpy(1:tree%max_id) = tree%boxes(1:tree%max_id)

    ! Deallocate current storage
    deallocate(tree%boxes)

    ! Use new array
    call move_alloc(boxes_cpy, tree%boxes)
  end subroutine a$D_resize_box_storage

  !> Adjust the refinement of a tree using the user-supplied set_ref_flags. If the
  !> argument n_changes is present, it contains the number of boxes that were
  !> (de)refined.
  !>
  !> This routine sets the bit a5_bit_new_children for each box that is refined.
  !> On input, the tree should be balanced. On output, the tree is still
  !> balanced, and its refinement is updated (with at most one level per call).
  subroutine a$D_adjust_refinement(tree, set_ref_flags, ref_info)
    type(a$D_t), intent(inout)           :: tree          !< Tree
    procedure(a$D_subr_ref)              :: set_ref_flags !< Refinement function
    type(ref_info_t), intent(inout)     :: ref_info !< Information about refinement
    integer                             :: lvl, id, i, c_ids(a$D_num_children)
    integer                             :: max_id_prev, max_id_req
    integer, allocatable                :: ref_flags(:)

    if (.not. tree%ready) stop "Tree not ready"
    max_id_prev = tree%max_id
    allocate(ref_flags(max_id_prev))

    ! Set refinement values for all boxes
    call consistent_ref_flags(tree, ref_flags, set_ref_flags)

    ! Check whether there is enough free space, otherwise extend the list
    max_id_req = max_id_prev + a$D_num_children * count(ref_flags == a5_refine)
    if (max_id_req > size(tree%boxes)) then
       print *, "Resizing box storage for refinement", max_id_req
       call a$D_resize_box_storage(tree, max_id_req)
    end if

    do lvl = 1, tree%lvls_max-1
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)

          if (id > max_id_prev) then
             cycle              ! This is a newly added box
          else if (ref_flags(id) == a5_refine) then
             ! Add children. First need to get num_children free id's
             call get_free_ids(tree, c_ids)
             call add_children(tree%boxes, id, c_ids, &
                  tree%n_var_cell, tree%n_var_face)
          else if (ref_flags(id) == a5_derefine) then
             ! Remove children
             call remove_children(tree%boxes, id)
          end if
       end do

       ! Update leaves / parents
       call set_leaves_parents(tree%boxes, tree%lvls(lvl))

       ! Set next level ids to children of this level
       call set_child_ids(tree%lvls(lvl)%parents, &
            tree%lvls(lvl+1)%ids, tree%boxes)

       ! Update connectivity of new children (id > max_id_prev)
       do i = 1, size(tree%lvls(lvl+1)%ids)
          id = tree%lvls(lvl+1)%ids(i)
          if (id > max_id_prev) call set_nbs_$Dd(tree%boxes, id)
       end do

       if (size(tree%lvls(lvl+1)%ids) == 0) exit
    end do

    tree%max_lvl = lvl

    ! Update leaves and parents for the last level, because we might have
    ! removed a refinement lvl.
    call set_leaves_parents(tree%boxes, tree%lvls(tree%max_lvl+1))

    ! Set information about the refinement
    call set_ref_info(tree, ref_flags, ref_info)

  end subroutine a$D_adjust_refinement

  !> Set information about the refinement for all "normal" levels (>= 1)
  subroutine set_ref_info(tree, ref_flags, ref_info)
    type(a$D_t), intent(in)         :: tree
    integer, intent(in)             :: ref_flags(:)
    type(ref_info_t), intent(inout) :: ref_info
    integer                         :: id, lvl, n, n_ch
    integer, allocatable            :: ref_count(:), drf_count(:)

    n_ch           = a$D_num_children
    ref_info%n_add = n_ch * count(ref_flags == a5_refine)
    ref_info%n_rm  = n_ch * count(ref_flags == a5_derefine)

    ! Use max_lvl+1 here because this lvl might have been completely removed
    if (allocated(ref_info%lvls)) deallocate(ref_info%lvls)
    allocate(ref_info%lvls(tree%max_lvl+1))
    allocate(ref_count(tree%max_lvl+1))
    allocate(drf_count(tree%max_lvl+1))

    ! Find the number of (de)refined boxes per level
    ref_count = 0
    drf_count = 0

    do id = 1, size(ref_flags)
       lvl = tree%boxes(id)%lvl

       if (ref_flags(id) == a5_refine) then
          ref_count(lvl) = ref_count(lvl) + 1
       else if (ref_flags(id) == a5_derefine) then
          drf_count(lvl) = drf_count(lvl) + 1
       end if
    end do

    ! Allocate storage per level
    ! There can be no new children at level 1
    allocate(ref_info%lvls(1)%add(0))
    allocate(ref_info%lvls(1)%rm(0))

    do lvl = 2, tree%max_lvl+1
       n = ref_count(lvl-1) * n_ch
       allocate(ref_info%lvls(lvl)%add(n))
       n = drf_count(lvl-1) * n_ch
       allocate(ref_info%lvls(lvl)%rm(n))
    end do

    ! Set the added and removed id's per level, these are the children of the
    ! (de)refined boxes
    ref_count = 0
    drf_count = 0

    do id = 1, size(ref_flags)
       lvl = tree%boxes(id)%lvl

       if (ref_flags(id) == a5_refine) then
          ref_count(lvl) = ref_count(lvl) + 1
          n = n_ch * (ref_count(lvl)-1) + 1
          ref_info%lvls(lvl+1)%add(n:n+n_ch-1) = tree%boxes(id)%children
       else if (ref_flags(id) == a5_derefine) then
          drf_count(lvl) = drf_count(lvl) + 1
          n = n_ch * (drf_count(lvl)-1) + 1
          ref_info%lvls(lvl+1)%rm(n:n+n_ch-1) = tree%boxes(id)%children
       end if
    end do
  end subroutine set_ref_info

  !> Get free ids from the boxes(:) array to store new boxes in. These ids are
  !> always consecutive.
  subroutine get_free_ids(tree, ids)
    type(a$D_t), intent(inout) :: tree
    integer, intent(out)      :: ids(:) !< Array which will be filled with free box ids
    integer                   :: i, max_id_prev, n_ids

    n_ids = size(ids)
    !$omp critical (crit_free_ids)
    max_id_prev = tree%max_id
    tree%max_id = tree%max_id + n_ids
    !$omp end critical (crit_free_ids)

    ids = [(max_id_prev + i, i=1,n_ids)]
  end subroutine get_free_ids

  !> Given the refinement function, return consistent refinement flags, that
  !> ensure that the tree is still balanced. Furthermore, it cannot derefine the
  !> base level, and it cannot refine above tree%lvls_max. The argument
  !> ref_flags is changed: for boxes that will be refined it holds a5_refine,
  !> for boxes that will be derefined it holds a5_derefine
  subroutine consistent_ref_flags(tree, ref_flags, set_ref_flags)
    type(a$D_t), intent(inout) :: tree         !< Tree for which we set refinement flags
    integer, intent(inout)    :: ref_flags(:) !< List of refinement flags for all boxes(:)
    procedure(a$D_subr_ref)    :: set_ref_flags     !< User-supplied refinement function.
    integer                   :: lvl, i, id, c_ids(a$D_num_children)
    integer                   :: nb, p_id, nb_id, p_nb_id
    integer                   :: lvls_max
    integer, allocatable      :: my_ref_flags(:)

    lvls_max = tree%lvls_max
    ref_flags(:) = -HUGE(1)
    my_ref_flags = ref_flags

    ! Set refinement flags for all boxes using set_ref_flags. Each thread first sets
    ! the flags on its own copy

    !$omp parallel private(lvl, i, id) firstprivate(my_ref_flags) shared(ref_flags)
    do lvl = 1, tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call set_ref_flags(tree%boxes, id, my_ref_flags)
       end do
       !$omp end do
    end do

    ! Now set refinement flags to the maximum (refine > keep ref > derefine)

    !$omp critical
    do i = 1, size(ref_flags)
       if (my_ref_flags(i) > ref_flags(i)) &
            ref_flags(i) = my_ref_flags(i)
    end do
    !$omp end critical
    !$omp end parallel

    deallocate(my_ref_flags)

    ! Set flags with unknown values to default (keep refinement)
    where (ref_flags > a5_do_ref) ref_flags = a5_kp_ref
    where (ref_flags < a5_rm_ref) ref_flags = a5_kp_ref

    ! Cannot refine beyond max level
    do i = 1, size(tree%lvls(lvls_max)%ids)
       id = tree%lvls(lvls_max)%ids(i)
       if (ref_flags(id) == a5_do_ref) ref_flags(id) = a5_kp_ref
    end do

    ! Ensure 2-1 balance
    do lvl = tree%max_lvl, 1, -1
       do i = 1, size(tree%lvls(lvl)%leaves) ! We only check leaf tree%boxes
          id = tree%lvls(lvl)%leaves(i)

          if (ref_flags(id) > a5_kp_ref) then ! This means refine
             ref_flags(id) = a5_refine ! Mark for actual refinement

             ! Ensure we will have the necessary neighbors
             do nb = 1, a$D_num_neighbors
                nb_id = tree%boxes(id)%neighbors(nb)
                if (nb_id == a5_no_box) then
                   ! Mark the parent containing neighbor for refinement
                   p_id = tree%boxes(id)%parent
                   p_nb_id = tree%boxes(p_id)%neighbors(nb)
                   ref_flags(p_nb_id) = a5_refine ! Mark for actual refinement
                end if
             end do

          else if (ref_flags(id) == a5_rm_ref) then
             ! Ensure we do not remove a required neighbor
             do nb = 1, a$D_num_neighbors
                nb_id = tree%boxes(id)%neighbors(nb)
                if (nb_id > a5_no_box) then
                   if (a$D_has_children(tree%boxes(nb_id)) .or. &
                        ref_flags(nb_id) > a5_kp_ref) then
                      ref_flags(id) = a5_kp_ref
                      exit
                   end if
                end if
             end do
          end if

       end do
    end do

    ! Make the (de)refinement flags consistent for blocks with children. Also
    ! ensure that at most one level can be removed at a time.
    do lvl = tree%max_lvl-1, 1, -1
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)

          ! Can only remove children if they are all marked for
          ! derefinement, and the box itself not for refinement.
          c_ids = tree%boxes(id)%children
          if (all(ref_flags(c_ids) == a5_rm_ref) .and. &
               ref_flags(id) <= a5_kp_ref) then
             ref_flags(id) = a5_derefine
          else
             ref_flags(id) = a5_kp_ref
          end if
       end do
    end do

  end subroutine consistent_ref_flags

  !> Remove the children of box id
  subroutine remove_children(boxes, id)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id       !< Id of box whose children will be removed
    integer                     :: ic, c_id, nb_id, nb_rev, nb

    do ic = 1, a$D_num_children
       c_id               = boxes(id)%children(ic)

       ! Remove from neighbors
       do nb = 1, a$D_num_neighbors
          nb_id = boxes(c_id)%neighbors(nb)
          if (nb_id > a5_no_box) then
             nb_rev = a$D_nb_rev(nb)
             boxes(nb_id)%neighbors(nb_rev) = a5_no_box
          end if
       end do

       call clear_box(boxes(c_id))
    end do

    boxes(id)%children = a5_no_box
  end subroutine remove_children

  !> Add children to box id, using the indices in c_ids
  subroutine add_children(boxes, id, c_ids, n_cc, n_fc)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id       !< Id of box that gets children
    integer, intent(in)         :: c_ids(a$D_num_children) !< Free ids for the children
    integer, intent(in)         :: n_cc                   !< Number of cell-centered variables
    integer, intent(in)         :: n_fc                   !< Number of face-centered variables
    integer                     :: i, nb, ch_nb(2**($D-1)), c_id, c_ix_base($D)

    boxes(id)%children = c_ids
    c_ix_base          = 2 * boxes(id)%ix - 1

    do i = 1, a$D_num_children
       c_id                  = c_ids(i)
       boxes(c_id)%ix        = c_ix_base + a$D_ch_dix(:,i)
       boxes(c_id)%lvl       = boxes(id)%lvl+1
       boxes(c_id)%parent    = id
       boxes(c_id)%tag       = a5_init_tag
       boxes(c_id)%children  = a5_no_box
       boxes(c_id)%neighbors = a5_no_box
       boxes(c_id)%n_cell    = boxes(id)%n_cell
       boxes(c_id)%coord_t   = boxes(id)%coord_t
       boxes(c_id)%dr        = 0.5_dp * boxes(id)%dr
       boxes(c_id)%r_min     = boxes(id)%r_min + 0.5_dp * boxes(id)%dr * &
            a$D_ch_dix(:,i) * boxes(id)%n_cell

       call init_box(boxes(c_id), boxes(id)%n_cell, n_cc, n_fc)
    end do

    ! Set boundary conditions at children
    do nb = 1, a$D_num_neighbors
       if (boxes(id)%neighbors(nb) < a5_no_box) then
          ch_nb = c_ids(a$D_ch_adj_nb(:, nb)) ! Neighboring children
          boxes(ch_nb)%neighbors(nb) = boxes(id)%neighbors(nb)
       end if
    end do
  end subroutine add_children

  !> Create a list c_ids(:) of all the children of p_ids(:). This is used after
  !> a level has been refined.
  subroutine set_child_ids(p_ids, c_ids, boxes)
    integer, intent(in)                 :: p_ids(:) !< All the parents ids
    integer, allocatable, intent(inout) :: c_ids(:) !< Output: all the children's ids
    type(box$D_t), intent(in)            :: boxes(:) !< List of all the boxes
    integer                             :: i, i0, i1, n_children

    n_children = a$D_num_children * size(p_ids)
    if (n_children /= size(c_ids)) then
       deallocate(c_ids)
       allocate(c_ids(n_children))
    end if

    do i = 1, size(p_ids)
       i1 = i * a$D_num_children
       i0 = i1 - a$D_num_children + 1
       c_ids(i0:i1) = boxes(p_ids(i))%children
    end do
  end subroutine set_child_ids

  !> Test if a box has children
  elemental logical function a$D_has_children(box)
    type(box$D_t), intent(in) :: box
    a$D_has_children = (box%children(1) /= a5_no_box)
  end function a$D_has_children

  !> Restrict fluxes from children to parents on refinement boundaries.
  subroutine a$D_consistent_fluxes(tree, f_ixs)
    type(a$D_t), intent(inout)     :: tree         !< Tree to operate on
    integer, intent(in)           :: f_ixs(:)     !< Indices of the fluxes
    integer                       :: lvl, i, id, nb, nb_id

    if (.not. tree%ready) stop "Tree not ready"
    !$omp parallel private(lvl, i, id, nb, nb_id)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl-1
       !$omp do
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          do nb = 1, a$D_num_neighbors
             nb_id = tree%boxes(id)%neighbors(nb)

             ! If the neighbor exists and has no children, set flux
             if (nb_id > a5_no_box) then
                if (.not. a$D_has_children(tree%boxes(nb_id))) then
                   call flux_from_children(tree%boxes, id, nb, f_ixs)
                end if
             end if
          end do
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine a$D_consistent_fluxes

  !> The neighbor nb has no children and id does, so set flux on the neighbor
  !> from our children. This ensures flux consistency at refinement boundary.
  subroutine flux_from_children(boxes, id, nb, f_ixs)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all the boxes
    integer, intent(in)         :: id        !< Id of box for which we set fluxes
    integer, intent(in)         :: nb        !< Direction in which fluxes are set
    integer, intent(in)         :: f_ixs(:)  !< Indices of the fluxes
    integer                     :: nc, nch, c_id, i_ch, i, ic, d, ioff($D)
    integer                     :: n_chnb, nb_id, i_nb

    nc     = boxes(id)%n_cell
    nch    = ishft(nc, -1) ! nc/2
    d      = a$D_nb_dim(nb)
    n_chnb = 2**($D-1)
    nb_id  = boxes(id)%neighbors(nb)

    if (a$D_nb_low(nb)) then
       i = 1
       i_nb = nc+1
    else
       i = nc+1
       i_nb = 1
    end if

    select case (d)
#if $D == 2
    case (1)
       do ic = 1, n_chnb
          ! Get index of child adjacent to neighbor
          i_ch = a2_ch_adj_nb(ic, nb)
          c_id = boxes(id)%children(i_ch)
          ! Index offset of child w.r.t. parent
          ioff = nch*a2_ch_dix(:, i_ch)
          boxes(nb_id)%fx(i_nb, ioff(2)+1:ioff(2)+nch, f_ixs) = 0.5_dp * ( &
               boxes(c_id)%fx(i, 1:nc:2, f_ixs) + &
               boxes(c_id)%fx(i, 2:nc:2, f_ixs))
       end do
    case (2)
       do ic = 1, n_chnb
          i_ch = a2_ch_adj_nb(ic, nb)
          c_id = boxes(id)%children(i_ch)
          ioff = nch*a2_ch_dix(:, i_ch)
          boxes(nb_id)%fy(ioff(1)+1:ioff(1)+nch, i_nb, f_ixs) = 0.5_dp * ( &
               boxes(c_id)%fy(1:nc:2, i, f_ixs) + &
               boxes(c_id)%fy(2:nc:2, i, f_ixs))
       end do
#elif $D == 3
    case (1)
       do ic = 1, n_chnb
          i_ch = a3_ch_adj_nb(ic, nb)
          c_id = boxes(id)%children(i_ch)
          ioff = nch*a3_ch_dix(:, i_ch)
          boxes(nb_id)%fx(i_nb, ioff(2)+1:ioff(2)+nch, &
               ioff(3)+1:ioff(3)+nch, f_ixs) = 0.25_dp * ( &
               boxes(c_id)%fx(i, 1:nc:2, 1:nc:2, f_ixs) + &
               boxes(c_id)%fx(i, 2:nc:2, 1:nc:2, f_ixs) + &
               boxes(c_id)%fx(i, 1:nc:2, 2:nc:2, f_ixs) + &
               boxes(c_id)%fx(i, 2:nc:2, 2:nc:2, f_ixs))
       end do
    case (2)
       do ic = 1, n_chnb
          i_ch = a3_ch_adj_nb(ic, nb)
          c_id = boxes(id)%children(i_ch)
          ioff = nch*a3_ch_dix(:, i_ch)
          boxes(nb_id)%fy(ioff(1)+1:ioff(1)+nch, i_nb, &
               ioff(3)+1:ioff(3)+nch, f_ixs) = 0.25_dp * ( &
               boxes(c_id)%fy(1:nc:2, i, 1:nc:2, f_ixs) + &
               boxes(c_id)%fy(2:nc:2, i, 1:nc:2, f_ixs) + &
               boxes(c_id)%fy(1:nc:2, i, 2:nc:2, f_ixs) + &
               boxes(c_id)%fy(2:nc:2, i, 2:nc:2, f_ixs))
       end do
    case (3)
       do ic = 1, n_chnb
          i_ch = a3_ch_adj_nb(ic, nb)
          c_id = boxes(id)%children(i_ch)
          ioff = nch*a3_ch_dix(:, i_ch)
          boxes(nb_id)%fz(ioff(1)+1:ioff(1)+nch, &
               ioff(2)+1:ioff(2)+nch, i_nb, f_ixs) = 0.25_dp * ( &
               boxes(c_id)%fz(1:nc:2, 1:nc:2, i, f_ixs) + &
               boxes(c_id)%fz(2:nc:2, 1:nc:2, i, f_ixs) + &
               boxes(c_id)%fz(1:nc:2, 2:nc:2, i, f_ixs) + &
               boxes(c_id)%fz(2:nc:2, 2:nc:2, i, f_ixs))
       end do
#endif
    end select
  end subroutine flux_from_children

end module m_a$D_core
