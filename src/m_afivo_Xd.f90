!> AFiVO code for $D-dimensional simulations
!> \author Jannis Teunissen
!> \copyright GPLv3

! The following replacements take place on this code:
! 1. $D -> 2 or 3 (dimension of code)
! 2. preprocess file with cpp
! 3. cat -s (merge multiple blank lines)
module m_afivo_$Dd
  use m_afivo_constants

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

#if $D == 2
  ! Children (same location as **corners**)
  integer, parameter :: a2_num_children = 4
  integer, parameter :: a2_ch_lxly = 1
  integer, parameter :: a2_ch_hxly = 2
  integer, parameter :: a2_ch_lxhy = 3
  integer, parameter :: a2_ch_hxhy = 4

  ! Neighboring indices for each child
  integer, parameter :: a2_ch_nbs(2, 4) = reshape([1,3,2,3,1,4,2,4], [2,4])
  ! Index offset for each child
  integer, parameter :: a2_ch_dix(2, 4) = reshape([0,0,1,0,0,1,1,1], [2,4])
  ! Reverse child index in each direction
  integer, parameter :: a2_ch_rev(4, 2) = reshape([2,1,4,3,3,4,1,2], [4,2])
  ! Children adjacent to a neighbor
  integer, parameter :: a2_ch_adj_nb(2, 4) = reshape([1,3,2,4,1,2,3,4], [2,4])
  ! Which children have a low index per dimension
  logical, parameter :: a2_ch_low(4, 2) = reshape([.true., .false., .true., &
       .false., .true., .true., .false., .false.], [4,2])

  ! Neighbor topology information
  integer, parameter :: a2_num_neighbors = 4
  integer, parameter :: a2_nb_lx = 1
  integer, parameter :: a2_nb_hx = 2
  integer, parameter :: a2_nb_ly = 3
  integer, parameter :: a2_nb_hy = 4

  ! Index offsets of neighbors
  integer, parameter :: a2_nb_dix(2, 4) = reshape([-1,0,1,0,0,-1,0,1], [2,4])
  ! Which neighbors have a lower index
  logical, parameter :: a2_nb_low(4) = [.true., .false., .true., .false.]
  ! Reverse neighbors
  integer, parameter :: a2_nb_rev(4) = [2, 1, 4, 3]
  ! Direction (dimension) for a neighbor
  integer, parameter :: a2_nb_dim(4) = [1, 1, 2, 2]
#elif $D == 3
  ! Children (same location as **corners**)
  integer, parameter :: a3_num_children = 8
  integer, parameter :: a3_ch_lxlylz = 1
  integer, parameter :: a3_ch_hxlylz = 2
  integer, parameter :: a3_ch_lxhylz = 3
  integer, parameter :: a3_ch_hxhylz = 4
  integer, parameter :: a3_ch_lxlyhz = 5
  integer, parameter :: a3_ch_hxlyhz = 6
  integer, parameter :: a3_ch_lxhyhz = 7
  integer, parameter :: a3_ch_hxhyhz = 8

  ! Neighboring indices for each child
  integer, parameter :: a3_ch_nbs(3, 8) = reshape( &
       [1,3,5, 2,3,5, 1,4,5, 2,4,5, &
       1,3,6, 2,3,6, 1,4,6, 2,4,6], [3,8])
  ! Index offset for each child
  integer, parameter :: a3_ch_dix(3, 8) = reshape( &
       [0,0,0, 1,0,0, 0,1,0, 1,1,0, &
       0,0,1, 1,0,1, 0,1,1, 1,1,1], [3,8])
  ! Reverse child index in each direction
  integer, parameter :: a3_ch_rev(8, 3) = reshape( &
       [2,1,4,3,6,5,8,7, 3,4,1,2,7,8,5,6, 5,6,7,8,1,2,3,4], [8,3])
  ! Children adjacent to a neighbor
  integer, parameter :: a3_ch_adj_nb(4, 6) = reshape( &
       [1,3,5,7, 2,4,6,8, 1,2,5,6, 3,4,7,8, 1,2,3,4, 5,6,7,8], [4,6])
  ! Which children have a low index per dimension
  logical, parameter :: a3_ch_low(8, 3) = reshape([ &
       .true., .false., .true., .false., .true., .false., .true., .false., &
       .true., .true., .false., .false., .true., .true., .false., .false., &
       .true., .true., .true., .true., .false., .false., .false., .false.], &
       [8,3])

  ! Neighbor topology information
  integer, parameter :: a3_num_neighbors = 6
  integer, parameter :: a3_nb_lx = 1
  integer, parameter :: a3_nb_hx = 2
  integer, parameter :: a3_nb_ly = 3
  integer, parameter :: a3_nb_hy = 4
  integer, parameter :: a3_nb_lz = 5
  integer, parameter :: a3_nb_hz = 6
  ! Index offsets of neighbors
  integer, parameter :: a3_nb_dix(3, 6) = reshape( &
       [-1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1], [3,6])
  ! Which neighbors have a lower index
  logical, parameter :: a3_nb_low(6) = &
       [.true., .false., .true., .false., .true., .false.]
  ! Reverse neighbors
  integer, parameter :: a3_nb_rev(6) = [2, 1, 4, 3, 6, 5]
  ! Direction (dimension) for a neighbor
  integer, parameter :: a3_nb_dim(6) = [1, 1, 2, 2, 3, 3]
#endif

  !> The basic building block of afivo: a box with cell-centered and face
  !> centered data, and information about its position, neighbors, children etc.
  type box$D_t
     integer               :: lvl    !< level of the box
     integer               :: tag    !< for setting tag bits
     integer               :: ix($D) !< index in the domain
     integer               :: parent !< index of parent in box list
     !> index of children in box list
     integer               :: children(a$D_num_children)
     !> index of neighbors in box list
     integer               :: neighbors(a$D_num_neighbors)
     integer               :: n_cell    !< number of cells per dimension
     real(dp)              :: dr        !< width/height of a cell
     real(dp)              :: r_min($D) !< min coords. of box
     integer, allocatable  :: ud(:)     !< User data (can be anything)
#if $D   == 2
     real(dp), allocatable :: cc(:, :, :) !< cell centered variables
     real(dp), allocatable :: fx(:, :, :) !< x-face centered variables
     real(dp), allocatable :: fy(:, :, :) !< y-face centered variables
#elif $D == 3
     real(dp), allocatable :: cc(:, :, :, :) !< cell centered variables
     real(dp), allocatable :: fx(:, :, :, :) !< x-face centered variables
     real(dp), allocatable :: fy(:, :, :, :) !< y-face centered variables
     real(dp), allocatable :: fz(:, :, :, :) !< z-face centered variables
#endif
  end type box$D_t

  !> Type which contains the indices of all boxes at a refinement level, as well
  !> as a list with all the "leaf" boxes and non-leaf (parent) boxes
  type lvl_t
     integer, allocatable :: ids(:)     !< indices of boxes of level
     integer, allocatable :: leaves(:)  !< all ids(:) that are leaves
     integer, allocatable :: parents(:) !< all ids(:) that have children
  end type lvl_t

  !> Type which stores all the boxes and levels, as well as some information
  !> about the number of boxes, variables and levels.
  type a$D_t
     integer                    :: lvls_max   !< maximum allowed level
     integer                    :: max_lvl    !< current maximum level
     integer                    :: max_id     !< max index in box list
     integer                    :: n_cell     !< number of cells per dimension
     integer                    :: n_var_cell !< number of cc variables
     integer                    :: n_var_face !< number of fc variables
     real(dp)                   :: r_base($D) !< coords of box at index (1,1)
     real(dp)                   :: dr_base    !< cell spacing at lvl 1
     type(lvl_t), allocatable   :: lvls(:)    !< list storing the tree levels
     type(box$D_t), allocatable :: boxes(:)   !< list of all boxes
  end type a$D_t

  abstract interface
     !> Function that gets a box and returns an int
     integer function a$D_to_int_f(box)
       import
       type(box$D_t), intent(in) :: box
     end function a$D_to_int_f

     !> Subroutine that gets a box
     subroutine a$D_subr(box)
       import
       type(box$D_t), intent(inout) :: box
     end subroutine a$D_subr

     !> Subroutine that gets a box and an array of reals
     subroutine a$D_subr_arg(box, rarg)
       import
       type(box$D_t), intent(inout) :: box
       real(dp), intent(in)        :: rarg(:)
     end subroutine a$D_subr_arg

     !> Subroutine that gets a list of boxes and a box id
     subroutine a$D_subr_boxes(boxes, id)
       import
       type(box$D_t), intent(inout) :: boxes(:)
       integer, intent(in)         :: id
     end subroutine a$D_subr_boxes

     !> Subroutine that gets a list of boxes, an id and an array of reals
     subroutine a$D_subr_boxes_arg(boxes, id, rarg)
       import
       type(box$D_t), intent(inout) :: boxes(:)
       integer, intent(in)         :: id
       real(dp), intent(in)        :: rarg(:)
     end subroutine a$D_subr_boxes_arg

     !> Subroutine for filling ghost cells on the sides of boxes, near
     !> refinement or physical boundaries.
     subroutine a$D_subr_gc(boxes, id, nb, iv)
       import
       type(box$D_t), intent(inout) :: boxes(:) !< Array with all boxes
       integer, intent(in)         :: id       !< Id of the box that needs to have ghost cells filled
       integer, intent(in)         :: nb       !< Neighbor direction in which ghost cells need to be filled,
                                               !< see a$D_nb_lx etc.
       integer, intent(in)         :: iv       !< Variable for which ghost cells are filled
     end subroutine a$D_subr_gc
  end interface

  private :: set_leaves_parents
  private :: alloc_box
  private :: dealloc_box
  private :: set_nbs_$Dd
  private :: find_nb_$Dd
  private :: get_free_ids
  private :: set_ref_flags
  private :: remove_children
  private :: add_children
  private :: set_child_ids
  private :: l2i

contains

  !> Initialize a $Dd tree type.
  subroutine a$D_init(tree, n_cell, n_var_cell, n_var_face, &
       dr, r_min, lvls_max, n_boxes, coarsen_to)
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

    integer                        :: lvls_max_a, n_boxes_a, coarsen_to_a
    real(dp)                       :: r_min_a($D)
    integer                        :: lvl, min_lvl

    ! Set default arguments if not present
    lvls_max_a = 30;   if (present(lvls_max)) lvls_max_a = lvls_max
    n_boxes_a = 100;   if (present(n_boxes)) n_boxes_a = n_boxes
    coarsen_to_a = -1; if (present(coarsen_to)) coarsen_to_a = coarsen_to
    r_min_a = 0.0_dp;  if (present(r_min)) r_min_a = r_min

    if (n_cell < 2)       stop "a$D_init: n_cell should be >= 2"
    if (btest(n_cell, 0)) stop "a$D_init: n_cell should be even"
    if (n_var_cell <= 0)  stop "a$D_init: n_var_cell should be > 0"
    if (n_boxes_a <= 0)     stop "a$D_init: n_boxes should be > 0"
    if (lvls_max_a <= 0)     stop "a$D_init: lvls_max should be > 0"

    allocate(tree%boxes(n_boxes_a))

    if (coarsen_to_a > 0) then
       ! Determine number of lvls for subtree
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

    tree%boxes(:)%tag = 0
    tree%n_cell       = n_cell
    tree%n_var_cell   = n_var_cell
    tree%n_var_face   = n_var_face
    tree%r_base       = r_min_a
    tree%dr_base      = dr
    tree%lvls_max     = lvls_max_a
    tree%max_id       = 0
    tree%max_lvl      = 0
  end subroutine a$D_init

  !> "Destroy" the data in a tree. Since we don't use pointers, you can also
  !> just let a tree get out of scope
  subroutine a$D_destroy(tree)
    type(a$D_t), intent(inout) :: tree
    integer                   :: lvl

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

    ! Neighbors only have to be specified from one side, mirror them
    n_boxes = size(ix_list, 2)
    do i = 1, n_boxes
       do nb = 1, a$D_num_neighbors
          nb_id = nb_list(nb, i)
          if (nb_id > a5_no_box) nb_list(a$D_nb_rev(nb), nb_id) = i
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
          tree%boxes(id)%tag         = ibset(0, a5_bit_in_use)
          tree%boxes(id)%dr          = tree%dr_base * 0.5_dp**(lvl-1)
          tree%boxes(id)%r_min       = (ix - 1) * tree%dr_base * tree%n_cell
          tree%boxes(id)%n_cell      = tree%n_cell / (2**(1-lvl))

          tree%boxes(id)%parent      = a5_no_box
          tree%boxes(id)%children(:) = a5_no_box ! Gets overwritten, see below

          ! Connectivity is the same for all lvls
          where (nb_list(:, i) > a5_no_box)
             tree%boxes(id)%neighbors = nb_list(:, i) + offset
          elsewhere
             tree%boxes(id)%neighbors = nb_list(:, i)
          end where

          call alloc_box(tree%boxes(id), tree%boxes(id)%n_cell, &
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

  end subroutine a$D_set_base

  !> Call procedure for each box in tree
  subroutine a$D_loop_box(tree, my_procedure, leaves_only)
    type(a$D_t), intent(inout)     :: tree
    procedure(a$D_subr)           :: my_procedure
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                       :: lvl, i, id

    leaves = .false.; if (present(leaves_only)) leaves = leaves_only

    !$omp parallel private(lvl, i, id)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       if (leaves) then
          !$omp do
          do i = 1, size(tree%lvls(lvl)%leaves)
             id = tree%lvls(lvl)%leaves(i)
             call my_procedure(tree%boxes(id))
          end do
          !$omp end do
       else
          !$omp do
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             call my_procedure(tree%boxes(id))
          end do
          !$omp end do
       end if
    end do
    !$omp end parallel
  end subroutine a$D_loop_box

  !> Call procedure for each box in tree, with argument rarg
  subroutine a$D_loop_box_arg(tree, my_procedure, rarg, leaves_only)
    type(a$D_t), intent(inout)     :: tree
    procedure(a$D_subr_arg)        :: my_procedure
    real(dp), intent(in)          :: rarg(:)
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                       :: lvl, i, id

    leaves = .false.; if (present(leaves_only)) leaves = leaves_only

    !$omp parallel private(lvl, i, id)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       if (leaves) then
          !$omp do
          do i = 1, size(tree%lvls(lvl)%leaves)
             id = tree%lvls(lvl)%leaves(i)
             call my_procedure(tree%boxes(id), rarg)
          end do
          !$omp end do
       else
          !$omp do
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             call my_procedure(tree%boxes(id), rarg)
          end do
          !$omp end do
       end if
    end do
    !$omp end parallel
  end subroutine a$D_loop_box_arg

  !> Call procedure for each id in tree, giving the list of boxes
  subroutine a$D_loop_boxes(tree, my_procedure, leaves_only)
    type(a$D_t), intent(inout)     :: tree
    procedure(a$D_subr_boxes)      :: my_procedure
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                       :: lvl, i, id

    leaves = .false.; if (present(leaves_only)) leaves = leaves_only

    !$omp parallel private(lvl, i, id)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       if (leaves) then
          !$omp do
          do i = 1, size(tree%lvls(lvl)%leaves)
             id = tree%lvls(lvl)%leaves(i)
             call my_procedure(tree%boxes, id)
          end do
          !$omp end do
       else
          !$omp do
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             call my_procedure(tree%boxes, id)
          end do
          !$omp end do
       end if
    end do
    !$omp end parallel
  end subroutine a$D_loop_boxes

  !> Call procedure for each id in tree, giving the list of boxes
  subroutine a$D_loop_boxes_arg(tree, my_procedure, rarg)
    type(a$D_t), intent(inout)    :: tree
    procedure(a$D_subr_boxes_arg) :: my_procedure
    real(dp), intent(in)         :: rarg(:)
    integer                      :: lvl, i, id

    !$omp parallel private(lvl, i, id)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call my_procedure(tree%boxes, id, rarg)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine a$D_loop_boxes_arg

  !> Clear "bit" from all the tags in the tree
  subroutine a$D_clear_tagbit(tree, bit)
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)       :: bit
    tree%boxes(1:tree%max_id)%tag = ibclr(tree%boxes(1:tree%max_id)%tag, bit)
  end subroutine a$D_clear_tagbit

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

    if (goal_frac_used > max_frac_used) &
         stop "a$D_tidy_up: need goal_frac_used < max_frac_used"
    if (max_frac_used > 1.0_dp) stop "a$D_tidy_up: need max_frac_used < 1"
    if (n_clean_min < 1)        stop "a$D_tidy_up: need n_clean_min > 0"

    max_id      = tree%max_id
    n_used      = count(btest(tree%boxes(1:max_id)%tag, a5_bit_in_use))
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
             call dealloc_box(tree%boxes(n)) ! Remove unused data
             tree%boxes(n)%tag = 0
          end do
       else
          deallocate(tree%boxes)
          call move_alloc(boxes_cpy, tree%boxes)
          tree%boxes(n_used+1:)%tag = 0
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

  !> Allocate data storage for a box, for its cell- and face-centered data
  subroutine alloc_box(box, n_cell, n_cc, n_fc)
    type(box$D_t), intent(inout) :: box !< Box for which we allocate memory
    integer, intent(in)         :: n_cell !< Number of cells per dimension in the box
    integer, intent(in)         :: n_cc   !< Number of cell-centered variables
    integer, intent(in)         :: n_fc   !< Number of face-centered variables

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
  end subroutine alloc_box

  !> Deallocate data storage for a box
  subroutine dealloc_box(box)
    type(box$D_t), intent(inout) :: box
    deallocate(box%cc)
    deallocate(box%fx)
    deallocate(box%fy)
#if $D == 3
    deallocate(box%fz)
#endif
    if (allocated(box%ud)) deallocate(box%ud)
  end subroutine dealloc_box

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

  !> Get the offset of a box with respect to its parent (e.g. in 2d, there can
  !> be a child at offset 0,0, one at n_cell/2,0, one at 0,n_cell/2 etc.)
  function a$D_get_child_offset(box) result(ix_offset)
    type(box$D_t), intent(in) :: box !< A child box
    integer                  :: ix_offset($D)

    ! Where the box index is even, set offset to n_cell/2, elsewhere set to 0
    ix_offset =  iand(box%ix-1, 1) * ishft(box%n_cell, -1)
  end function a$D_get_child_offset

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

    ! Store boxes in larger array boxes_cpy
    allocate(boxes_cpy(new_size))
    boxes_cpy(1:tree%max_id)      = tree%boxes(1:tree%max_id)
    boxes_cpy(tree%max_id+1:)%tag = 0 ! empty tag

    ! Deallocate current storage
    deallocate(tree%boxes)

    ! Use new array
    call move_alloc(boxes_cpy, tree%boxes)
  end subroutine a$D_resize_box_storage

  !> Adjust the refinement of a tree using the user-supplied ref_func. If the
  !> argument n_changes is present, it contains the number of boxes that were
  !> (de)refined.
  !>
  !> This routine sets the bit a5_bit_new_children for each box that is refined.
  !> On input, the tree should be balanced. On output, the tree is still
  !> balanced, and its refinement is updated (with at most one level per call).
  subroutine a$D_adjust_refinement(tree, ref_func, n_changes)
    type(a$D_t), intent(inout)      :: tree !< Tree whose refinement will be adjusted
    procedure(a$D_to_int_f)         :: ref_func !< User supplied refinement function
    integer, intent(out), optional :: n_changes !< Number of (de)refined boxes
    integer                        :: lvl, id, i, c_ids(a$D_num_children), i_c
    integer                        :: max_id_prev, max_id_req
    integer, allocatable           :: ref_flags(:)

    max_id_prev = tree%max_id
    allocate(ref_flags(max_id_prev))

    ! Clear refinement tags from previous calls
    call a$D_clear_tagbit(tree, a5_bit_new_children)

    ! Set refinement values for all boxes
    call set_ref_flags(tree, ref_flags, ref_func)

    if (present(n_changes)) n_changes = count(ref_flags /= a5_kp_ref)

    ! Check whether there is enough free space, otherwise extend the list
    max_id_req = max_id_prev + a$D_num_children * count(ref_flags == a5_do_ref)
    if (max_id_req > size(tree%boxes)) then
       print *, "Resizing box storage for refinement", max_id_req
       call a$D_resize_box_storage(tree, max_id_req)
    end if

    do lvl = 1, tree%lvls_max-1
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)

          if (id > max_id_prev) then
             cycle              ! This is a newly added box
          else if (ref_flags(id) == a5_do_ref) then
             ! Add children. First need to get num_children free id's
             call get_free_ids(tree, c_ids)
             call add_children(tree%boxes, id, c_ids, &
                  tree%n_var_cell, tree%n_var_face)
          else if (ref_flags(id) == a5_rm_children) then
             ! Remove children
             call remove_children(tree%boxes, id)
          end if
       end do

       ! Set next level ids to children of this level
       call set_child_ids(tree%lvls(lvl)%ids, &
            tree%lvls(lvl+1)%ids, tree%boxes)

       ! Update connectivity of children
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          if (a$D_has_children(tree%boxes(id))) then
             do i_c = 1, a$D_num_children
                call set_nbs_$Dd(tree%boxes, tree%boxes(id)%children(i_c))
             end do
          end if
       end do

       if (size(tree%lvls(lvl+1)%ids) == 0) exit
    end do

    tree%max_lvl = lvl

    ! Update leaves and parents. max_lvl+1 because we might have removed a
    ! refinement lvl.
    do lvl = 1, tree%max_lvl+1
       call set_leaves_parents(tree%boxes, tree%lvls(lvl))
    end do
  end subroutine a$D_adjust_refinement

  !> Get free ids from the boxes(:) array to store new boxes in. These ids are
  !> always consecutive.
  subroutine get_free_ids(tree, ids)
    type(a$D_t), intent(inout) :: tree
    integer, intent(out)      :: ids(:) !< Array which will be filled with free box ids
    integer                   :: i, max_id_prev, n_ids

    n_ids = size(ids)
    !$omp critical
    max_id_prev = tree%max_id
    tree%max_id = tree%max_id + n_ids
    !$omp end critical

    ids = [(max_id_prev + i, i=1,n_ids)]
  end subroutine get_free_ids

  !> Given the refinement function, return consistent refinement flags, that
  !> ensure that the tree is still balanced. Furthermore, it cannot derefine the
  !> base level, and it cannot refine above tree%lvls_max.
  subroutine set_ref_flags(tree, ref_flags, ref_func)
    type(a$D_t), intent(inout) :: tree         !< Tree for which we set refinement flags
    integer, intent(inout)    :: ref_flags(:) !< List of refinement flags for all boxes(:)
    procedure(a$D_to_int_f)    :: ref_func     !< User-supplied refinement function.
    integer                   :: lvl, i, id, c_ids(a$D_num_children)
    integer                   :: nb, p_id, nb_id, p_nb_id
    integer                   :: lvls_max

    lvls_max = tree%lvls_max
    ref_flags(:) = a5_kp_ref      ! Used indices are overwritten

    ! Set refinement flags for all boxes using ref_func
    do lvl = 1, tree%max_lvl
       do i = 1, size(tree%lvls(lvl)%ids)
          id            = tree%lvls(lvl)%ids(i)
          ref_flags(id) = ref_func(tree%boxes(id))
       end do
    end do

    ! Cannot derefine lvl 1
    do i = 1, size(tree%lvls(1)%ids)
       id = tree%lvls(1)%ids(i)
       if (ref_flags(id) == a5_rm_ref) ref_flags(id) = a5_kp_ref
    end do

    ! Cannot refine beyond max level
    do i = 1, size(tree%lvls(lvls_max)%ids)
       id = tree%lvls(lvls_max)%ids(i)
       if (ref_flags(id) == a5_do_ref) ref_flags(id) = a5_kp_ref
    end do

    ! Ensure 2-1 balance.
    do lvl = tree%max_lvl, 2, -1
       do i = 1, size(tree%lvls(lvl)%leaves) ! We only check leaf tree%boxes
          id = tree%lvls(lvl)%leaves(i)

          if (ref_flags(id) == a5_do_ref) then
             ! Ensure we will have the necessary neighbors
             do nb = 1, a$D_num_neighbors
                nb_id = tree%boxes(id)%neighbors(nb)
                if (nb_id == a5_no_box) then
                   ! Mark the parent containing neighbor for refinement
                   p_id = tree%boxes(id)%parent
                   p_nb_id = tree%boxes(p_id)%neighbors(nb)
                   ref_flags(p_nb_id) = a5_do_ref
                end if
             end do
          else if (ref_flags(id) == a5_rm_ref) then
             ! Ensure we do not remove a required neighbor
             do nb = 1, a$D_num_neighbors
                nb_id = tree%boxes(id)%neighbors(nb)
                if (nb_id > a5_no_box) then
                   if (a$D_has_children(tree%boxes(nb_id)) .or. &
                        ref_flags(nb_id) == a5_do_ref) then
                      ref_flags(id) = a5_kp_ref
                   end if
                end if
             end do
          end if
       end do
    end do

    ! Make the (de)refinement flags consistent for blocks with children.
    do lvl = tree%max_lvl-1, 1, -1
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)

          ! Can only remove children if they are all marked for
          ! derefinement, and the box itself not for refinement.
          c_ids = tree%boxes(id)%children
          if (all(ref_flags(c_ids) == a5_rm_ref) .and. &
               ref_flags(id) /= a5_do_ref) then
             ref_flags(id) = a5_rm_children
          else
             ref_flags(id) = a5_kp_ref
             where (ref_flags(c_ids) == a5_rm_ref)
                ref_flags(c_ids) = a5_kp_ref
             end where
          end if
       end do
    end do

  end subroutine set_ref_flags

  !> Remove the children of box id
  subroutine remove_children(boxes, id)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id       !< Id of box whose children will be removed
    integer                     :: ic, c_id, nb_id, nb_rev, nb

    do ic = 1, a$D_num_children
       c_id            = boxes(id)%children(ic)
       boxes(c_id)%tag = 0      ! clear tag

       do nb = 1, a$D_num_neighbors             ! Remove from neighbors
          nb_id = boxes(c_id)%neighbors(nb)
          if (nb_id > a5_no_box) then
             nb_rev = a$D_nb_rev(nb)
             boxes(nb_id)%neighbors(nb_rev) = a5_no_box
          end if
       end do
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
    boxes(id)%tag      = ibset(boxes(id)%tag, a5_bit_new_children)

    do i = 1, a$D_num_children
       c_id                  = c_ids(i)
       boxes(c_id)%ix        = c_ix_base + a$D_ch_dix(:,i)
       boxes(c_id)%lvl       = boxes(id)%lvl+1
       boxes(c_id)%parent    = id
       boxes(c_id)%children  = a5_no_box
       boxes(c_id)%neighbors = a5_no_box
       boxes(c_id)%n_cell    = boxes(id)%n_cell
       boxes(c_id)%dr        = 0.5_dp * boxes(id)%dr
       boxes(c_id)%r_min     = boxes(id)%r_min + 0.5_dp * boxes(id)%dr * &
            a$D_ch_dix(:,i) * boxes(id)%n_cell
       boxes(c_id)%tag       = ibset(0, a5_bit_in_use)

       call alloc_box(boxes(c_id), boxes(id)%n_cell, n_cc, n_fc)
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
    integer                             :: i, ip, i_c, n_children

    ! Count a$D_num_children times the number of refined parent blocks
    n_children = a$D_num_children * count(a$D_has_children(boxes(p_ids)))
    if (n_children /= size(c_ids)) then
       deallocate(c_ids)
       allocate(c_ids(n_children))
    end if

    i_c = 0
    do i = 1, size(p_ids)
       ip = p_ids(i)
       if (a$D_has_children(boxes(ip))) then
          c_ids(i_c+1:i_c+a$D_num_children) = boxes(ip)%children
          i_c = i_c + a$D_num_children
       end if
    end do
  end subroutine set_child_ids

  !> Test if a box has children
  elemental logical function a$D_has_children(box)
    type(box$D_t), intent(in) :: box
    a$D_has_children = (box%children(1) /= a5_no_box)
  end function a$D_has_children

  !> Return dr at lvl
  pure function a$D_lvl_dr(tree, lvl) result(dr)
    type(a$D_t), intent(in) :: tree
    integer, intent(in)    :: lvl
    real(dp)               :: dr !< Output: dr at the finest lvl of the tree
    dr = tree%dr_base * 0.5_dp**(lvl-1)
  end function a$D_lvl_dr

  !> Return finest dr that is used in the tree
  pure function a$D_min_dr(tree) result(dr)
    type(a$D_t), intent(in) :: tree
    real(dp)               :: dr !< Output: dr at the finest lvl of the tree
    dr = a$D_lvl_dr(tree, tree%max_lvl)
  end function a$D_min_dr

  !> Returns whether r is inside or within a distance d from box
  pure function a$D_r_inside(box, r, d) result(inside)
    type(box$D_t), intent(in)       :: box
    real(dp), intent(in)           :: r($D)
    real(dp), intent(in), optional :: d
    real(dp)                       :: r_max($D)
    logical                        :: inside

    r_max = box%r_min + box%dr * box%n_cell
    if (present(d)) then
       inside = all(r+d >= box%r_min) .and. all(r-d <= r_max)
    else
       inside = all(r >= box%r_min) .and. all(r <= r_max)
    end if
  end function a$D_r_inside

  !> Return the coordinate of the center of a box
  pure function a$D_r_center(box) result(r_center)
    type(box$D_t), intent(in) :: box
    real(dp)                 :: r_center($D)
    r_center = box%r_min + 0.5_dp * box%n_cell * box%dr
  end function a$D_r_center

  !> Get the location of the cell center with index cc_ix
  pure function a$D_r_cc(box, cc_ix) result(r)
    type(box$D_t), intent(in) :: box
    integer, intent(in)      :: cc_ix($D)
    real(dp)                 :: r($D)
    r = box%r_min + (cc_ix-0.5_dp) * box%dr
  end function a$D_r_cc

  !> Get a general location with real index cc_ix (like a$D_r_cc).
  pure function a$D_rr_cc(box, cc_ix) result(r)
    type(box$D_t), intent(in) :: box
    real(dp), intent(in)     :: cc_ix($D)
    real(dp)                 :: r($D)
    r = box%r_min + (cc_ix-0.5_dp) * box%dr
  end function a$D_rr_cc

  !> Get the location of a node (cell-corner) with index nd_ix
  pure function a$D_r_node(box, nd_ix) result(r)
    type(box$D_t), intent(in) :: box
    integer, intent(in)      :: nd_ix($D)
    real(dp)                 :: r($D)
    r = box%r_min + (nd_ix-1) * box%dr
  end function a$D_r_node

  !> Map true -> 1, false -> 0
  elemental integer function l2i(my_bool)
    logical, intent(in) :: my_bool
    if (my_bool) then
       l2i = 1
    else
       l2i = 0
    end if
  end function l2i

  !> Set cc(..., iv) = 0
  subroutine a$D_box_clear_cc(box, iv)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: iv
#if $D == 2
    box%cc(:,:, iv) = 0
#elif $D == 3
    box%cc(:,:,:, iv) = 0
#endif
  end subroutine a$D_box_clear_cc

  !> Add cc(..., iv_from) to box%cc(..., iv_to)
  subroutine a$D_box_add_cc(box, iv_from, iv_to)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: iv_from, iv_to
#if $D == 2
    box%cc(:,:, iv_to) = box%cc(:,:, iv_to) + box%cc(:,:, iv_from)
#elif $D == 3
    box%cc(:,:,:, iv_to) = box%cc(:,:,:, iv_to) + box%cc(:,:,:, iv_from)
#endif
  end subroutine a$D_box_add_cc

  !> Subtract cc(..., iv_from) from box%cc(..., iv_to)
  subroutine a$D_box_sub_cc(box, iv_from, iv_to)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: iv_from, iv_to
#if $D == 2
    box%cc(:,:, iv_to) = box%cc(:,:, iv_to) - box%cc(:,:, iv_from)
#elif $D == 3
    box%cc(:,:,:, iv_to) = box%cc(:,:,:, iv_to) - box%cc(:,:,:, iv_from)
#endif
  end subroutine a$D_box_sub_cc

  !> Multipy cc(..., iv) with a
  subroutine a$D_box_times_cc(box, a, iv)
    type(box$D_t), intent(inout) :: box
    real(dp), intent(in)        :: a
    integer, intent(in)         :: iv
#if $D == 2
    box%cc(:,:, iv) = a * box%cc(:,:, iv)
#elif $D == 3
    box%cc(:,:,:, iv) = a * box%cc(:,:,:, iv)
#endif
  end subroutine a$D_box_times_cc

  !> Set cc(..., iv_b) = a * cc(..., iv_a) + b * cc(..., iv_b)
  subroutine a$D_box_lincomb_cc(box, a, iv_a, b, iv_b)
    type(box$D_t), intent(inout) :: box
    real(dp), intent(in)        :: a, b
    integer, intent(in)         :: iv_a, iv_b
#if $D == 2
    box%cc(:,:, iv_b) = a * box%cc(:,:, iv_a) + b * box%cc(:,:, iv_b)
#elif $D == 3
    box%cc(:,:,:, iv_b) = a * box%cc(:,:,:, iv_a) + b * box%cc(:,:,:, iv_b)
#endif
  end subroutine a$D_box_lincomb_cc

  !> Copy cc(..., iv_from) to box%cc(..., iv_to)
  subroutine a$D_box_copy_cc(box, iv_from, iv_to)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: iv_from, iv_to
#if $D == 2
    box%cc(:,:, iv_to) = box%cc(:,:, iv_from)
#elif $D == 3
    box%cc(:,:,:, iv_to) = box%cc(:,:,:, iv_from)
#endif
  end subroutine a$D_box_copy_cc

  !> Copy cc(..., iv_from) to box%cc(..., iv_to) for all ids
  subroutine a$D_boxes_copy_cc(boxes, ids, iv_from, iv_to)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:), iv_from, iv_to
    integer                     :: i

    !$omp parallel do
    do i = 1, size(ids)
       call a$D_box_copy_cc(boxes(ids(i)), iv_from, iv_to)
    end do
    !$omp end parallel do
  end subroutine a$D_boxes_copy_cc

  !> Copy cc(..., iv_from) to box%cc(..., iv_to) for full tree
  subroutine a$D_tree_copy_cc(tree, iv_from, iv_to)
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)       :: iv_from, iv_to
    integer                   :: lvl

    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       call a$D_boxes_copy_cc(tree%boxes, tree%lvls(lvl)%ids, iv_from, iv_to)
    end do
  end subroutine a$D_tree_copy_cc

  !> Find maximum value of cc(..., iv). Ghost cells are not used.
  subroutine a$D_tree_max_cc(tree, iv, cc_max)
    type(a$D_t), intent(in) :: tree
    integer, intent(in)    :: iv
    real(dp), intent(out)  :: cc_max
    real(dp)               :: tmp, my_max
    integer                :: i, id, lvl, nc

    my_max = -huge(1.0_dp)

    !$omp parallel reduction(max: my_max) private(lvl, i, id, nc, tmp)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          nc = tree%boxes(id)%n_cell
#if $D == 2
          tmp = maxval(tree%boxes(id)%cc(1:nc, 1:nc, iv))
#elif $D == 3
          tmp = maxval(tree%boxes(id)%cc(1:nc, 1:nc, 1:nc, iv))
#endif
          if (tmp > my_max) my_max = tmp
       end do
       !$omp end do
    end do
    !$omp end parallel

    cc_max = my_max
  end subroutine a$D_tree_max_cc

  !> Copy fx/fy/fz(..., iv_from) to fx/fy/fz(..., iv_to)
  subroutine a$D_box_copy_fc(box, iv_from, iv_to)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: iv_from, iv_to
#if $D == 2
    box%fx(:,:, iv_to) = box%fx(:,:, iv_from)
    box%fy(:,:, iv_to) = box%fy(:,:, iv_from)
#elif $D == 3
    box%fx(:,:,:, iv_to) = box%fx(:,:,:, iv_from)
    box%fy(:,:,:, iv_to) = box%fy(:,:,:, iv_from)
    box%fz(:,:,:, iv_to) = box%fz(:,:,:, iv_from)
#endif
  end subroutine a$D_box_copy_fc

  !> Copy fx/fy/fz(..., iv_from) to fx/fy/fz(..., iv_to) for all ids
  subroutine a$D_boxes_copy_fc(boxes, ids, iv_from, iv_to)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:), iv_from, iv_to
    integer                     :: i

    !$omp parallel do
    do i = 1, size(ids)
       call a$D_box_copy_fc(boxes(ids(i)), iv_from, iv_to)
    end do
    !$omp end parallel do
  end subroutine a$D_boxes_copy_fc

  !> Copy fx/fy/fz(..., iv_from) to fx/fy/fz(..., iv_to) for full tree
  subroutine a$D_tree_copy_fc(tree, iv_from, iv_to)
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)       :: iv_from, iv_to
    integer                   :: lvl

    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       call a$D_boxes_copy_fc(tree%boxes, tree%lvls(lvl)%ids, iv_from, iv_to)
    end do
  end subroutine a$D_tree_copy_fc

  !> Zeroth-order prolongation to children.
  subroutine a$D_prolong0_from(boxes, id, iv, fill_gc)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Box whose children we will fill
    integer, intent(in)         :: iv        !< Variable that is filled
    logical, intent(in)         :: fill_gc   !< Also fill ghost cells?
    integer                     :: dgc, nc, i_c, c_id, ix_offset($D)
    integer                     :: i, j, i_c1, j_c1
#if $D == 3
    integer                     :: k, k_c1
#endif

    nc  = boxes(id)%n_cell
    dgc = l2i(fill_gc)

    do i_c = 1, a$D_num_children
       c_id = boxes(id)%children(i_c)
       if (c_id == a5_no_box) cycle

       ! Offset of child w.r.t. parent
       ix_offset = a$D_ch_dix(:, i_c) * ishft(nc, -1)

       ! In these loops, we calculate the closest coarse index (_c1)
#if $D == 2
       do j = 1-dgc, nc+dgc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          do i = 1-dgc, nc+dgc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2

             boxes(c_id)%cc(i, j, iv) = boxes(id)%cc(i_c1, j_c1, iv)
          end do
       end do
#elif $D == 3
       do k = 1-dgc, nc+dgc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          do j = 1-dgc, nc+dgc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             do i = 1-dgc, nc+dgc
                i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2

                boxes(c_id)%cc(i, j, k, iv) = &
                     boxes(id)%cc(i_c1, j_c1, k_c1, iv)
             end do
          end do
       end do
#endif
    end do
  end subroutine a$D_prolong0_from

  !> Linear prolongation to children. We use 2-1-1 interpolation (2d) and
  !> 1-1-1-1 interpolation (3D), which do not require corner ghost cells.
  subroutine a$D_prolong1_from(boxes, id, iv, fill_gc)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Box whose children we will fill
    integer, intent(in)         :: iv        !< Variable that is filled
    logical, intent(in)         :: fill_gc   !< Also fill ghost cells?
    integer                     :: dgc, nc, i_c, c_id, ix_offset($D)
    integer                     :: i, j, i_c1, i_c2, j_c1, j_c2
#if $D == 3
    integer                     :: k, k_c1, k_c2
#endif

    nc  = boxes(id)%n_cell
    dgc = l2i(fill_gc)


    do i_c = 1, a$D_num_children
       c_id = boxes(id)%children(i_c)

       ! Offset of child w.r.t. parent
       ix_offset = a$D_ch_dix(:, i_c) * ishft(nc, -1)

       ! In these loops, we calculate the closest coarse index (_c1), and the
       ! one-but-closest (_c2). The fine cell lies in between.
#if $D == 2
       do j = 1-dgc, nc+dgc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          do i = 1-dgc, nc+dgc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1

             boxes(c_id)%cc(i, j, iv) = &
                  0.5_dp * boxes(id)%cc(i_c1, j_c1, iv) + &
                  0.25_dp * boxes(id)%cc(i_c2, j_c1, iv) + &
                  0.25_dp * boxes(id)%cc(i_c1, j_c2, iv)
          end do
       end do
#elif $D == 3
       do k = 1-dgc, nc+dgc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
          do j = 1-dgc, nc+dgc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
             do i = 1-dgc, nc+dgc
                i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
                i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1

                boxes(c_id)%cc(i, j, k, iv) = 0.25_dp * ( &
                     boxes(id)%cc(i_c1, j_c1, k_c1, iv) + &
                     boxes(id)%cc(i_c2, j_c1, k_c1, iv) + &
                     boxes(id)%cc(i_c1, j_c2, k_c1, iv) + &
                     boxes(id)%cc(i_c1, j_c1, k_c2, iv))
             end do
          end do
       end do
#endif
    end do
  end subroutine a$D_prolong1_from

  !> Partial prolongation to a child (from parent) using injection (simply copy value)
  subroutine a$D_prolong0_to(boxes, id, iv, lo_a, hi_a)
    type(box$D_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Id of child
    integer, intent(in)           :: iv       !< Variable to fill
    integer, intent(in), optional :: lo_a($D) !< Min cell index at child
    integer, intent(in), optional :: hi_a($D) !< Max cell index at child
    integer                       :: nc, p_id, ix_offset($D)
    integer                       :: i, j, i_c1, j_c1, lo($D), hi($D)
#if $D == 3
    integer                       :: k, k_c1
#endif

    nc   = boxes(id)%n_cell
    p_id = boxes(id)%parent
    lo   = 1; if (present(lo_a)) lo = lo_a
    hi   = nc; if (present(hi_a)) hi = hi_a

    ! Offset of child w.r.t. parent
    ix_offset = (boxes(id)%ix - 2*boxes(p_id)%ix + 1) * ishft(nc, -1)

#if $D == 2
    do j = lo(2), hi(2)
       j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
       do i = lo(1), hi(1)
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          boxes(id)%cc(i, j, iv) = boxes(p_id)%cc(i_c1, j_c1, iv)
       end do
    end do
#elif $D == 3
    do k = lo(3), hi(3)
       k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
       do j = lo(2), hi(2)
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          do i = lo(1), hi(1)
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             boxes(id)%cc(i, j, k, iv) = boxes(p_id)%cc(i_c1, j_c1, k_c1, iv)
          end do
       end do
    end do
#endif
  end subroutine a$D_prolong0_to

  !> Partial prolongation to a child (from parent) using linear interpolation.
  !> We use 2-1-1 interpolation (2D) and 1-1-1-1 interpolation (3D) which do not
  !> need corner ghost cells.
  subroutine a$D_prolong1_to(boxes, id, iv, lo_a, hi_a)
    type(box$D_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Id of child
    integer, intent(in), optional :: lo_a($D) !< Min cell index at child
    integer, intent(in), optional :: hi_a($D) !< Max cell index at child
    integer, intent(in)           :: iv       !< Variable to fill
    integer                       :: nc, p_id, ix_offset($D), lo($D), hi($D)
    integer                       :: i, j, i_c1, i_c2, j_c1, j_c2
#if $D == 3
    integer                       :: k, k_c1, k_c2
#endif

    nc        = boxes(id)%n_cell
    p_id      = boxes(id)%parent
    lo   = 1; if (present(lo_a)) lo = lo_a
    hi   = nc; if (present(hi_a)) hi = hi_a

    ! Offset of child w.r.t. parent
    ix_offset = (boxes(id)%ix - 2 * boxes(p_id)%ix + 1) * ishft(nc, -1)

    ! In these loops, we calculate the closest coarse index (i_c1, j_c1), and
    ! the one-but-closest (i_c2, j_c2). The fine cell lies in between.
#if $D == 2
    do j = lo(2), hi(2)
       j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
       do i = lo(1), hi(1)
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
          boxes(id)%cc(i, j, iv) = &
               0.5_dp * boxes(p_id)%cc(i_c1, j_c1, iv) + &
               0.25_dp * boxes(p_id)%cc(i_c2, j_c1, iv) + &
               0.25_dp * boxes(p_id)%cc(i_c1, j_c2, iv)
       end do
    end do
#elif $D == 3
    do k = lo(3), hi(3)
       k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
       k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
       do j = lo(2), hi(2)
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          do i = lo(1), hi(1)
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             boxes(id)%cc(i, j, k, iv) = 0.25_dp * ( &
                  boxes(p_id)%cc(i_c1, j_c1, k_c1, iv) + &
                  boxes(p_id)%cc(i_c2, j_c1, k_c1, iv) + &
                  boxes(p_id)%cc(i_c1, j_c2, k_c1, iv) + &
                  boxes(p_id)%cc(i_c1, j_c1, k_c2, iv))
          end do
       end do
    end do
#endif
  end subroutine a$D_prolong1_to

  !> Restrict the children of a box to the box (e.g., in 2D, average the values
  !> at the four children to get the value for the parent)
  subroutine a$D_restrict_to_box(boxes, id, iv, i_to)
    type(box$D_t), intent(inout)   :: boxes(:) !< List of all the boxes
    integer, intent(in)           :: id       !< Box whose children will be restricted to it
    integer, intent(in)           :: iv       !< Variable to restrict
    integer, intent(in), optional :: i_to    !< Destination (if /= iv)
    integer                       :: nc, i_c, c_id

    nc = boxes(id)%n_cell
    do i_c = 1, a$D_num_children
       c_id = boxes(id)%children(i_c)
       if (c_id == a5_no_box) cycle
       call a$D_restrict_box(boxes(c_id), boxes(id), iv, i_to)
    end do
  end subroutine a$D_restrict_to_box

  !> Restrict the children of boxes ids(:) to them.
  subroutine a$D_restrict_to_boxes(boxes, ids, iv, i_to)
    type(box$D_t), intent(inout)   :: boxes(:) !< List of all the boxes
    integer, intent(in)           :: ids(:)   !< Boxes whose children will be restricted to it
    integer, intent(in)           :: iv       !< Variable to restrict
    integer, intent(in), optional :: i_to    !< Destination (if /= iv)
    integer                       :: i

    !$omp parallel do
    do i = 1, size(ids)
       call a$D_restrict_to_box(boxes, ids(i), iv, i_to)
    end do
    !$omp end parallel do
  end subroutine a$D_restrict_to_boxes

  !> Restrict variables iv to all parent boxes, from the highest to the lowest level
  subroutine a$D_restrict_tree(tree, iv, i_to)
    type(a$D_t), intent(inout)     :: tree  !< Tree to restrict on
    integer, intent(in)           :: iv    !< Variable to restrict
    integer, intent(in), optional :: i_to !< Destination (if /= iv)
    integer                       :: lvl

    do lvl = tree%max_lvl-1, lbound(tree%lvls, 1), -1
       call a$D_restrict_to_boxes(tree%boxes, tree%lvls(lvl)%parents, iv, i_to)
    end do
  end subroutine a$D_restrict_tree

  !> Restriction of child box (box_c) to another box (box_p), typically its
  !> parent. Note that ix_offset is used to specify to which part of box_p we
  !> should restrict box_c.
  subroutine a$D_restrict_box(box_c, box_p, iv, i_to)
    type(box$D_t), intent(in)      :: box_c         !< Child box to restrict
    type(box$D_t), intent(inout)   :: box_p         !< Parent box to restrict to
    integer, intent(in)           :: iv            !< Variable to restrict
    integer, intent(in), optional :: i_to         !< Destination (if /= iv)
    integer                       :: i, j, i_c1, j_c1, i_dest
    integer                       :: nc, ix_offset($D)
#if $D == 3
    integer                       :: k, k_c1
#endif

    nc        = box_c%n_cell
    ix_offset = a$D_get_child_offset(box_c)

    if (present(i_to)) then
       i_dest = i_to
    else
       i_dest = iv
    end if

#if $D == 2
    do j = 1, nc, 2
       j_c1 = ix_offset(2) + ishft(j+1, -1)  ! (j+1)/2
       do i = 1, nc, 2
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          box_p%cc(i_c1, j_c1, i_dest) = 0.25_dp * sum(box_c%cc(i:i+1, j:j+1, iv))
       end do
    end do
#elif $D == 3
    do k = 1, nc, 2
       k_c1 = ix_offset(3) + ishft(k+1, -1)  ! (k+1)/2
       do j = 1, nc, 2
          j_c1 = ix_offset(2) + ishft(j+1, -1)  ! (j+1)/2
          do i = 1, nc, 2
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             box_p%cc(i_c1, j_c1, k_c1, i_dest) = 0.125_dp * &
                  sum(box_c%cc(i:i+1, j:j+1, k:k+1, iv))
          end do
       end do
    end do
#endif
  end subroutine a$D_restrict_box

  !> Fill ghost cells for variables iv on the sides of all boxes, using
  !> subr_rb on refinement boundaries and subr_bc on physical boundaries
  subroutine a$D_gc_sides(tree, iv, subr_rb, subr_bc)
    type(a$D_t), intent(inout) :: tree !< Tree to fill ghost cells on
    integer, intent(in)       :: iv !< Variable for which ghost cells are set
    procedure(a$D_subr_gc)     :: subr_rb !< Procedure called at refinement boundaries
    procedure(a$D_subr_gc)     :: subr_bc    !< Procedure called at physical boundaries
    integer                   :: lvl, i, id

    !$omp parallel private(lvl, i, id)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call a$D_gc_box_sides(tree%boxes, id, iv, subr_rb, subr_bc)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine a$D_gc_sides

  !> Fill ghost cells for variables iv on the sides of a box, using
  !> subr_rb on refinement boundaries and subr_bc on physical boundaries
  subroutine a$D_gc_box_sides(boxes, id, iv, subr_rb, subr_bc)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all the boxes
    integer, intent(in)         :: id !< Id of box for which we set ghost cells
    integer, intent(in)         :: iv !< Variable for which ghost cells are set
    procedure(a$D_subr_gc)      :: subr_rb !< Procedure called at refinement boundaries
    procedure(a$D_subr_gc)      :: subr_bc    !< Procedure called at physical boundaries
    integer                     :: nb

    do nb = 1, a$D_num_neighbors
       if (boxes(id)%neighbors(nb) > a5_no_box) then
          call a$D_gc_side_from_nb(boxes, id, nb, iv)
       else if (boxes(id)%neighbors(nb) == a5_no_box) then
          call subr_rb(boxes, id, nb, iv)
       else
          call subr_bc(boxes, id, nb, iv)
       end if
    end do
  end subroutine a$D_gc_box_sides

  !> Linear interpolation (2-1-1 type) from parent to fill ghost cells on the
  !> side of a box
  subroutine a$D_sides_prolong1(boxes, id, nb, iv)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc

    nc = boxes(id)%n_cell

    select case (nb)
#if $D == 2
    case (a2_nb_lx)
       call a2_prolong1_to(boxes, id, iv, [0, 1], [0, nc])
    case (a2_nb_hx)
       call a2_prolong1_to(boxes, id, iv, [nc+1, 1], [nc+1, nc])
    case (a2_nb_ly)
       call a2_prolong1_to(boxes, id, iv, [1, 0], [nc, 0])
    case (a2_nb_hy)
       call a2_prolong1_to(boxes, id, iv, [1, nc+1], [nc, nc+1])
#elif $D == 3
    case (a3_nb_lx)
       call a3_prolong1_to(boxes, id, iv, [0, 1, 1], [0, nc, nc])
    case (a3_nb_hx)
       call a3_prolong1_to(boxes, id, iv, [nc+1, 1, 1], [nc+1, nc, nc])
    case (a3_nb_ly)
       call a3_prolong1_to(boxes, id, iv, [1, 0, 1], [nc, 0, nc])
    case (a3_nb_hy)
       call a3_prolong1_to(boxes, id, iv, [1, nc+1, 1], [nc, nc+1, nc])
    case (a3_nb_lz)
       call a3_prolong1_to(boxes, id, iv, [1, 1, 0], [nc, nc, 0])
    case (a3_nb_hz)
       call a3_prolong1_to(boxes, id, iv, [1, 1, nc+1], [nc, nc, nc+1])
#endif
    end select
  end subroutine a$D_sides_prolong1

  !> Interpolate between fine points and coarse neighbors to fill ghost cells
  !> near refinement boundaries
  subroutine a$D_sides_interp(boxes, id, nb, iv)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc, ix, ix_c, dix, i, di, j, dj
    integer                     :: i_c1, i_c2, j_c1, j_c2, p_nb_id
    integer                     :: p_id, ix_offset($D)
#if $D == 3
    integer                     :: k_c1, k_c2, k, dk
    real(dp), parameter         :: f1=1/32.0_dp, f3=3*f1, f9=9*f1
#endif

    nc        = boxes(id)%n_cell
    p_id      = boxes(id)%parent
    p_nb_id   = boxes(p_id)%neighbors(nb)
    ix_offset = (boxes(id)%ix - 2 * boxes(p_id)%ix + 1) * ishft(nc, -1)

    if (a$D_nb_low(nb)) then
       ix = 0
       dix = 1
       ix_c = nc
    else
       ix = nc+1
       dix = -1
       ix_c = 1
    end if

    select case (a$D_nb_dim(nb))
#if $D == 2
    case (1)
       i = ix
       di = dix
       i_c1 = ix_c

       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          boxes(id)%cc(i, j, iv) = &
               0.375_dp * boxes(p_nb_id)%cc(i_c1, j_c1, iv) + &
               0.125_dp * boxes(p_nb_id)%cc(i_c1, j_c2, iv) + &
               0.5_dp * boxes(id)%cc(i+di, j, iv)
       end do
    case (2)
       j = ix
       dj = dix
       j_c1 = ix_c

       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
          boxes(id)%cc(i, j, iv) = &
               0.375_dp * boxes(p_nb_id)%cc(i_c1, j_c1, iv) + &
               0.125_dp * boxes(p_nb_id)%cc(i_c2, j_c1, iv) + &
               0.5_dp * boxes(id)%cc(i, j+dj, iv)
       end do
#elif $D==3
    case (1)
       i = ix
       di = dix
       i_c1 = ix_c

       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)          ! even: +1, odd: -1
          do j = 1, nc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
             boxes(id)%cc(i, j, k, iv) = &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c1, iv) + &
                  0.125_dp * boxes(p_nb_id)%cc(i_c1, j_c2, k_c1, iv) + &
                  0.125_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c2, iv) + &
                  0.5_dp * boxes(id)%cc(i+di, j, k, iv)
          end do
       end do
    case (2)
       j = ix
       dj = dix
       j_c1 = ix_c

       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             boxes(id)%cc(i, j, k, iv) = &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c1, iv) + &
                  0.125_dp * boxes(p_nb_id)%cc(i_c2, j_c1, k_c1, iv) + &
                  0.125_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c2, iv) + &
                  0.5_dp * boxes(id)%cc(i, j+dj, k, iv)
          end do
       end do
    case (3)
       k = ix
       dk = dix
       k_c1 = ix_c

       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             boxes(id)%cc(i, j, k, iv) = &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c1, iv) + &
                  0.125_dp * boxes(p_nb_id)%cc(i_c1, j_c2, k_c1, iv) + &
                  0.125_dp * boxes(p_nb_id)%cc(i_c2, j_c1, k_c1, iv) + &
                  0.5_dp * boxes(id)%cc(i, j, k+dk, iv)
          end do
       end do
#endif
    end select

  end subroutine a$D_sides_interp

  !> Fill ghost cells near refinement boundaries which preserves diffusive fluxes.
  !>
  !> Basically, we extrapolate from the fine cells to a corner point, and then
  !> take the average between this corner point and a coarse neighbor to fill
  !> ghost cells for the fine cells.
  !> @todo Remove corner points in 3D
  subroutine a$D_sides_extrap(boxes, id, nb, iv)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc, ix, dix, i, di, j, dj
#if $D == 3
    integer                     :: k, dk
#endif

    nc = boxes(id)%n_cell

    if (a$D_nb_low(nb)) then
       ix = 1
       dix = 1
    else
       ix = nc
       dix = -1
    end if

    select case (a$D_nb_dim(nb))
#if $D == 2
    case (1)
       i = ix
       di = dix
       call a2_prolong0_to(boxes, id, iv, [i-di, 1], [i-di, nc])

       do j = 1, nc
          dj = -1 + 2 * iand(j, 1)
          boxes(id)%cc(i-di, j, iv) = 0.5_dp * boxes(id)%cc(i-di, j, iv) + &
               boxes(id)%cc(i, j, iv) - 0.25_dp * (boxes(id)%cc(i+di, j, iv) &
               + boxes(id)%cc(i, j+dj, iv))
       end do

    case (2)
       j = ix
       dj = dix
       call a2_prolong0_to(boxes, id, iv, [1, j-dj], [nc, j-dj])

       do i = 1, nc
          di = -1 + 2 * iand(i, 1)
          boxes(id)%cc(i, j-dj, iv) = 0.5_dp * boxes(id)%cc(i, j-dj, iv) + &
               boxes(id)%cc(i, j, iv) - 0.25_dp * (boxes(id)%cc(i+di, j, iv) &
               + boxes(id)%cc(i, j+dj, iv))
       end do
#elif $D == 3
    case (1)
       i = ix
       di = dix
       call a3_prolong0_to(boxes, id, iv, [i-di, 1, 1], [i-di, nc, nc])

       do k = 1, nc
          dk = -1 + 2 * iand(k, 1)
          do j = 1, nc
             dj = -1 + 2 * iand(j, 1)

             boxes(id)%cc(i-di, j, k, iv) = &
                  0.5_dp * boxes(id)%cc(i-di, j, k, iv) + 0.0625_dp * &
                  (27 * boxes(id)%cc(i, j, k, iv) - 9 * boxes(id)%cc(i+di, j, k, iv) &
                  - 9 * boxes(id)%cc(i, j+dj, k, iv) - 9 * boxes(id)%cc(i, j, k+dk, iv) &
                  + 3 * boxes(id)%cc(i+di, j+dj, k, iv) + 3 * boxes(id)%cc(i+di, j, k+dk, iv) &
                  + 3 * boxes(id)%cc(i, j+dj, k+dk, iv) - boxes(id)%cc(i+di, j+dj, k+dk, iv))
          end do
       end do

    case (2)
       j = ix
       dj = dix
       call a3_prolong0_to(boxes, id, iv, [1, j-dj, 1], [nc, j-dj, nc])

       do k = 1, nc
          dk = -1 + 2 * iand(k, 1)
          do i = 1, nc
             di = -1 + 2 * iand(i, 1)

             boxes(id)%cc(i, j-dj, k, iv) = &
                  0.5_dp * boxes(id)%cc(i, j-dj, k, iv) + 0.0625_dp * &
                  (27 * boxes(id)%cc(i, j, k, iv) - 9 * boxes(id)%cc(i+di, j, k, iv) &
                  - 9 * boxes(id)%cc(i, j+dj, k, iv) - 9 * boxes(id)%cc(i, j, k+dk, iv) &
                  + 3 * boxes(id)%cc(i+di, j+dj, k, iv) + 3 * boxes(id)%cc(i+di, j, k+dk, iv) &
                  + 3 * boxes(id)%cc(i, j+dj, k+dk, iv) - boxes(id)%cc(i+di, j+dj, k+dk, iv))
          end do
       end do

    case (3)
       k = ix
       dk = dix
       call a3_prolong0_to(boxes, id, iv, [1, 1, k-dk], [nc, nc, k-dk])

       do j = 1, nc
          dj = -1 + 2 * iand(j, 1)
          do i = 1, nc
             di = -1 + 2 * iand(i, 1)

             boxes(id)%cc(i, j, k-dk, iv) = &
                  0.5_dp * boxes(id)%cc(i, j, k-dk, iv) + 0.0625_dp * &
                  (27 * boxes(id)%cc(i, j, k, iv) - 9 * boxes(id)%cc(i+di, j, k, iv) &
                  - 9 * boxes(id)%cc(i, j+dj, k, iv) - 9 * boxes(id)%cc(i, j, k+dk, iv) &
                  + 3 * boxes(id)%cc(i+di, j+dj, k, iv) + 3 * boxes(id)%cc(i+di, j, k+dk, iv) &
                  + 3 * boxes(id)%cc(i, j+dj, k+dk, iv) - boxes(id)%cc(i+di, j+dj, k+dk, iv))
          end do
       end do
#endif
    end select
  end subroutine a$D_sides_extrap

  !> Fill values on the side of a box from a neighbor nb
  subroutine a$D_gc_side_from_nb(boxes, id, nb, iv)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell / neighbor direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc, nb_id

    nc    = boxes(id)%n_cell
    nb_id = boxes(id)%neighbors(nb)

    select case (nb)
#if $D == 2
    case (a2_nb_lx)
       boxes(id)%cc(0, 1:nc, iv)    = boxes(nb_id)%cc(nc, 1:nc, iv)
    case (a2_nb_hx)
       boxes(id)%cc(nc+1, 1:nc, iv) = boxes(nb_id)%cc(1, 1:nc, iv)
    case (a2_nb_ly)
       boxes(id)%cc(1:nc, 0, iv)    = boxes(nb_id)%cc(1:nc, nc, iv)
    case (a2_nb_hy)
       boxes(id)%cc(1:nc, nc+1, iv) = boxes(nb_id)%cc(1:nc, 1, iv)
#elif $D == 3
    case (a3_nb_lx)
       boxes(id)%cc(0, 1:nc, 1:nc, iv)    = boxes(nb_id)%cc(nc, 1:nc, 1:nc, iv)
    case (a3_nb_hx)
       boxes(id)%cc(nc+1, 1:nc, 1:nc, iv) = boxes(nb_id)%cc(1, 1:nc, 1:nc, iv)
    case (a3_nb_ly)
       boxes(id)%cc(1:nc, 0, 1:nc, iv)    = boxes(nb_id)%cc(1:nc, nc, 1:nc, iv)
    case (a3_nb_hy)
       boxes(id)%cc(1:nc, nc+1, 1:nc, iv) = boxes(nb_id)%cc(1:nc, 1, 1:nc, iv)
    case (a3_nb_lz)
       boxes(id)%cc(1:nc, 1:nc, 0, iv)    = boxes(nb_id)%cc(1:nc, 1:nc, nc, iv)
    case (a3_nb_hz)
       boxes(id)%cc(1:nc, 1:nc, nc+1, iv) = boxes(nb_id)%cc(1:nc, 1:nc, 1, iv)
#endif
    end select
  end subroutine a$D_gc_side_from_nb

  !> Restrict fluxes from children to parents on refinement boundaries
  subroutine a$D_consistent_fluxes(tree, f_ixs)
    use omp_lib
    type(a$D_t), intent(inout) :: tree    !< Tree to operate on
    integer, intent(in)       :: f_ixs(:) !< Indices of the fluxes
    integer                   :: lvl, i, id, nb, nb_id

    !$omp parallel private(lvl, i, id, nb, nb_id)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl-1
       !$omp do
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          do nb = 1, a$D_num_neighbors
             nb_id = tree%boxes(id)%neighbors(nb)
             if (nb_id < a5_no_box) cycle ! Boundary condition
             if (.not. a$D_has_children(tree%boxes(nb_id))) &
                  call a$D_flux_from_children(tree%boxes, id, nb, f_ixs)
          end do
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine a$D_consistent_fluxes

  !> The neighbor nb has no children and id does, so get flux from children for
  !> consisency at refinement boundary.
  !> @todo: TEST
  subroutine a$D_flux_from_children(boxes, id, nb, f_ixs)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all the boxes
    integer, intent(in)         :: id        !< Id of box for which we set fluxes
    integer, intent(in)         :: nb        !< Direction in which fluxes are set
    integer, intent(in)         :: f_ixs(:)  !< Indices of the fluxes
    integer                     :: nc, nch, c_id, i_ch, i, ic, d, ioff($D)
    integer                     :: n_chnb

    nc     = boxes(id)%n_cell
    nch    = ishft(nc, -1) ! nc/2
    d      = a$D_nb_dim(nb)
    n_chnb = 2**($D-1)

    if (a$D_nb_low(nb)) then
       i = 1
    else
       i = nc
    end if

    select case (d)
#if $D == 2
    case (1)
       do ic = 1, n_chnb
          i_ch = a2_ch_adj_nb(ic, nb)
          c_id = boxes(id)%children(i_ch)
          ioff = nch*a2_ch_dix(:, i_ch)
          boxes(id)%fx(i, ioff(2)+1:ioff(2)+nch, f_ixs) = 0.5_dp * ( &
               boxes(c_id)%fx(i, 1:nc:2, f_ixs) + &
               boxes(c_id)%fx(i, 2:nc:2, f_ixs))
       end do
    case (2)
       do ic = 1, n_chnb
          i_ch = a2_ch_adj_nb(ic, nb)
          c_id = boxes(id)%children(i_ch)
          ioff = nch*a2_ch_dix(:, i_ch)
          boxes(id)%fx(ioff(1)+1:ioff(1)+nch, i, f_ixs) = 0.5_dp * ( &
               boxes(c_id)%fx(1:nc:2, i, f_ixs) + &
               boxes(c_id)%fx(2:nc:2, i, f_ixs))
       end do
#elif $D == 3
    case (1)
       do ic = 1, n_chnb
          i_ch = a3_ch_adj_nb(ic, nb)
          c_id = boxes(id)%children(i_ch)
          ioff = nch*a3_ch_dix(:, i_ch)
          boxes(id)%fx(i, ioff(2)+1:ioff(2)+nch, &
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
          boxes(id)%fx(ioff(1)+1:ioff(1)+nch, i, &
               ioff(3)+1:ioff(3)+nch, f_ixs) = 0.25_dp * ( &
               boxes(c_id)%fx(1:nc:2, i, 1:nc:2, f_ixs) + &
               boxes(c_id)%fx(2:nc:2, i, 1:nc:2, f_ixs) + &
               boxes(c_id)%fx(1:nc:2, i, 2:nc:2, f_ixs) + &
               boxes(c_id)%fx(2:nc:2, i, 2:nc:2, f_ixs))
       end do
    case (3)
       do ic = 1, n_chnb
          i_ch = a3_ch_adj_nb(ic, nb)
          c_id = boxes(id)%children(i_ch)
          ioff = nch*a3_ch_dix(:, i_ch)
          boxes(id)%fx(ioff(1)+1:ioff(1)+nch, &
               ioff(2)+1:ioff(2)+nch, i, f_ixs) = 0.25_dp * ( &
               boxes(c_id)%fx(1:nc:2, 1:nc:2, i, f_ixs) + &
               boxes(c_id)%fx(2:nc:2, 1:nc:2, i, f_ixs) + &
               boxes(c_id)%fx(1:nc:2, 2:nc:2, i, f_ixs) + &
               boxes(c_id)%fx(2:nc:2, 2:nc:2, i, f_ixs))
       end do
#endif
    end select
  end subroutine a$D_flux_from_children

  !> Write the cell centered data of a tree to a vtk unstructured file
  subroutine a$D_write_vtk(tree, filename, cc_names, n_cycle, time, ivs)
    use m_vtk
    type(a$D_t), intent(in)       :: tree        !< Tree to write out
    character(len=*)              :: filename    !< Filename for the vtk file
    character(len=*)              :: cc_names(:) !< Names of the cell-centered variables
    integer, intent(in)           :: n_cycle     !< Cycle-number for vtk file (counter)
    real(dp), intent(in)          :: time        !< Time for output file
    integer, intent(in), optional :: ivs(:)      !< Oncly include these variables
    integer                       :: lvl, bc, bn, n, n_cells, n_nodes
    integer                       :: ig, i, j, id, n_ix, c_ix, n_grids
    integer                       :: cell_ix, node_ix
    integer, parameter            :: n_ch = a$D_num_children
    integer                       :: nodes_per_box, cells_per_box
    real(dp), allocatable         :: coords(:), cc_vars(:,:)
    integer, allocatable          :: offsets(:), connects(:)
    integer, allocatable          :: cell_types(:), ivs_used(:)
    type(vtk_t)                   :: vtkf
#if $D == 3
    integer                       :: k, bn2
#endif

    if (present(ivs)) then
       if (maxval(ivs) > tree%n_var_cell .or. &
            minval(ivs) < 1) stop "a$D_write_vtk: wrong indices given (ivs)"
       if (size(ivs) /= size(cc_names)) &
            stop "a$D_write_vtk: size(cc_names) /= size(ivs)"
       ivs_used = ivs
    else
       if (size(cc_names) /= tree%n_var_cell) &
            stop "a$D_write_vtk: size(cc_names) /= n_var_cell"
       ivs_used = [(i, i = 1, tree%n_var_cell)]
    end if

    bc            = tree%n_cell     ! number of Box Cells
    bn            = tree%n_cell + 1 ! number of Box Nodes
    nodes_per_box = bn**$D
    cells_per_box = bc**$D

    n_grids = 0
    do lvl = 1, tree%max_lvl
       n_grids = n_grids + size(tree%lvls(lvl)%leaves)
    end do
    n_nodes = nodes_per_box * n_grids
    n_cells = cells_per_box * n_grids

    allocate(coords($D * n_nodes))
    allocate(cc_vars(n_cells, size(ivs_used)))
    allocate(offsets(cells_per_box * n_grids))
    allocate(cell_types(cells_per_box * n_grids))
    allocate(connects(n_ch * cells_per_box * n_grids))

#if $D == 2
    cell_types = 8  ! VTK pixel type
#elif $D       == 3
    bn2        = bn**2
    cell_types = 11 ! VTK voxel type
#endif

    ig = 0
    do lvl = 1, tree%max_lvl
       do n = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(n)

          ig = ig + 1
          cell_ix = (ig-1) * cells_per_box
          node_ix = (ig-1) * nodes_per_box

#if $D == 2
          do j = 1, bn
             do i = 1, bn
                n_ix = 2 * (node_ix + (j-1) * bn + i)
                coords(n_ix-1:n_ix) = tree%boxes(id)%r_min + &
                     [i-1,j-1] * tree%boxes(id)%dr
             end do
          end do

          do j = 1, bc
             do i = 1, bc
                ! In vtk, indexing starts at 0, so subtract 1
                n_ix                      = node_ix + (j-1) * bn + i - 1
                c_ix                      = cell_ix + (j-1) * bc + i
                cc_vars(c_ix, :)          = tree%boxes(id)%cc(i, j, ivs_used)
                offsets(c_ix)             = a$D_num_children * c_ix
                connects(n_ch*(c_ix-1)+1:n_ch*c_ix) = [n_ix, n_ix+1, n_ix+bn, n_ix+bn+1]
             end do
          end do
#elif $D == 3
          do k = 1, bn
             do j = 1, bn
                do i = 1, bn
                   n_ix = 3 * (node_ix + (k-1) * bn2 + (j-1) * bn + i)
                   coords(n_ix-2:n_ix) = tree%boxes(id)%r_min + &
                        [i-1,j-1,k-1] * tree%boxes(id)%dr
                end do
             end do
          end do

          do k = 1, bc
             do j = 1, bc
                do i = 1, bc
                   ! In vtk, indexing starts at 0, so subtract 1
                   n_ix                      = node_ix + (k-1) * bn2 + &
                        (j-1) * bn + i - 1
                   c_ix                      = cell_ix + (k-1) * bc**2 + &
                        (j-1) * bc + i
                   cc_vars(c_ix, :)          = tree%boxes(id)%cc(i, j, k, ivs_used)
                   offsets(c_ix)             = 8 * c_ix
                   connects(n_ch*(c_ix-1)+1:n_ch*c_ix) = &
                        [n_ix, n_ix+1, n_ix+bn, n_ix+bn+1, &
                        n_ix+bn2, n_ix+bn2+1, n_ix+bn2+bn, n_ix+bn2+bn+1]
                end do
             end do
          end do
#endif
       end do
    end do

    call vtk_ini_xml(vtkf, trim(filename), 'UnstructuredGrid')
    call vtk_dat_xml(vtkf, "UnstructuredGrid", .true.)
    call vtk_geo_xml(vtkf, coords, n_nodes, n_cells, $D, n_cycle, time)
    call vtk_con_xml(vtkf, connects, offsets, cell_types, n_cells)
    call vtk_dat_xml(vtkf, "CellData", .true.)

    do n = 1, size(ivs_used)
       call vtk_var_r8_xml(vtkf, trim(cc_names(n)), cc_vars(:, n), n_cells)
    end do

    call vtk_dat_xml(vtkf, "CellData", .false.)
    call vtk_geo_xml_close(vtkf)
    call vtk_dat_xml(vtkf, "UnstructuredGrid", .false.)
    call vtk_end_xml(vtkf)
    print *, "Written ", trim(filename), ", n_grids", n_grids
  end subroutine a$D_write_vtk

  subroutine a$D_write_silo(tree, filename, cc_names, n_cycle, time, ivs)
    use m_write_silo
    type(a$D_t), intent(inout)       :: tree
    character(len=*)                :: filename, cc_names(:)
    integer, intent(in)             :: n_cycle
    real(dp), intent(in)            :: time
    integer, intent(in), optional   :: ivs(:)
    character(len=*), parameter     :: grid_name = "gg", var_name  = "vv"
    character(len=*), parameter     :: amr_name  = "amr", dir_name = "data"
    character(len=100), allocatable :: grid_list(:), var_list(:, :)
    integer                         :: lvl, i, id, i_grid, iv, nc, n_grids_max
    integer                         :: n_vars, i0, j0, dbix
    integer                         :: nx, ny, nx_prev, ny_prev, ix, iy
    integer, allocatable            :: ids(:), nb_ids(:), ivs_used(:)
    logical, allocatable            :: box_done(:)
    real(dp)                        :: dr($D), r_min($D)
#if $D == 2
    integer, allocatable            :: box_list(:,:), new_box_list(:, :)
    real(dp), allocatable           :: var_data(:,:,:)
#elif $D == 3
    integer, allocatable            :: box_list(:,:,:), new_box_list(:,:,:)
    real(dp), allocatable           :: var_data(:,:,:,:)
    integer                         :: k0, nz, nz_prev, iz
#endif

    if (present(ivs)) then
       if (maxval(ivs) > tree%n_var_cell .or. &
            minval(ivs) < 1) stop "a$D_write_silo: wrong indices given (ivs)"
       if (size(ivs) /= size(cc_names)) &
            stop "a$D_write_silo: size(cc_names) /= size(ivs)"
       ivs_used = ivs
    else
       if (size(cc_names) /= tree%n_var_cell) &
            stop "a$D_write_silo: size(cc_names) /= n_var_cell"
       ivs_used = [(i, i = 1, tree%n_var_cell)]
    end if

    nc = tree%n_cell
    n_vars = size(ivs_used)
    n_grids_max = 0
    do lvl = 1, tree%max_lvl
       n_grids_max = n_grids_max + size(tree%lvls(lvl)%leaves)
    end do

    allocate(grid_list(n_grids_max))
    allocate(var_list(n_vars, n_grids_max))
    allocate(box_done(tree%max_id))
    box_done = .false.

    call SILO_create_file(filename, dbix)
    call SILO_mkdir(dbix, dir_name)
    i_grid = 0

    do lvl = 1, tree%max_lvl
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          if (box_done(id)) cycle

          i_grid = i_grid + 1

          ! Find largest rectangular box including id and other leaves that
          ! haven't been written yet
#if $D == 2
          allocate(box_list(1,1))
          box_list(1,1) = id
          nx = 1
          ny = 1

          do
             nx_prev = nx
             ny_prev = ny

             ! Check whether we can extend to the -x direction
             ids = box_list(1, :)
             nb_ids = tree%boxes(ids)%neighbors(a2_nb_lx)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     .not. any(a2_has_children(tree%boxes(nb_ids)))) then
                   nx = nx + 1
                   allocate(new_box_list(nx, ny))
                   new_box_list(1, :) = nb_ids
                   new_box_list(2:, :) = box_list
                   box_list = new_box_list
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the +x direction
             ids = box_list(nx, :)
             nb_ids = tree%boxes(ids)%neighbors(a2_nb_hx)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     .not. any(a2_has_children(tree%boxes(nb_ids)))) then
                   nx = nx + 1
                   allocate(new_box_list(nx, ny))
                   new_box_list(1:nx-1, :) = box_list
                   new_box_list(nx, :) = nb_ids
                   box_list = new_box_list
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the -y direction
             ids = box_list(:, 1)
             nb_ids = tree%boxes(ids)%neighbors(a2_nb_ly)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     .not. any(a2_has_children(tree%boxes(nb_ids)))) then
                   ny = ny + 1
                   allocate(new_box_list(nx, ny))
                   new_box_list(:, 1) = nb_ids
                   new_box_list(:, 2:) = box_list
                   box_list = new_box_list
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the +y direction
             ids = box_list(:, ny)
             nb_ids = tree%boxes(ids)%neighbors(a2_nb_hy)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     .not. any(a2_has_children(tree%boxes(nb_ids)))) then
                   ny = ny + 1
                   allocate(new_box_list(nx, ny))
                   new_box_list(:, 1:ny-1) = box_list
                   new_box_list(:, ny) = nb_ids
                   box_list = new_box_list
                   deallocate(new_box_list)
                end if
             end if

             if (nx == nx_prev .and. ny == ny_prev) exit
          end do

          allocate(var_data(nx * nc, ny * nc, n_vars))
          do ix = 1, nx
             do iy = 1, ny
                id = box_list(ix, iy)
                box_done(id) = .true.
                i0 = 1 + (ix-1) * nc
                j0 = 1 + (iy-1) * nc
                var_data(i0:i0+nc-1, j0:j0+nc-1, :) = &
                     tree%boxes(id)%cc(1:nc, 1:nc, ivs_used)
             end do
          end do

          id = box_list(1, 1)
          dr = tree%boxes(id)%dr
          r_min = tree%boxes(id)%r_min

          write(grid_list(i_grid), "(A,I0)") dir_name // '/' // grid_name, i_grid
          call SILO_add_grid(dbix, grid_list(i_grid), 2, &
               [nx*nc + 1, ny*nc + 1], r_min, dr)
          do iv = 1, n_vars
             write(var_list(iv, i_grid), "(A,I0)") dir_name // '/' // &
                  trim(cc_names(iv)) // "_", i_grid
             call SILO_add_var(dbix, var_list(iv, i_grid), grid_list(i_grid), &
                  pack(var_data(:, :, iv), .true.), [nx*nc, ny*nc])
          end do

          deallocate(var_data)
          deallocate(box_list)
#elif $D == 3
                    allocate(box_list(1,1,1))
          box_list(1,1,1) = id
          nx = 1
          ny = 1
          nz = 1

          do
             nx_prev = nx
             ny_prev = ny
             nz_prev = nz

             ! Check whether we can extend to the -x direction
             ids = pack(box_list(1, :, :), .true.)

             nb_ids = tree%boxes(ids)%neighbors(a3_nb_lx)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   nx = nx + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(1, :, :) = reshape(nb_ids, [ny, nz])
                   new_box_list(2:, :, :) = box_list
                   box_list = new_box_list
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the +x direction
             ids = pack(box_list(nx, :, :), .true.)
             nb_ids = tree%boxes(ids)%neighbors(a3_nb_hx)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   nx = nx + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(1:nx-1, :, :) = box_list
                   new_box_list(nx, :, :) = reshape(nb_ids, [ny, nz])
                   box_list = new_box_list
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the -y direction
             ids = pack(box_list(:, 1, :), .true.)

             nb_ids = tree%boxes(ids)%neighbors(a3_nb_ly)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   ny = ny + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(:, 1, :) = reshape(nb_ids, [nx, nz])
                   new_box_list(:, 2:, :) = box_list
                   box_list = new_box_list
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the +y direction
             ids = pack(box_list(:, ny, :), .true.)

             nb_ids = tree%boxes(ids)%neighbors(a3_nb_hy)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   ny = ny + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(:, 1:ny-1, :) = box_list
                   new_box_list(:, ny, :) = reshape(nb_ids, [nx, nz])
                   box_list = new_box_list
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the -z direction
             ids = pack(box_list(:, :, 1), .true.)

             nb_ids = tree%boxes(ids)%neighbors(a3_nb_lz)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   nz = nz + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(:, :, 1) = reshape(nb_ids, [nx, ny])
                   new_box_list(:, :, 2:) = box_list
                   box_list = new_box_list
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the +z direction
             ids = pack(box_list(:, :, nz), .true.)

             nb_ids = tree%boxes(ids)%neighbors(a3_nb_hz)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   nz = nz + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(:, :, 1:nz-1) = box_list
                   new_box_list(:, :, nz) = reshape(nb_ids, [nx, ny])
                   box_list = new_box_list
                   deallocate(new_box_list)
                end if
             end if

             if (nx == nx_prev .and. ny == ny_prev .and. nz == nz_prev) exit
          end do

          allocate(var_data(nx * nc, ny * nc, nz * nc, n_vars))
          do iz = 1, nz
             do ix = 1, nx
                do iy = 1, ny
                   id = box_list(ix, iy, iz)
                   box_done(id) = .true.
                   i0 = 1 + (ix-1) * nc
                   j0 = 1 + (iy-1) * nc
                   k0 = 1 + (iz-1) * nc
                   var_data(i0:i0+nc-1, j0:j0+nc-1, k0:k0+nc-1, :) = &
                        tree%boxes(id)%cc(1:nc, 1:nc, 1:nc, ivs_used)
                end do
             end do
          end do

          id = box_list(1, 1, 1)
          dr = tree%boxes(id)%dr
          r_min = tree%boxes(id)%r_min

          write(grid_list(i_grid), "(A,I0)") dir_name // '/' // grid_name, i_grid
          call SILO_add_grid(dbix, grid_list(i_grid), 3, &
               [nx*nc + 1, ny*nc + 1, nz*nc + 1], r_min, dr)
          do iv = 1, n_vars
             write(var_list(iv, i_grid), "(A,I0)") dir_name // '/' // &
                  trim(cc_names(iv)) // "_", i_grid
             call SILO_add_var(dbix, var_list(iv, i_grid), grid_list(i_grid), &
                  pack(var_data(:, :, :, iv), .true.), [nx*nc, ny*nc, nz*nc])
          end do

          deallocate(var_data)
          deallocate(box_list)
#endif

       end do
    end do

    call SILO_set_mmesh_grid(dbix, amr_name, grid_list(1:i_grid), n_cycle, time)
    do iv = 1, n_vars
       call SILO_set_mmesh_var(dbix, trim(cc_names(iv)), amr_name, &
            var_list(iv, 1:i_grid), n_cycle, time)
    end do
    call SILO_close_file(dbix)

    print *, "SILO: Number of grids", i_grid
  end subroutine a$D_write_silo

end module m_afivo_$Dd
