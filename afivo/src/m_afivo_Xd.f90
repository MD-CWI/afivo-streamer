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
  ! For cylindrical coordinates
  integer, parameter :: a2_r_dim = 1
  integer, parameter :: a2_z_dim = 2

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
  ! Opposite of nb_low, but now as integers
  integer, parameter :: a2_nb_hi01(4) = [0, 1, 0, 1]

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
  ! Opposite of nb_low, but now as integers
  integer, parameter :: a3_nb_hi01(6) = [0, 1, 0, 1, 0, 1]
  ! Reverse neighbors
  integer, parameter :: a3_nb_rev(6) = [2, 1, 4, 3, 6, 5]
  ! Direction (dimension) for a neighbor
  integer, parameter :: a3_nb_dim(6) = [1, 1, 2, 2, 3, 3]
#endif

  !> The basic building block of afivo: a box with cell-centered and face
  !> centered data, and information about its position, neighbors, children etc.
  type box$D_t
     integer               :: lvl    !< level of the box
     logical               :: in_use=.false.  !< is the box in use?
     integer               :: tag=a5_init_tag !< for the user
     integer               :: ix($D) !< index in the domain
     integer               :: parent !< index of parent in box list
     !> index of children in box list
     integer               :: children(a$D_num_children)
     !> index of neighbors in box list
     integer               :: neighbors(a$D_num_neighbors)
     integer               :: n_cell    !< number of cells per dimension
     real(dp)              :: dr        !< width/height of a cell
     real(dp)              :: r_min($D) !< min coords. of box
     integer               :: coord_t   !< Coordinate type (e.g. Cartesian)
     class(*), pointer     :: ud=>null() !< User data (can be anything)
#if $D == 2
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
     integer                    :: coord_t    !< Type of coordinates
     real(dp)                   :: r_base($D) !< min. coords of box at index (1,1)
     real(dp)                   :: dr_base    !< cell spacing at lvl 1
     type(lvl_t), allocatable   :: lvls(:)    !< list storing the tree levels
     type(box$D_t), allocatable :: boxes(:)   !< list of all boxes
     logical                    :: ready = .false. !< Is tree ready for use?
  end type a$D_t

  !> Type specifying the location of a cell
  type a$D_loc_t
     integer :: id = -1         !< Id of the box that the cell is in
     integer :: ix($D) = -1     !< Index inside the box
  end type a$D_loc_t

  !> Type that contains the refinement changes in a level
  type ref_lvl_t
     integer, allocatable :: add(:) !< Id's of newly added boxes
     integer, allocatable :: rm(:) !< Id's of removed boxes
  end type ref_lvl_t

  !> Type that contains the refinement changes in a tree
  type ref_info_t
     integer :: n_add = 0                    !< Total number of added boxes
     integer :: n_rm = 0                     !< Total number removed boxes
     type(ref_lvl_t), allocatable :: lvls(:) !< Information per level
  end type ref_info_t

  abstract interface
     !> Function for setting refinement flags
     subroutine a$D_subr_ref(boxes, id, ref_flags)
       import
       type(box$D_t), intent(in) :: boxes(:) ! List of boxes
       integer, intent(in)       :: id  ! Id (index) of box
       integer, intent(inout)    :: ref_flags(:) ! Refinement flags
     end subroutine a$D_subr_ref

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
       integer, intent(in)         :: iv       !< Variable for which ghost cells are filled
     end subroutine a$D_subr_gc

     !> Subroutine for getting extra ghost cell data (> 1) near physical boundaries
     subroutine a$D_subr_egc(boxes, id, nb, iv, gc_data, nc)
       import
       type(box$D_t), intent(inout) :: boxes(:) !< Array with all boxes
       integer, intent(in)         :: id       !< Id of the box that needs to have ghost cells filled
       integer, intent(in)         :: nb       !< Neighbor direction
       integer, intent(in)         :: iv       !< Variable for which ghost cells are filled
       integer, intent(in)         :: nc       !< box%n_cell (this is purely for convenience)
#if $D == 2
       real(dp), intent(out)       :: gc_data(nc) !< The requested ghost cells
#elif $D == 3
       real(dp), intent(out)       :: gc_data(nc, nc) !< The requested ghost cells
#endif
     end subroutine a$D_subr_egc
  end interface

  private :: set_leaves_parents
  private :: child_that_contains
  private :: set_nbs_$Dd
  private :: find_nb_$Dd
  private :: get_free_ids
  private :: set_ref_info
  private :: consistent_ref_flags
  private :: remove_children
  private :: add_children
  private :: set_child_ids
  private :: sides_from_nb
  private :: sides2_from_nb

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

  !> Call procedure for each box in tree
  subroutine a$D_loop_box(tree, my_procedure, leaves_only)
    type(a$D_t), intent(inout)     :: tree
    procedure(a$D_subr)           :: my_procedure
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                       :: lvl, i, id

    leaves = .false.; if (present(leaves_only)) leaves = leaves_only
    if (.not. tree%ready) stop "a$D_loop_box: set_base has not been called"

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

  !> Get the location of the finest cell containing rr. If max_lvl is present,
  !> do not go to a finer level than max_lvl. If there is no box containing rr,
  !> return a location of -1
  pure function a$D_get_loc(tree, rr, max_lvl) result(loc)
    type(a$D_t), intent(in)       :: tree   !< Tree
    real(dp), intent(in)          :: rr($D) !< Coordinate
    integer, intent(in), optional :: max_lvl !< Maximum level of box
    type(a$D_loc_t)               :: loc    !< Location of cell

    integer :: i, id, i_ch, lvl_max

    lvl_max = tree%lvls_max
    if (present(max_lvl)) lvl_max = max_lvl

    ! Find lvl 1 box that includes rr
    do i = 1, size(tree%lvls(1)%ids)
       id = tree%lvls(1)%ids(i)
       if (a$D_r_inside(tree%boxes(id), rr)) exit
    end do

    ! If not inside any box, return
    if (i > size(tree%lvls(1)%ids)) then
       loc%id = -1
       loc%ix = -1
       return
    end if

    ! Jump into children for as long as possible
    do
       if (tree%boxes(id)%lvl < lvl_max .and. &
            a$D_has_children(tree%boxes(id))) then
          i_ch = child_that_contains(tree%boxes(id), rr)
          id = tree%boxes(id)%children(i_ch)
       else
          exit
       end if
    end do

    loc%id = id
    loc%ix = a$D_cc_ix(tree%boxes(id), rr)
  end function a$D_get_loc

  !> For a box with children that contains rr, find in which child rr lies
  pure function child_that_contains(box, rr) result(i_ch)
    type(box$D_t), intent(in) :: box    !< A box with children
    real(dp), intent(in)      :: rr($D) !< Location inside the box
    integer                   :: i_ch   !< Index of child containing rr
    real(dp)                  :: cntr($D)

    i_ch = 1
    cntr = box%r_min + box%dr * ishft(box%n_cell, -1)

    if (rr(1) > cntr(1)) i_ch = i_ch + 1
    if (rr(2) > cntr(2)) i_ch = i_ch + 2
#if $D==3
    if (rr(3) > cntr(3)) i_ch = i_ch + 4
#endif
  end function child_that_contains

  !> Call procedure for each box in tree, with argument rarg
  subroutine a$D_loop_box_arg(tree, my_procedure, rarg, leaves_only)
    type(a$D_t), intent(inout)     :: tree
    procedure(a$D_subr_arg)        :: my_procedure
    real(dp), intent(in)          :: rarg(:)
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                       :: lvl, i, id

    if (.not. tree%ready) stop "Tree not ready"
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

    if (.not. tree%ready) stop "Tree not ready"
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
  subroutine a$D_loop_boxes_arg(tree, my_procedure, rarg, leaves_only)
    type(a$D_t), intent(inout)    :: tree
    procedure(a$D_subr_boxes_arg) :: my_procedure
    real(dp), intent(in)         :: rarg(:)
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                      :: lvl, i, id

    if (.not. tree%ready) stop "Tree not ready"
    leaves = .false.; if (present(leaves_only)) leaves = leaves_only

    !$omp parallel private(lvl, i, id)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       if (leaves) then
          !$omp do
          do i = 1, size(tree%lvls(lvl)%leaves)
             id = tree%lvls(lvl)%leaves(i)
             call my_procedure(tree%boxes, id, rarg)
          end do
          !$omp end do
       else
          !$omp do
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             call my_procedure(tree%boxes, id, rarg)
          end do
          !$omp end do
       end if
    end do
    !$omp end parallel
  end subroutine a$D_loop_boxes_arg

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

  !> Get the offset of a box with respect to its parent (e.g. in 2d, there can
  !> be a child at offset 0,0, one at n_cell/2,0, one at 0,n_cell/2 etc.)
  function a$D_get_child_offset(box, nb) result(ix_offset)
    type(box$D_t), intent(in)           :: box   !< A child box
    integer, intent(in), optional      :: nb     !< Optional: get index on parent neighbor
    integer                            :: ix_offset($D)
    if (box%lvl > 1) then
       ix_offset = iand(box%ix-1, 1) * ishft(box%n_cell, -1) ! * n_cell / 2
       if (present(nb)) ix_offset = ix_offset - a$D_nb_dix(:, nb) * box%n_cell
    else                        ! In the subtree, parents are half the size
       ix_offset = 0
       if (present(nb)) ix_offset = ix_offset - &
            a$D_nb_dix(:, nb) * ishft(box%n_cell, -1) ! n_cell / 2
    endif
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

  !> Return n_cell at lvl. For all lvls >= 1, n_cell has the same value, but
  !> for lvls <= 0, n_cell changes.
  pure function a$D_n_cell(tree, lvl) result(n_cell)
    type(a$D_t), intent(in) :: tree
    integer, intent(in)    :: lvl
    integer                :: n_cell !< Output: n_cell at lvl

    if (lvl >= 1) then
       n_cell = tree%n_cell
    else
       n_cell = tree%n_cell / (2**(1-lvl))
    end if
  end function a$D_n_cell

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

  !> Get the index of the cell that includes point r
  pure function a$D_cc_ix(box, r) result(cc_ix)
    type(box$D_t), intent(in) :: box
    real(dp), intent(in)     :: r($D)
    integer                  :: cc_ix($D)
    cc_ix = ceiling((r - box%r_min) / box%dr)
  end function a$D_cc_ix

  !> Get the location of the cell center with index cc_ix
  pure function a$D_r_cc(box, cc_ix) result(r)
    type(box$D_t), intent(in) :: box
    integer, intent(in)      :: cc_ix($D)
    real(dp)                 :: r($D)
    r = box%r_min + (cc_ix-0.5_dp) * box%dr
  end function a$D_r_cc

  !> Get the location of "loc"
  pure function a$D_r_loc(tree, loc) result(r)
    type(a$D_t), intent(in)     :: tree
    type(a$D_loc_t), intent(in) :: loc
    real(dp)                   :: r($D)
    r = tree%boxes(loc%id)%r_min + &
         (loc%ix-0.5_dp) * tree%boxes(loc%id)%dr
  end function a$D_r_loc

#if $D == 2
  !> Get the radius of the cell center with index cc_ix
  pure function a$D_cyl_radius_cc(box, cc_ix) result(r)
    type(box$D_t), intent(in) :: box
    integer, intent(in)      :: cc_ix($D)
    real(dp)                 :: r
    r = box%r_min(a2_r_dim) + (cc_ix(a2_r_dim)-0.5_dp) * box%dr
  end function a$D_cyl_radius_cc
#endif

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

  !> Copy cc(..., iv_from) from box_in to cc(..., iv_to) on box_out
  subroutine a$D_box_copy_cc_to(box_from, iv_from, box_to, iv_to)
    type(box$D_t), intent(in)    :: box_from
    type(box$D_t), intent(inout) :: box_to
    integer, intent(in)         :: iv_from, iv_to
#if $D == 2
    box_to%cc(:,:, iv_to) = box_from%cc(:,:, iv_from)
#elif $D == 3
    box_to%cc(:,:,:, iv_to) = box_from%cc(:,:,:, iv_from)
#endif
  end subroutine a$D_box_copy_cc_to

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

  !> A general scalar reduction method
  !> TODO: test
  subroutine a$D_reduction(tree, box_func, reduction, init_val, out_val)
    type(a$D_t), intent(in) :: tree    !< Tree to do the reduction on
    real(dp), intent(in)   :: init_val !< Initial value for the reduction
    real(dp), intent(out)  :: out_val  !< Result of the reduction
    real(dp)               :: tmp, my_val
    integer                :: i, id, lvl

    interface
       real(dp) function box_func(box)
         import
         type(box$D_t), intent(in) :: box
       end function box_func

       real(dp) function reduction(a, b)
         import
         real(dp), intent(in) :: a, b
       end function reduction
    end interface

    if (.not. tree%ready) stop "Tree not ready"
    out_val = init_val
    my_val  = init_val

    !$omp parallel private(lvl, i, id, tmp) firstprivate(my_val)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          tmp = box_func(tree%boxes(id))
          my_val = reduction(tmp, my_val)
       end do
       !$omp end do
    end do

    !$omp critical
    out_val = reduction(my_val, out_val)
    !$omp end critical
    !$omp end parallel
  end subroutine a$D_reduction

  !> A general scalar reduction method, that returns the location of the
  !> minimum/maximum value found
  !> TODO: test
  subroutine a$D_reduction_loc(tree, box_subr, reduction, &
       init_val, out_val, out_loc)
    type(a$D_t), intent(in)      :: tree     !< Tree to do the reduction on
    real(dp), intent(in)        :: init_val !< Initial value for the reduction
    real(dp), intent(out)       :: out_val  !< Result of the reduction
    type(a$D_loc_t), intent(out) :: out_loc  !< Location
    real(dp)                    :: tmp, new_val, my_val
    integer                     :: i, id, lvl, tmp_ix($D)
    type(a$D_loc_t)             :: my_loc

    interface
       subroutine box_subr(box, val, ix)
         import
         type(box$D_t), intent(in) :: box
         real(dp), intent(out)    :: val
         integer, intent(out)     :: ix($D)
       end subroutine box_subr

       real(dp) function reduction(a, b)
         import
         real(dp), intent(in) :: a, b
       end function reduction
    end interface

    if (.not. tree%ready) stop "Tree not ready"
    out_val   = init_val
    my_val    = init_val
    my_loc%id = -1
    my_loc%ix = -1

    !$omp parallel private(lvl, i, id, tmp, tmp_ix, new_val) &
    !$omp firstprivate(my_val, my_loc)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call box_subr(tree%boxes(id), tmp, tmp_ix)
          new_val = reduction(tmp, my_val)
          if (abs(new_val - my_val) > 0) then
             my_loc%id = id
             my_loc%ix = tmp_ix
             my_val = tmp
          end if
       end do
       !$omp end do
    end do

    !$omp critical
    new_val = reduction(my_val, out_val)
    if (abs(new_val - out_val) > 0) then
       out_loc%id = my_loc%id
       out_loc%ix = my_loc%ix
       out_val = my_val
    end if
    !$omp end critical
    !$omp end parallel
  end subroutine a$D_reduction_loc

  !> Find maximum value of cc(..., iv). Only loop over leaves, and ghost cells
  !> are not used.
  subroutine a$D_tree_max_cc(tree, iv, cc_max)
    type(a$D_t), intent(in) :: tree
    integer, intent(in)    :: iv
    real(dp), intent(out)  :: cc_max
    real(dp)               :: tmp, my_max
    integer                :: i, id, lvl, nc

    if (.not. tree%ready) stop "Tree not ready"
    my_max = -huge(1.0_dp)

    !$omp parallel reduction(max: my_max) private(lvl, i, id, nc, tmp)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
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

  !> Find minimum value of cc(..., iv). Only loop over leaves, and ghost cells
  !> are not used.
  subroutine a$D_tree_min_cc(tree, iv, cc_min)
    type(a$D_t), intent(in) :: tree
    integer, intent(in)    :: iv
    real(dp), intent(out)  :: cc_min
    real(dp)               :: tmp, my_min
    integer                :: i, id, lvl, nc

    if (.not. tree%ready) stop "Tree not ready"
    my_min = huge(1.0_dp)

    !$omp parallel reduction(min: my_min) private(lvl, i, id, nc, tmp)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          nc = tree%boxes(id)%n_cell
#if $D == 2
          tmp = minval(tree%boxes(id)%cc(1:nc, 1:nc, iv))
#elif $D == 3
          tmp = minval(tree%boxes(id)%cc(1:nc, 1:nc, 1:nc, iv))
#endif
          if (tmp < my_min) my_min = tmp
       end do
       !$omp end do
    end do
    !$omp end parallel

    cc_min = my_min
  end subroutine a$D_tree_min_cc

  !> Find weighted sum of cc(..., iv). Only loop over leaves, and ghost cells
  !> are not used.
  subroutine a$D_tree_sum_cc(tree, iv, cc_sum)
    type(a$D_t), intent(in) :: tree
    integer, intent(in)    :: iv
    real(dp), intent(out)  :: cc_sum
    real(dp)               :: tmp, my_sum, fac
    integer                :: i, id, lvl, nc

    if (.not. tree%ready) stop "Tree not ready"
    my_sum = 0

    !$omp parallel reduction(+: my_sum) private(lvl, i, id, nc, tmp, fac)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       fac = a$D_lvl_dr(tree, lvl)**$D

       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          nc = tree%boxes(id)%n_cell
#if $D == 2
          if (tree%coord_t == a5_cyl) then
             tmp = sum_2pr_box(tree%boxes(id), iv)
          else
             tmp = sum(tree%boxes(id)%cc(1:nc, 1:nc, iv))
          end if
#elif $D == 3
          tmp = sum(tree%boxes(id)%cc(1:nc, 1:nc, 1:nc, iv))
#endif
          my_sum = my_sum + fac * tmp
       end do
       !$omp end do
    end do
    !$omp end parallel

    cc_sum = my_sum

#if $D == 2
  contains

    ! Sum of 2 * pi * r * values
    pure function sum_2pr_box(box, iv) result(res)
      type(box2_t), intent(in) :: box
      integer, intent(in)      :: iv
      real(dp), parameter      :: twopi = 2 * acos(-1.0_dp)
      real(dp)                 :: res
      integer                  :: i, j, nc

      res = 0
      nc  = box%n_cell

      do j = 1, nc
         do i = 1, nc
            res = res + box%cc(i, j, iv) * a2_cyl_radius_cc(box, [i, j])
         end do
      end do
      res = res * twopi
    end function sum_2pr_box
#endif
  end subroutine a$D_tree_sum_cc

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

    if (.not. tree%ready) stop "Tree not ready"
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       call a$D_boxes_copy_fc(tree%boxes, tree%lvls(lvl)%ids, iv_from, iv_to)
    end do
  end subroutine a$D_tree_copy_fc

  !> Zeroth-order prolongation to children.
  subroutine a$D_prolong0_from(boxes, id, iv)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Box whose children we will fill
    integer, intent(in)         :: iv        !< Variable that is filled
    integer                     :: i_c, c_id

    do i_c = 1, a$D_num_children
       c_id = boxes(id)%children(i_c)
       if (c_id == a5_no_box) cycle
       call a$D_prolong0_to(boxes, c_id, iv)
    end do
  end subroutine a$D_prolong0_from

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
    ix_offset = a$D_get_child_offset(boxes(id))

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

  !> Partial prolongation to the ghost cells of box id from parent
  subroutine a$D_sides_prolong0(boxes, id, nb, iv)
    type(box$D_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Id of child
    integer, intent(in)           :: iv       !< Variable to fill
    integer, intent(in)           :: nb       !< Neighbor to get data from
    integer                       :: nb_dim, lo($D), hi($D)

    nb_dim     = a$D_nb_dim(nb)
    lo(:)      = 1
    hi(:)      = boxes(id)%n_cell
    lo(nb_dim) = a$D_nb_hi01(nb) * (boxes(id)%n_cell+1)
    hi(nb_dim) = a$D_nb_hi01(nb) * (boxes(id)%n_cell+1)

    call a$D_prolong0_to(boxes, id, iv, lo, hi)
  end subroutine a$D_sides_prolong0

  !> Linear prolongation to children. We use 2-1-1 interpolation (2d) and
  !> 1-1-1-1 interpolation (3D), which do not require corner ghost cells.
  subroutine a$D_prolong1_from(boxes, id, iv, add)
    type(box$D_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Box whose children we will fill
    integer, intent(in)           :: iv       !< Variable that is filled
    logical, intent(in), optional :: add      !< Add to old values
    integer                       :: i_c, c_id

    do i_c = 1, a$D_num_children
       c_id = boxes(id)%children(i_c)
       if (c_id == a5_no_box) cycle
       call a$D_prolong1_to(boxes, c_id, iv, add)
    end do
  end subroutine a$D_prolong1_from

  !> Prolongation to a child (from parent) using linear interpolation. We use
  !> 2-1-1 interpolation (2D) and 1-1-1-1 interpolation (3D) which do not need
  !> corner ghost cells.
  subroutine a$D_prolong1_to(boxes, id, iv, add)
    type(box$D_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Id of child
    integer, intent(in)           :: iv       !< Variable to fill
    logical, intent(in), optional :: add      !< Add to old values
    integer                       :: hnc, nc, p_id, ix_offset($D)
    integer                       :: i, j, i_c, i_f, j_c, j_f
    real(dp)                      :: f0, flx, fhx, fly, fhy
    logical                       :: add_to
#if $D == 3
    real(dp)                      :: flz, fhz
    integer                       :: k, k_c, k_f
#endif

    nc        = boxes(id)%n_cell
    hnc       = ishft(boxes(id)%n_cell, -1)
    p_id      = boxes(id)%parent
    ix_offset = a$D_get_child_offset(boxes(id))
    add_to    = .false.; if (present(add)) add_to = add

    if (.not. add_to) then
#if $D == 2
       boxes(id)%cc(1:nc, 1:nc, iv) = 0
#elif $D == 3
       boxes(id)%cc(1:nc, 1:nc, 1:nc, iv) = 0
#endif
    end if

#if $D == 2
    do j = 1, hnc
       j_c = j + ix_offset(2)
       j_f = 2 * j - 1
       do i = 1, hnc
          i_c = i + ix_offset(1)
          i_f = 2 * i - 1

          f0 = 0.5_dp * boxes(p_id)%cc(i_c, j_c, iv)
          flx = 0.25_dp * boxes(p_id)%cc(i_c-1, j_c, iv)
          fhx = 0.25_dp * boxes(p_id)%cc(i_c+1, j_c, iv)
          fly = 0.25_dp * boxes(p_id)%cc(i_c, j_c-1, iv)
          fhy = 0.25_dp * boxes(p_id)%cc(i_c, j_c+1, iv)

          boxes(id)%cc(i_f,   j_f,   iv) = f0 + flx + fly &
               + boxes(id)%cc(i_f,   j_f,   iv)
          boxes(id)%cc(i_f+1, j_f,   iv) = f0 + fhx + fly &
               + boxes(id)%cc(i_f+1, j_f,   iv)
          boxes(id)%cc(i_f,   j_f+1, iv) = f0 + flx + fhy &
               + boxes(id)%cc(i_f,   j_f+1, iv)
          boxes(id)%cc(i_f+1, j_f+1, iv) = f0 + fhx + fhy &
               + boxes(id)%cc(i_f+1, j_f+1, iv)
       end do
    end do
#elif $D == 3
    do k = 1, hnc
       k_c = k + ix_offset(3)
       k_f = 2 * k - 1
       do j = 1, hnc
          j_c = j + ix_offset(2)
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = i + ix_offset(1)
             i_f = 2 * i - 1

             f0  = 0.25_dp * boxes(p_id)%cc(i_c,   j_c,   k_c,   iv)
             flx = 0.25_dp * boxes(p_id)%cc(i_c-1, j_c,   k_c,   iv)
             fhx = 0.25_dp * boxes(p_id)%cc(i_c+1, j_c,   k_c,   iv)
             fly = 0.25_dp * boxes(p_id)%cc(i_c,   j_c-1, k_c,   iv)
             fhy = 0.25_dp * boxes(p_id)%cc(i_c,   j_c+1, k_c,   iv)
             flz = 0.25_dp * boxes(p_id)%cc(i_c,   j_c,   k_c-1, iv)
             fhz = 0.25_dp * boxes(p_id)%cc(i_c,   j_c,   k_c+1, iv)

             boxes(id)%cc(i_f,   j_f,   k_f,   iv) = f0 + flx + &
                  fly + flz + boxes(id)%cc(i_f,   j_f,   k_f,   iv)
             boxes(id)%cc(i_f+1, j_f,   k_f,   iv) = f0 + fhx + &
                  fly + flz + boxes(id)%cc(i_f+1, j_f,   k_f,   iv)
             boxes(id)%cc(i_f,   j_f+1, k_f,   iv) = f0 + flx + &
                  fhy + flz + boxes(id)%cc(i_f,   j_f+1, k_f,   iv)
             boxes(id)%cc(i_f+1, j_f+1, k_f,   iv) = f0 + fhx + &
                  fhy + flz + boxes(id)%cc(i_f+1, j_f+1, k_f,   iv)
             boxes(id)%cc(i_f,   j_f,   k_f+1, iv) = f0 + flx + &
                  fly + fhz + boxes(id)%cc(i_f,   j_f,   k_f+1, iv)
             boxes(id)%cc(i_f+1, j_f,   k_f+1, iv) = f0 + fhx + &
                  fly + fhz + boxes(id)%cc(i_f+1, j_f,   k_f+1, iv)
             boxes(id)%cc(i_f,   j_f+1, k_f+1, iv) = f0 + flx + &
                  fhy + fhz + boxes(id)%cc(i_f,   j_f+1, k_f+1, iv)
             boxes(id)%cc(i_f+1, j_f+1, k_f+1, iv) = f0 + fhx + &
                  fhy + fhz + boxes(id)%cc(i_f+1, j_f+1, k_f+1, iv)
          end do
       end do
    end do
#endif
  end subroutine a$D_prolong1_to

  !> Quadratic prolongation to children. We use stencils that do not require
  !> corner ghost cells.
  subroutine a$D_prolong2_from(boxes, id, iv)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Box whose children we will fill
    integer, intent(in)         :: iv        !< Variable that is filled
    integer                     :: i_c, c_id

    do i_c = 1, a$D_num_children
       c_id = boxes(id)%children(i_c)
       if (c_id == a5_no_box) cycle
       call a$D_prolong2_to(boxes, c_id, iv)
    end do
  end subroutine a$D_prolong2_from

  !> Prolongation to a child (from parent) using quadratic interpolation. We use
  !> 5 / 7 point stencils which do not need corner ghost cells.
  !> @TODO 3D version
  subroutine a$D_prolong2_to(boxes, id, iv)
    type(box$D_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)          :: id       !< Id of child
    integer, intent(in)          :: iv       !< Variable to fill
    integer                      :: hnc, p_id, ix_offset($D)
    integer                      :: i, j
    integer                      :: i_c, i_f, j_c, j_f
    real(dp)                     :: f0, fx, fy, fxx, fyy, f2
#if $D == 3
    real(dp)                     :: fz, fzz
    integer                      :: k, k_c, k_f
#endif

    hnc       = ishft(boxes(id)%n_cell, -1)
    p_id      = boxes(id)%parent
    ix_offset = a$D_get_child_offset(boxes(id))

#if $D == 2
    do j = 1, hnc
       j_c = j + ix_offset(2)
       j_f = 2 * j - 1
       do i = 1, hnc
          i_c = i + ix_offset(1)
          i_f = 2 * i - 1

          f0 = boxes(p_id)%cc(i_c, j_c, iv)
          fx = 0.125_dp * (boxes(p_id)%cc(i_c+1, j_c, iv) - &
               boxes(p_id)%cc(i_c-1, j_c, iv))
          fy = 0.125_dp * (boxes(p_id)%cc(i_c, j_c+1, iv) - &
               boxes(p_id)%cc(i_c, j_c-1, iv))
          fxx = 0.03125_dp * (boxes(p_id)%cc(i_c-1, j_c, iv) - &
               2 * f0 + boxes(p_id)%cc(i_c+1, j_c, iv))
          fyy = 0.03125_dp * (boxes(p_id)%cc(i_c, j_c-1, iv) - &
               2 * f0 + boxes(p_id)%cc(i_c, j_c+1, iv))
          f2 = fxx + fyy

          boxes(id)%cc(i_f,   j_f,   iv) = f0 - fx - fy + f2
          boxes(id)%cc(i_f+1, j_f,   iv) = f0 + fx - fy + f2
          boxes(id)%cc(i_f,   j_f+1, iv) = f0 - fx + fy + f2
          boxes(id)%cc(i_f+1, j_f+1, iv) = f0 + fx + fy + f2
       end do
    end do
#elif $D == 3
    do k = 1, hnc
       k_c = k + ix_offset(3)
       k_f = 2 * k - 1
       do j = 1, hnc
          j_c = j + ix_offset(2)
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = i + ix_offset(1)
             i_f = 2 * i - 1

             f0 = boxes(p_id)%cc(i_c, j_c, k_c, iv)
             fx = 0.125_dp * (boxes(p_id)%cc(i_c+1, j_c, k_c, iv) - &
                  boxes(p_id)%cc(i_c-1, j_c, k_c, iv))
             fy = 0.125_dp * (boxes(p_id)%cc(i_c, j_c+1, k_c, iv) - &
                  boxes(p_id)%cc(i_c, j_c-1, k_c, iv))
             fz = 0.125_dp * (boxes(p_id)%cc(i_c, j_c, k_c+1, iv) - &
                  boxes(p_id)%cc(i_c, j_c, k_c-1, iv))
             fxx = 0.03125_dp * (boxes(p_id)%cc(i_c-1, j_c, k_c, iv) - &
                  2 * f0 + boxes(p_id)%cc(i_c+1, j_c, k_c, iv))
             fyy = 0.03125_dp * (boxes(p_id)%cc(i_c, j_c-1, k_c, iv) - &
                  2 * f0 + boxes(p_id)%cc(i_c, j_c+1, k_c, iv))
             fzz = 0.03125_dp * (boxes(p_id)%cc(i_c, j_c, k_c-1, iv) - &
                  2 * f0 + boxes(p_id)%cc(i_c, j_c, k_c+1, iv))
             f2 = fxx + fyy + fzz

             boxes(id)%cc(i_f,   j_f,   k_f,   iv) = f0 - fx - fy - fz + f2
             boxes(id)%cc(i_f+1, j_f,   k_f,   iv) = f0 + fx - fy - fz + f2
             boxes(id)%cc(i_f,   j_f+1, k_f,   iv) = f0 - fx + fy - fz + f2
             boxes(id)%cc(i_f+1, j_f+1, k_f,   iv) = f0 + fx + fy - fz + f2
             boxes(id)%cc(i_f,   j_f,   k_f+1, iv) = f0 - fx - fy + fz + f2
             boxes(id)%cc(i_f+1, j_f,   k_f+1, iv) = f0 + fx - fy + fz + f2
             boxes(id)%cc(i_f,   j_f+1, k_f+1, iv) = f0 - fx + fy + fz + f2
             boxes(id)%cc(i_f+1, j_f+1, k_f+1, iv) = f0 + fx + fy + fz + f2
          end do
       end do
    end do
#endif
  end subroutine a$D_prolong2_to

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

    if (.not. tree%ready) stop "Tree not ready"
    do lvl = tree%max_lvl-1, lbound(tree%lvls, 1), -1
       call a$D_restrict_to_boxes(tree%boxes, tree%lvls(lvl)%parents, iv, i_to)
    end do
  end subroutine a$D_restrict_tree

  !> Restriction of child box (box_c) to its parent (box_p)
  subroutine a$D_restrict_box(box_c, box_p, iv, i_to)
    type(box$D_t), intent(in)      :: box_c         !< Child box to restrict
    type(box$D_t), intent(inout)   :: box_p         !< Parent box to restrict to
    integer, intent(in)           :: iv            !< Variable to restrict
    integer, intent(in), optional :: i_to         !< Destination (if /= iv)
    integer                       :: i, j, i_f, j_f, i_c, j_c, i_dest
    integer                       :: hnc, ix_offset($D)
#if $D == 2
    real(dp)                      :: r, dr16, rfac
#elif $D == 3
    integer                       :: k, k_f, k_c
#endif

    hnc       = ishft(box_c%n_cell, -1) ! n_cell / 2
    ix_offset = a$D_get_child_offset(box_c)

    if (present(i_to)) then
       i_dest = i_to
    else
       i_dest = iv
    end if

#if $D == 2
    if (box_p%coord_t == a5_cyl) then
       dr16 = 0.0625_dp * box_p%dr   ! (dr / 4) / 4

       do j = 1, hnc
          j_c = ix_offset(2) + j
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = ix_offset(1) + i
             i_f = 2 * i - 1

             ! The weight of cells is proportional to their radius.
             r = a2_cyl_radius_cc(box_p, [i, j])
             rfac = dr16 / r

             box_p%cc(i_c, j_c, i_dest) = &
                  (0.25_dp - rfac) * sum(box_c%cc(i_f, j_f:j_f+1, iv)) + &
                  (0.25_dp + rfac) * sum(box_c%cc(i_f+1, j_f:j_f+1, iv))
          end do
       end do
    else
       do j = 1, hnc
          j_c = ix_offset(2) + j
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = ix_offset(1) + i
             i_f = 2 * i - 1
             box_p%cc(i_c, j_c, i_dest) = 0.25_dp * &
                  sum(box_c%cc(i_f:i_f+1, j_f:j_f+1, iv))
          end do
       end do
    endif
#elif $D == 3
    do k = 1, hnc
       k_c = ix_offset(3) + k
       k_f = 2 * k - 1
       do j = 1, hnc
          j_c = ix_offset(2) + j
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = ix_offset(1) + i
             i_f = 2 * i - 1
             box_p%cc(i_c, j_c, k_c, i_dest) = 0.125_dp * &
                  sum(box_c%cc(i_f:i_f+1, j_f:j_f+1, k_f:k_f+1, iv))
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

    if (.not. tree%ready) stop "Tree not ready"
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
    procedure(a$D_subr_gc)      :: subr_bc !< Procedure called at physical boundaries
    integer                     :: nb, nb_id

    do nb = 1, a$D_num_neighbors
       nb_id = boxes(id)%neighbors(nb)
       if (nb_id > a5_no_box) then
          call sides_from_nb(boxes(id), boxes(nb_id), nb, iv)
       else if (nb_id == a5_no_box) then
          call subr_rb(boxes, id, nb, iv)
       else
          call subr_bc(boxes, id, nb, iv)
       end if
    end do
  end subroutine a$D_gc_box_sides

  !> Get a second layer of ghost cell data (the 'normal' routines give just one
  !> layer of ghost cells). Use subr_rb > on refinement boundaries and subr_bc
  !> on physical boundaries.
  subroutine a$D_gc2_box_sides(boxes, id, iv, subr_rb, subr_bc, gc_data, nc)
    type(box$D_t), intent(inout) :: boxes(:)        !< List of all the boxes
    integer, intent(in)          :: id              !< Id of box for which we set ghost cells
    integer, intent(in)          :: iv              !< Variable for which ghost cells are set
    procedure(a$D_subr_egc)      :: subr_rb         !< Procedure called at refinement boundaries
    procedure(a$D_subr_egc)      :: subr_bc         !< Procedure called at physical boundaries
    integer, intent(in)          :: nc              !< box%n_cell
#if $D   == 2
    real(dp), intent(out)        :: gc_data(nc, 2*$D)     !< The requested ghost cells
#elif $D == 3
    real(dp), intent(out)        :: gc_data(nc, nc, 2*$D) !< The requested ghost cells
#endif
    integer                      :: nb, nb_id

    do nb = 1, a$D_num_neighbors
       nb_id = boxes(id)%neighbors(nb)
       if (nb_id > a5_no_box) then
#if $D == 2
          call sides2_from_nb(boxes(nb_id), nb, iv, gc_data(:, nb), nc)
#elif $D == 3
          call sides2_from_nb(boxes(nb_id), nb, iv, gc_data(:, :, nb), nc)
#endif
       else if (nb_id == a5_no_box) then
#if $D == 2
          call subr_rb(boxes, id, nb, iv, gc_data(:, nb), nc)
#elif $D == 3
          call subr_rb(boxes, id, nb, iv, gc_data(:, :, nb), nc)
#endif
       else
#if $D == 2
          call subr_bc(boxes, id, nb, iv, gc_data(:, nb), nc)
#elif $D == 3
          call subr_bc(boxes, id, nb, iv, gc_data(:, :, nb), nc)
#endif
       end if
    end do
  end subroutine a$D_gc2_box_sides

  !> Interpolate between fine points and coarse neighbors to fill ghost cells
  !> near refinement boundaries
  subroutine a$D_sides_interp(boxes, id, nb, iv)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc, ix, ix_c, ix_f, i, j
    integer                     :: i_c1, i_c2, j_c1, j_c2, p_nb_id
    integer                     :: p_id, ix_offset($D)
    real(dp), parameter         :: sixth=1/6.0_dp, third=1/3.0_dp
#if $D == 3
    integer                     :: k_c1, k_c2, k
#endif

    nc        = boxes(id)%n_cell
    p_id      = boxes(id)%parent
    p_nb_id   = boxes(p_id)%neighbors(nb)
    ix_offset = a$D_get_child_offset(boxes(id), nb)

    if (a$D_nb_low(nb)) then
       ix = 0
       ix_f = 1
       ix_c = nc
    else
       ix = nc+1
       ix_f = nc
       ix_c = 1
    end if

    select case (a$D_nb_dim(nb))
#if $D == 2
    case (1)
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          boxes(id)%cc(ix, j, iv) = &
               0.5_dp * boxes(p_nb_id)%cc(ix_c, j_c1, iv) + &
               sixth * boxes(p_nb_id)%cc(ix_c, j_c2, iv) + &
               third * boxes(id)%cc(ix_f, j, iv)
       end do
    case (2)
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
          boxes(id)%cc(i, ix, iv) = &
               0.5_dp * boxes(p_nb_id)%cc(i_c1, ix_c, iv) + &
               sixth * boxes(p_nb_id)%cc(i_c2, ix_c, iv) + &
               third * boxes(id)%cc(i, ix_f, iv)
       end do
#elif $D==3
    case (1)
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)          ! even: +1, odd: -1
          do j = 1, nc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
             boxes(id)%cc(ix, j, k, iv) = &
                  third * boxes(p_nb_id)%cc(ix_c, j_c1, k_c1, iv) + &
                  sixth * boxes(p_nb_id)%cc(ix_c, j_c2, k_c1, iv) + &
                  sixth * boxes(p_nb_id)%cc(ix_c, j_c1, k_c2, iv) + &
                  third * boxes(id)%cc(ix_f, j, k, iv)
          end do
       end do
    case (2)
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             boxes(id)%cc(i, ix, k, iv) = &
                  third * boxes(p_nb_id)%cc(i_c1, ix_c, k_c1, iv) + &
                  sixth * boxes(p_nb_id)%cc(i_c2, ix_c, k_c1, iv) + &
                  sixth * boxes(p_nb_id)%cc(i_c1, ix_c, k_c2, iv) + &
                  third * boxes(id)%cc(i, ix_f, k, iv)
          end do
       end do
    case (3)
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             boxes(id)%cc(i, j, ix, iv) = &
                  third * boxes(p_nb_id)%cc(i_c1, j_c1, ix_c, iv) + &
                  sixth * boxes(p_nb_id)%cc(i_c1, j_c2, ix_c, iv) + &
                  sixth * boxes(p_nb_id)%cc(i_c2, j_c1, ix_c, iv) + &
                  third * boxes(id)%cc(i, j, ix_f, iv)
          end do
       end do
#endif
    end select

  end subroutine a$D_sides_interp

  !> Interpolate between fine points and coarse neighbors to fill ghost cells
  !> near refinement boundaries. The ghost values are less than twice the coarse
  !> values.
  subroutine a$D_sides_interp_lim(boxes, id, nb, iv)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc, ix, ix_c, ix_f, i, j
    integer                     :: i_c1, i_c2, j_c1, j_c2, p_nb_id
    integer                     :: p_id, ix_offset($D)
    real(dp)                    :: c1, c2
    real(dp), parameter         :: sixth=1/6.0_dp, third=1/3.0_dp
#if $D == 3
    integer                     :: k_c1, k_c2, k
    real(dp)                    :: c3
#endif

    nc        = boxes(id)%n_cell
    p_id      = boxes(id)%parent
    p_nb_id   = boxes(p_id)%neighbors(nb)
    ix_offset = a$D_get_child_offset(boxes(id), nb)

    if (a$D_nb_low(nb)) then
       ix = 0
       ix_f = 1
       ix_c = nc
    else
       ix = nc+1
       ix_f = nc
       ix_c = 1
    end if

    select case (a$D_nb_dim(nb))
#if $D == 2
    case (1)
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          c1 = boxes(p_nb_id)%cc(ix_c, j_c1, iv)
          c2 = boxes(p_nb_id)%cc(ix_c, j_c2, iv)
          boxes(id)%cc(ix, j, iv) = 0.5_dp * c1 + sixth * c2 + &
               third * boxes(id)%cc(ix_f, j, iv)
          if (boxes(id)%cc(ix, j, iv) > 2 * c1) boxes(id)%cc(ix, j, iv) = 2 * c1
       end do
    case (2)
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
          c1 = boxes(p_nb_id)%cc(i_c1, ix_c, iv)
          c2 = boxes(p_nb_id)%cc(i_c2, ix_c, iv)
          boxes(id)%cc(i, ix, iv) = 0.5_dp * c1 + sixth * c2 + &
               third * boxes(id)%cc(i, ix_f, iv)
          if (boxes(id)%cc(i, ix, iv) > 2 * c1) boxes(id)%cc(i, ix, iv) = 2 * c1
       end do
#elif $D==3
    case (1)
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)          ! even: +1, odd: -1
          do j = 1, nc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
             c1 = boxes(p_nb_id)%cc(ix_c, j_c1, k_c1, iv)
             c2 = boxes(p_nb_id)%cc(ix_c, j_c2, k_c1, iv)
             c3 = boxes(p_nb_id)%cc(ix_c, j_c1, k_c2, iv)
             boxes(id)%cc(ix, j, k, iv) = third * c1 + sixth * c2 + &
                  sixth * c3 + third * boxes(id)%cc(ix_f, j, k, iv)
             if (boxes(id)%cc(ix, j, k, iv) > 2 * c1) &
                  boxes(id)%cc(ix, j, k, iv) = 2 * c1
          end do
       end do
    case (2)
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             c1 = boxes(p_nb_id)%cc(i_c1, ix_c, k_c1, iv)
             c2 = boxes(p_nb_id)%cc(i_c2, ix_c, k_c1, iv)
             c3 = boxes(p_nb_id)%cc(i_c1, ix_c, k_c2, iv)
             boxes(id)%cc(i, ix, k, iv) = third * c1 + sixth * c2 + &
                  sixth * c3 + third * boxes(id)%cc(i, ix_f, k, iv)
             if (boxes(id)%cc(i, ix, k, iv) > 2 * c1) &
                  boxes(id)%cc(i, ix, k, iv) = 2 * c1
          end do
       end do
    case (3)
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             c1 = boxes(p_nb_id)%cc(i_c1, j_c1, ix_c, iv)
             c2 = boxes(p_nb_id)%cc(i_c1, j_c2, ix_c, iv)
             c3 = boxes(p_nb_id)%cc(i_c2, j_c1, ix_c, iv)
             boxes(id)%cc(i, j, ix, iv) = third * c1 + sixth * c2 + &
                  sixth * c3 + third * boxes(id)%cc(i, j, ix_f, iv)
             if (boxes(id)%cc(i, j, ix, iv) > 2 * c1) &
                  boxes(id)%cc(i, j, ix, iv) = 2 * c1
          end do
       end do
#endif
    end select

  end subroutine a$D_sides_interp_lim

  ! This fills ghost cells near physical boundaries using Neumann zero
  subroutine a$D_bc_neumann(boxes, id, nb, iv)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    integer                     :: nc

    nc = boxes(id)%n_cell

    select case (nb)
#if $D == 2
    case (a2_nb_lx)
       boxes(id)%cc(0, 1:nc, iv) = boxes(id)%cc(1, 1:nc, iv)
    case (a2_nb_hx)
       boxes(id)%cc(nc+1, 1:nc, iv) = boxes(id)%cc(nc, 1:nc, iv)
    case (a2_nb_ly)
       boxes(id)%cc(1:nc, 0, iv) = boxes(id)%cc(1:nc, 1, iv)
    case (a2_nb_hy)
       boxes(id)%cc(1:nc, nc+1, iv) = boxes(id)%cc(1:nc, nc, iv)
#elif $D == 3
    case (a3_nb_lx)
       boxes(id)%cc(0, 1:nc, 1:nc, iv) = boxes(id)%cc(1, 1:nc, 1:nc, iv)
    case (a3_nb_hx)
       boxes(id)%cc(nc+1, 1:nc, 1:nc, iv) = boxes(id)%cc(nc, 1:nc, 1:nc, iv)
    case (a3_nb_ly)
       boxes(id)%cc(1:nc, 0, 1:nc, iv) = boxes(id)%cc(1:nc, 1, 1:nc, iv)
    case (a3_nb_hy)
       boxes(id)%cc(1:nc, nc+1, 1:nc, iv) = boxes(id)%cc(1:nc, nc, 1:nc, iv)
    case (a3_nb_lz)
       boxes(id)%cc(1:nc, 1:nc, 0, iv) = boxes(id)%cc(1:nc, 1:nc, 1, iv)
    case (a3_nb_hz)
       boxes(id)%cc(1:nc, 1:nc, nc+1, iv) = boxes(id)%cc(1:nc, 1:nc, nc, iv)
#endif
    end select
  end subroutine a$D_bc_neumann

    !> This fills ghost cells near physical boundaries using Dirichlet zero
  subroutine a$D_bc_dirichlet(boxes, id, nb, iv)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    integer                     :: nc

    nc = boxes(id)%n_cell

    select case (nb)
#if $D == 2
    case (a2_nb_lx)
       boxes(id)%cc(0, 1:nc, iv) = -boxes(id)%cc(1, 1:nc, iv)
    case (a2_nb_hx)
       boxes(id)%cc(nc+1, 1:nc, iv) = -boxes(id)%cc(nc, 1:nc, iv)
    case (a2_nb_ly)
       boxes(id)%cc(1:nc, 0, iv) = -boxes(id)%cc(1:nc, 1, iv)
    case (a2_nb_hy)
       boxes(id)%cc(1:nc, nc+1, iv) = -boxes(id)%cc(1:nc, nc, iv)
#elif $D == 3
    case (a3_nb_lx)
       boxes(id)%cc(0, 1:nc, 1:nc, iv) = -boxes(id)%cc(1, 1:nc, 1:nc, iv)
    case (a3_nb_hx)
       boxes(id)%cc(nc+1, 1:nc, 1:nc, iv) = -boxes(id)%cc(nc, 1:nc, 1:nc, iv)
    case (a3_nb_ly)
       boxes(id)%cc(1:nc, 0, 1:nc, iv) = -boxes(id)%cc(1:nc, 1, 1:nc, iv)
    case (a3_nb_hy)
       boxes(id)%cc(1:nc, nc+1, 1:nc, iv) = -boxes(id)%cc(1:nc, nc, 1:nc, iv)
    case (a3_nb_lz)
       boxes(id)%cc(1:nc, 1:nc, 0, iv) = -boxes(id)%cc(1:nc, 1:nc, 1, iv)
    case (a3_nb_hz)
       boxes(id)%cc(1:nc, 1:nc, nc+1, iv) = -boxes(id)%cc(1:nc, 1:nc, nc, iv)
#endif
    end select
  end subroutine a$D_bc_dirichlet

  !> Fill values on the side of a box from a neighbor nb
  subroutine sides_from_nb(box, box_nb, nb, iv)
    type(box$D_t), intent(inout) :: box    !< Box on which to fill ghost cells
    type(box$D_t), intent(in)    :: box_nb !< Neighbouring box
    integer, intent(in)         :: nb        !< Ghost cell / neighbor direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc

    nc = box%n_cell

    select case (nb)
#if $D == 2
    case (a2_nb_lx)
       box%cc(0, 1:nc, iv)    = box_nb%cc(nc, 1:nc, iv)
    case (a2_nb_hx)
       box%cc(nc+1, 1:nc, iv) = box_nb%cc(1, 1:nc, iv)
    case (a2_nb_ly)
       box%cc(1:nc, 0, iv)    = box_nb%cc(1:nc, nc, iv)
    case (a2_nb_hy)
       box%cc(1:nc, nc+1, iv) = box_nb%cc(1:nc, 1, iv)
#elif $D == 3
    case (a3_nb_lx)
       box%cc(0, 1:nc, 1:nc, iv)    = box_nb%cc(nc, 1:nc, 1:nc, iv)
    case (a3_nb_hx)
       box%cc(nc+1, 1:nc, 1:nc, iv) = box_nb%cc(1, 1:nc, 1:nc, iv)
    case (a3_nb_ly)
       box%cc(1:nc, 0, 1:nc, iv)    = box_nb%cc(1:nc, nc, 1:nc, iv)
    case (a3_nb_hy)
       box%cc(1:nc, nc+1, 1:nc, iv) = box_nb%cc(1:nc, 1, 1:nc, iv)
    case (a3_nb_lz)
       box%cc(1:nc, 1:nc, 0, iv)    = box_nb%cc(1:nc, 1:nc, nc, iv)
    case (a3_nb_hz)
       box%cc(1:nc, 1:nc, nc+1, iv) = box_nb%cc(1:nc, 1:nc, 1, iv)
#endif
    end select
  end subroutine sides_from_nb

  !> Fill values on the side of a box from a neighbor nb
  subroutine sides2_from_nb(box_nb, nb, iv, gc_side, nc)
    type(box$D_t), intent(in) :: box_nb !< Neighbouring box
    integer, intent(in)       :: nb     !< Ghost cell / neighbor direction
    integer, intent(in)       :: iv     !< Ghost cell variable
    integer, intent(in)       :: nc
#if $D == 2
    real(dp), intent(out)     :: gc_side(nc)
#elif $D == 3
    real(dp), intent(out)     :: gc_side(nc, nc)
#endif

    select case (nb)
#if $D == 2
    case (a2_nb_lx)
       gc_side = box_nb%cc(nc-1, 1:nc, iv)
    case (a2_nb_hx)
       gc_side = box_nb%cc(2, 1:nc, iv)
    case (a2_nb_ly)
       gc_side = box_nb%cc(1:nc, nc-1, iv)
    case (a2_nb_hy)
       gc_side = box_nb%cc(1:nc, 2, iv)
#elif $D == 3
    case (a3_nb_lx)
       gc_side = box_nb%cc(nc-1, 1:nc, 1:nc, iv)
    case (a3_nb_hx)
       gc_side = box_nb%cc(2, 1:nc, 1:nc, iv)
    case (a3_nb_ly)
       gc_side = box_nb%cc(1:nc, nc-1, 1:nc, iv)
    case (a3_nb_hy)
       gc_side = box_nb%cc(1:nc, 2, 1:nc, iv)
    case (a3_nb_lz)
       gc_side = box_nb%cc(1:nc, 1:nc, nc-1, iv)
    case (a3_nb_hz)
       gc_side = box_nb%cc(1:nc, 1:nc, 2, iv)
#endif
    end select
  end subroutine sides2_from_nb

  !> Linear interpolation (using data from neighbor) to fill ghost cells
  subroutine a$D_sides2_prolong1(boxes, id, nb, iv, gc_side, nc)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer, intent(in)         :: nc
#if $D == 2
    real(dp), intent(out)       :: gc_side(nc)
#elif $D == 3
    real(dp), intent(out)       :: gc_side(nc, nc)
#endif
    integer                     :: ix, i, j
    integer                     :: i_c1, i_c2, j_c1, j_c2, p_nb_id
    integer                     :: p_id, ix_offset($D)
#if $D == 3
    integer                     :: k, k_c1, k_c2
#endif

    p_id      = boxes(id)%parent
    p_nb_id   = boxes(p_id)%neighbors(nb)
    ix_offset = a$D_get_child_offset(boxes(id), nb)
    ix        = a$D_nb_hi01(nb) * (nc+3) - 1 ! -1 or nc+2

    select case (a$D_nb_dim(nb))
#if $D == 2
    case (1)
       i_c1 = ix_offset(1) + ishft(ix+1, -1) ! (ix+1)/2
       i_c2 = i_c1 + 1 - 2 * iand(ix, 1)     ! even: +1, odd: -1
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          gc_side(j) = &
               0.5_dp * boxes(p_nb_id)%cc(i_c1, j_c1, iv) + &
               0.25_dp * boxes(p_nb_id)%cc(i_c2, j_c1, iv) + &
               0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c2, iv)
       end do
    case (2)
       j_c1 = ix_offset(2) + ishft(ix+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(ix, 1)     ! even: +1, odd: -1
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
          gc_side(i) = &
               0.5_dp * boxes(p_nb_id)%cc(i_c1, j_c1, iv) + &
               0.25_dp * boxes(p_nb_id)%cc(i_c2, j_c1, iv) + &
               0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c2, iv)
       end do
#elif $D==3
    case (1)
       i_c1 = ix_offset(1) + ishft(ix+1, -1) ! (ix+1)/2
       i_c2 = i_c1 + 1 - 2 * iand(ix, 1)     ! even: +1, odd: -1
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
          do j = 1, nc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
             gc_side(j, k) = &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c2, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c2, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c2, iv)
          end do
       end do
    case (2)
       j_c1 = ix_offset(2) + ishft(ix+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(ix, 1)     ! even: +1, odd: -1
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
             gc_side(i, k) = &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c2, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c2, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c2, iv)
          end do
       end do
    case (3)
       k_c1 = ix_offset(3) + ishft(ix+1, -1) ! (k+1)/2
       k_c2 = k_c1 + 1 - 2 * iand(ix, 1)     ! even: +1, odd: -1
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
             gc_side(i, j) = &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c2, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c2, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c2, iv)
          end do
       end do
#endif
    end select

  end subroutine a$D_sides2_prolong1

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
                   call a$D_flux_from_children(tree%boxes, id, nb, f_ixs)
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
  subroutine a$D_flux_from_children(boxes, id, nb, f_ixs)
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
  end subroutine a$D_flux_from_children

  !> Write the cell centered data of a tree to a vtk unstructured file. Only the
  !> leaves of the tree are used
  subroutine a$D_write_vtk(tree, filename, cc_names, n_cycle, time, ixs_cc, &
       fc_names, ixs_fc, dir)
    use m_vtk
    type(a$D_t), intent(in)       :: tree        !< Tree to write out
    character(len=*), intent(in)  :: filename    !< Filename for the vtk file
    character(len=*), intent(in)  :: cc_names(:) !< Names of the cell-centered variables
    integer, intent(in)           :: n_cycle     !< Cycle-number for vtk file (counter)
    real(dp), intent(in)          :: time        !< Time for output file
    integer, intent(in), optional :: ixs_cc(:)   !< Oncly include these cell variables
    character(len=*), optional, intent(in) :: dir !< Directory to place files in
    !> If present, output fluxes with these names
    character(len=*), optional, intent(in) :: fc_names(:)
    integer, intent(in), optional :: ixs_fc(:)   !< Oncly include these face variables

    integer                       :: lvl, bc, bn, n, n_cells, n_nodes
    integer                       :: ig, i, j, id, n_ix, c_ix, n_grids
    integer                       :: cell_ix, node_ix
    integer                       :: n_cc, n_fc
    integer, parameter            :: n_ch = a$D_num_children
    integer                       :: nodes_per_box, cells_per_box
    real(dp), allocatable         :: coords(:), cc_vars(:,:)
    integer, allocatable          :: offsets(:), connects(:)
    integer, allocatable          :: cell_types(:), icc_used(:), ifc_used(:)
    type(vtk_t)                   :: vtkf
    character(len=400)            :: fname
    character(len=100), allocatable :: var_names(:)
#if $D == 3
    integer                       :: k, bn2
#endif

    if (.not. tree%ready) stop "Tree not ready"
    if (present(ixs_cc)) then
       if (maxval(ixs_cc) > tree%n_var_cell .or. &
            minval(ixs_cc) < 1) stop "a$D_write_vtk: wrong indices given (ixs_cc)"
       if (size(ixs_cc) /= size(cc_names)) &
            stop "a$D_write_vtk: size(cc_names) /= size(ixs_cc)"
       icc_used = ixs_cc
    else
       if (size(cc_names) /= tree%n_var_cell) &
            stop "a$D_write_vtk: size(cc_names) /= n_var_cell"
       icc_used = [(i, i = 1, tree%n_var_cell)]
    end if

    if (present(fc_names)) then
       if (.not. present(ixs_fc)) then
          stop "a$D_write_vtk: ixs_fc not present (but fc_names is)"
       else
          if (size(ixs_fc) * $D /= size(fc_names)) then
             stop "a$D_write_vtk: size(fc_names) /= size(ixs_fc) * $D"
          end if
       end if
       ifc_used = ixs_fc
    else
       allocate(ifc_used(0))
    end if

    n_cc = size(icc_used)
    n_fc = size(ifc_used)

    allocate(var_names(n_cc + n_fc * $D))
    var_names(1:n_cc) = cc_names
    if (present(fc_names)) var_names(n_cc+1:) = fc_names

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
    allocate(cc_vars(n_cells, n_cc + n_fc * $D))
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
                cc_vars(c_ix, 1:n_cc)     = tree%boxes(id)%cc(i, j, icc_used)
                cc_vars(c_ix, n_cc+1::$D) = &
                     0.5_dp * sum(tree%boxes(id)%fx(i:i+1, j, ifc_used))
                cc_vars(c_ix, n_cc+2::$D) = &
                     0.5_dp * sum(tree%boxes(id)%fy(i, j:j+1, ifc_used))
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
                   cc_vars(c_ix, 1:n_cc)     = tree%boxes(id)%cc(i, j, k, icc_used)
                   cc_vars(c_ix, n_cc+1::$D) = &
                        0.5_dp * sum(tree%boxes(id)%fx(i:i+1, j, k, ifc_used))
                   cc_vars(c_ix, n_cc+2::$D) = &
                        0.5_dp * sum(tree%boxes(id)%fy(i, j:j+1, k, ifc_used))
                   cc_vars(c_ix, n_cc+3::$D) = &
                        0.5_dp * sum(tree%boxes(id)%fz(i, j, k:k+1, ifc_used))
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

    fname = trim(filename) // ".vtu"

    if (present(dir)) then
       i = len_trim(dir)
       if (i > 0) then
          if (dir(i:i) == "/") then ! Dir has trailing slash
             fname = trim(dir) // trim(fname)
          else
             fname = trim(dir) // "/" // trim(fname)
          end if
       end if
    end if

    call vtk_ini_xml(vtkf, trim(fname), 'UnstructuredGrid')
    call vtk_dat_xml(vtkf, "UnstructuredGrid", .true.)
    call vtk_geo_xml(vtkf, coords, n_nodes, n_cells, $D, n_cycle, time)
    call vtk_con_xml(vtkf, connects, offsets, cell_types, n_cells)
    call vtk_dat_xml(vtkf, "CellData", .true.)

    do n = 1, n_cc + n_fc * $D
       call vtk_var_r8_xml(vtkf, trim(var_names(n)), cc_vars(:, n), n_cells)
    end do

    call vtk_dat_xml(vtkf, "CellData", .false.)
    call vtk_geo_xml_close(vtkf)
    call vtk_dat_xml(vtkf, "UnstructuredGrid", .false.)
    call vtk_end_xml(vtkf)
    print *, "Written ", trim(fname), ", n_grids", n_grids
  end subroutine a$D_write_vtk

  !> Write the cell centered data of a tree to a Silo file. Only the
  !> leaves of the tree are used
  subroutine a$D_write_silo(tree, filename, cc_names, n_cycle, time, ixs_cc, &
       fc_names, ixs_fc, dir)
    use m_write_silo
    type(a$D_t), intent(in)       :: tree        !< Tree to write out
    character(len=*)              :: filename    !< Filename for the vtk file
    character(len=*)              :: cc_names(:) !< Names of the cell-centered variables
    integer, intent(in)           :: n_cycle     !< Cycle-number for vtk file (counter)
    real(dp), intent(in)          :: time        !< Time for output file
    integer, intent(in), optional :: ixs_cc(:)      !< Oncly include these cell variables
    character(len=*), optional, intent(in) :: dir !< Directory to place files in
    !> If present, output fluxes with these names
    character(len=*), optional, intent(in) :: fc_names(:)
    integer, intent(in), optional :: ixs_fc(:)      !< Oncly include these face variables

    character(len=*), parameter     :: grid_name = "gg"
    character(len=*), parameter     :: amr_name  = "mesh", meshdir = "data"
    character(len=100), allocatable :: grid_list(:), var_list(:, :), var_names(:)
    character(len=400)              :: fname
    integer                         :: lvl, i, id, i_grid, iv, nc, n_grids_max
    integer                         :: n_vars, i0, j0, dbix, n_cc, n_fc
    integer                         :: nx, ny, nx_prev, ny_prev, ix, iy
    integer, allocatable            :: ids(:), nb_ids(:), icc_used(:), ifc_used(:)
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

    if (present(ixs_cc)) then
       if (maxval(ixs_cc) > tree%n_var_cell .or. &
            minval(ixs_cc) < 1) stop "a$D_write_silo: wrong indices given (ixs_cc)"
       if (size(ixs_cc) /= size(cc_names)) &
            stop "a$D_write_silo: size(cc_names) /= size(ixs_cc)"
       icc_used = ixs_cc
    else
       if (size(cc_names) /= tree%n_var_cell) &
            stop "a$D_write_silo: size(cc_names) /= n_var_cell"
       icc_used = [(i, i = 1, tree%n_var_cell)]
    end if

    if (present(fc_names)) then
       if (.not. present(ixs_fc)) then
          stop "a$D_write_vtk: ixs_fc not present (but fc_names is)"
       else
          if (size(ixs_fc) * $D /= size(fc_names)) then
             stop "a$D_write_vtk: size(fc_names) /= size(ixs_fc) * $D"
          end if
       end if
       ifc_used = ixs_fc
    else
       allocate(ifc_used(0))
    end if

    n_cc = size(icc_used)
    n_fc = size(ifc_used)

    allocate(var_names(n_cc + n_fc * $D))
    var_names(1:n_cc) = cc_names
    if (present(fc_names)) var_names(n_cc+1:) = fc_names

    nc = tree%n_cell
    n_vars = n_cc + n_fc * $D
    n_grids_max = 0
    do lvl = 1, tree%max_lvl
       n_grids_max = n_grids_max + size(tree%lvls(lvl)%leaves)
    end do

    allocate(grid_list(n_grids_max))
    allocate(var_list(n_vars, n_grids_max))
    allocate(box_done(tree%max_id))
    box_done = .false.

    fname = trim(filename) // ".silo"

    if (present(dir)) then
       i = len_trim(dir)
       if (i > 0) then
          if (dir(i:i) == "/") then ! Dir has trailing slash
             fname = trim(dir) // trim(fname)
          else
             fname = trim(dir) // "/" // trim(fname)
          end if
       end if
    end if

    call SILO_create_file(trim(fname), dbix)
    call SILO_mkdir(dbix, meshdir)
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
          box_done(id) = .true.
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
                     tree%boxes(nb_ids(1))%ix(1) < tree%boxes(ids(1))%ix(1) .and. &
                     .not. any(a2_has_children(tree%boxes(nb_ids)))) then
                   nx = nx + 1
                   allocate(new_box_list(nx, ny))
                   new_box_list(1, :) = nb_ids
                   new_box_list(2:, :) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the +x direction
             ids = box_list(nx, :)
             nb_ids = tree%boxes(ids)%neighbors(a2_nb_hx)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(1) > tree%boxes(ids(1))%ix(1) .and. &
                     .not. any(a2_has_children(tree%boxes(nb_ids)))) then
                   nx = nx + 1
                   allocate(new_box_list(nx, ny))
                   new_box_list(nx, :) = nb_ids
                   new_box_list(1:nx-1, :) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the -y direction
             ids = box_list(:, 1)
             nb_ids = tree%boxes(ids)%neighbors(a2_nb_ly)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(2) < tree%boxes(ids(1))%ix(2) .and. &
                     .not. any(a2_has_children(tree%boxes(nb_ids)))) then
                   ny = ny + 1
                   allocate(new_box_list(nx, ny))
                   new_box_list(:, 1) = nb_ids
                   new_box_list(:, 2:) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the +y direction
             ids = box_list(:, ny)
             nb_ids = tree%boxes(ids)%neighbors(a2_nb_hy)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(2) > tree%boxes(ids(1))%ix(2) .and. &
                     .not. any(a2_has_children(tree%boxes(nb_ids)))) then
                   ny = ny + 1
                   allocate(new_box_list(nx, ny))
                   new_box_list(:, ny) = nb_ids
                   new_box_list(:, 1:ny-1) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             if (nx == nx_prev .and. ny == ny_prev) exit
          end do

          allocate(var_data(nx * nc, ny * nc, n_vars))
          do ix = 1, nx
             do iy = 1, ny
                id = box_list(ix, iy)
                i0 = 1 + (ix-1) * nc
                j0 = 1 + (iy-1) * nc
                var_data(i0:i0+nc-1, j0:j0+nc-1, 1:n_cc) = &
                     tree%boxes(id)%cc(1:nc, 1:nc, icc_used)
                var_data(i0:i0+nc-1, j0:j0+nc-1, n_cc+1::$D) = &
                     0.5_dp * (tree%boxes(id)%fx(1:nc, 1:nc, ifc_used) + &
                     tree%boxes(id)%fx(2:nc+1, 1:nc, ifc_used))
                var_data(i0:i0+nc-1, j0:j0+nc-1, n_cc+2::$D) = &
                     0.5_dp * (tree%boxes(id)%fy(1:nc, 1:nc, ifc_used) + &
                     tree%boxes(id)%fy(1:nc, 2:nc+1, ifc_used))
             end do
          end do

          id = box_list(1, 1)
          dr = tree%boxes(id)%dr
          r_min = tree%boxes(id)%r_min

          write(grid_list(i_grid), "(A,I0)") meshdir // '/' // grid_name, i_grid
          call SILO_add_grid(dbix, grid_list(i_grid), 2, &
               [nx*nc + 1, ny*nc + 1], r_min, dr)
          do iv = 1, n_vars
             write(var_list(iv, i_grid), "(A,I0)") meshdir // '/' // &
                  trim(var_names(iv)) // "_", i_grid
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
                     tree%boxes(nb_ids(1))%ix(1) < tree%boxes(ids(1))%ix(1) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   nx = nx + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(1, :, :) = reshape(nb_ids, [ny, nz])
                   new_box_list(2:, :, :) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the +x direction
             ids = pack(box_list(nx, :, :), .true.)
             nb_ids = tree%boxes(ids)%neighbors(a3_nb_hx)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(1) > tree%boxes(ids(1))%ix(1) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   nx = nx + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(nx, :, :) = reshape(nb_ids, [ny, nz])
                   new_box_list(1:nx-1, :, :) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the -y direction
             ids = pack(box_list(:, 1, :), .true.)
             nb_ids = tree%boxes(ids)%neighbors(a3_nb_ly)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(2) < tree%boxes(ids(1))%ix(2) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   ny = ny + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(:, 1, :) = reshape(nb_ids, [nx, nz])
                   new_box_list(:, 2:, :) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the +y direction
             ids = pack(box_list(:, ny, :), .true.)
             nb_ids = tree%boxes(ids)%neighbors(a3_nb_hy)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(2) > tree%boxes(ids(1))%ix(2) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   ny = ny + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(:, ny, :) = reshape(nb_ids, [nx, nz])
                   new_box_list(:, 1:ny-1, :) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the -z direction
             ids = pack(box_list(:, :, 1), .true.)
             nb_ids = tree%boxes(ids)%neighbors(a3_nb_lz)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(3) < tree%boxes(ids(1))%ix(3) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   nz = nz + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(:, :, 1) = reshape(nb_ids, [nx, ny])
                   new_box_list(:, :, 2:) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the +z direction
             ids = pack(box_list(:, :, nz), .true.)
             nb_ids = tree%boxes(ids)%neighbors(a3_nb_hz)
             if (all(nb_ids > a5_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(3) > tree%boxes(ids(1))%ix(3) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   nz = nz + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(:, :, nz) = reshape(nb_ids, [nx, ny])
                   new_box_list(:, :, 1:nz-1) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
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
                   i0 = 1 + (ix-1) * nc
                   j0 = 1 + (iy-1) * nc
                   k0 = 1 + (iz-1) * nc
                   var_data(i0:i0+nc-1, j0:j0+nc-1, k0:k0+nc-1, 1:n_cc) = &
                        tree%boxes(id)%cc(1:nc, 1:nc, 1:nc, icc_used)
                   var_data(i0:i0+nc-1, j0:j0+nc-1, k0:k0+nc-1, n_cc+1::$D) = &
                        0.5_dp * (tree%boxes(id)%fx(1:nc, 1:nc, 1:nc, ifc_used) + &
                        tree%boxes(id)%fx(2:nc+1, 1:nc, 1:nc, ifc_used))
                   var_data(i0:i0+nc-1, j0:j0+nc-1, k0:k0+nc-1, n_cc+2::$D) = &
                        0.5_dp * (tree%boxes(id)%fy(1:nc, 1:nc, 1:nc, ifc_used) + &
                        tree%boxes(id)%fy(1:nc, 2:nc+1, 1:nc, ifc_used))
                   var_data(i0:i0+nc-1, j0:j0+nc-1, k0:k0+nc-1, n_cc+3::$D) = &
                        0.5_dp * (tree%boxes(id)%fz(1:nc, 1:nc, 1:nc, ifc_used) + &
                        tree%boxes(id)%fz(1:nc, 1:nc, 2:nc+1, ifc_used))
                end do
             end do
          end do

          id = box_list(1, 1, 1)
          dr = tree%boxes(id)%dr
          r_min = tree%boxes(id)%r_min

          write(grid_list(i_grid), "(A,I0)") meshdir // '/' // grid_name, i_grid
          call SILO_add_grid(dbix, grid_list(i_grid), 3, &
               [nx*nc + 1, ny*nc + 1, nz*nc + 1], r_min, dr)
          do iv = 1, n_vars
             write(var_list(iv, i_grid), "(A,I0)") meshdir // '/' // &
                  trim(var_names(iv)) // "_", i_grid
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
       call SILO_set_mmesh_var(dbix, trim(var_names(iv)), amr_name, &
            var_list(iv, 1:i_grid), n_cycle, time)
    end do
    call SILO_close_file(dbix)

    print *, "Written ", trim(fname), ", n_grids", i_grid
  end subroutine a$D_write_silo

end module m_afivo_$Dd

