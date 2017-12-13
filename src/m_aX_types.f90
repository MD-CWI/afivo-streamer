!> This module contains the basic types and constants that are used in the
!> $D-dimensional version of Afivo, together with some basic routines. The
!> dimension-independent types and constant are place in m_afivo_types.
module m_a$D_types

  ! Import dimension-independent types
  use m_afivo_types

  implicit none
  public

#if $D == 2
  ! Numbering of children (same location as **corners**)
  integer, parameter :: a2_num_children = 4
  integer, parameter :: a2_child_lowx_lowy = 1
  integer, parameter :: a2_child_highx_lowy = 2
  integer, parameter :: a2_child_lowx_highy = 3
  integer, parameter :: a2_child_highx_highy = 4

  ! Neighboring indices for each child
  integer, parameter :: a2_child_neighbs(2, 4) = reshape([1,3,2,3,1,4,2,4], [2,4])
  ! Index offset for each child
  integer, parameter :: a2_child_dix(2, 4) = reshape([0,0,1,0,0,1,1,1], [2,4])
  ! Reverse child index in each direction
  integer, parameter :: a2_child_rev(4, 2) = reshape([2,1,4,3,3,4,1,2], [4,2])
  ! Children adjacent to a neighbor
  integer, parameter :: a2_child_adj_nb(2, 4) = reshape([1,3,2,4,1,2,3,4], [2,4])
  ! Neighbors adjacent to a child
  integer, parameter :: a2_nb_adj_child(2, 4) = reshape([1,3,2,3,1,4,2,4], [2,4])
  ! Which children have a low index per dimension
  logical, parameter :: a2_child_low(2, 4) = reshape([.true., .true., &
       .false., .true., .true., .false., .false., .false.], [2, 4])
  integer, parameter :: a2_child_high_01(2, 4) = &
       reshape([0, 0, 1, 0, 0, 1, 1, 1], [2, 4])

  ! Neighbor topology information
  integer, parameter :: a2_num_neighbors = 4
  integer, parameter :: a2_neighb_lowx = 1
  integer, parameter :: a2_neighb_highx = 2
  integer, parameter :: a2_neighb_lowy = 3
  integer, parameter :: a2_neighb_highy = 4

  ! Index offsets of neighbors
  integer, parameter :: a2_neighb_dix(2, 4) = reshape([-1,0,1,0,0,-1,0,1], [2,4])
  ! Which neighbors have a lower index
  logical, parameter :: a2_neighb_low(4) = [.true., .false., .true., .false.]
  ! The low neighbors
  integer, parameter :: a2_low_neighbs(2) = [1, 3]
  ! The high neighbors
  integer, parameter :: a2_high_neighbs(2) = [2, 4]
  ! Opposite of nb_low, but now as 0,1 integers
  integer, parameter :: a2_neighb_high_01(4) = [0, 1, 0, 1]
  ! Opposite of nb_low, but now as -1,1 integers
  integer, parameter :: a2_neighb_high_pm(4) = [-1, 1, -1, 1]

  ! Reverse neighbors
  integer, parameter :: a2_neighb_rev(4) = [2, 1, 4, 3]
  ! Direction (dimension) for a neighbor
  integer, parameter :: a2_neighb_dim(4) = [1, 1, 2, 2]
#elif $D == 3
  ! Numbering of children (same location as **corners**)
  integer, parameter :: a3_num_children = 8
  integer, parameter :: a3_child_lowx_lowy_lowz = 1
  integer, parameter :: a3_child_highx_lowy_lowz = 2
  integer, parameter :: a3_child_lowx_highy_lowz = 3
  integer, parameter :: a3_child_highx_highy_lowz = 4
  integer, parameter :: a3_child_lowx_lowy_highz = 5
  integer, parameter :: a3_child_highx_lowy_highz = 6
  integer, parameter :: a3_child_lowx_highy_highz = 7
  integer, parameter :: a3_child_highx_highy_highz = 8

  ! Neighboring indices for each child
  integer, parameter :: a3_child_neighbs(3, 8) = reshape( &
       [1,3,5, 2,3,5, 1,4,5, 2,4,5, &
       1,3,6, 2,3,6, 1,4,6, 2,4,6], [3,8])
  ! Index offset for each child
  integer, parameter :: a3_child_dix(3, 8) = reshape( &
       [0,0,0, 1,0,0, 0,1,0, 1,1,0, &
       0,0,1, 1,0,1, 0,1,1, 1,1,1], [3,8])
  ! Reverse child index in each direction
  integer, parameter :: a3_child_rev(8, 3) = reshape( &
       [2,1,4,3,6,5,8,7, 3,4,1,2,7,8,5,6, 5,6,7,8,1,2,3,4], [8,3])
  ! Children adjacent to a neighbor
  integer, parameter :: a3_child_adj_nb(4, 6) = reshape( &
       [1,3,5,7, 2,4,6,8, 1,2,5,6, 3,4,7,8, 1,2,3,4, 5,6,7,8], [4,6])
  ! Neighbors adjacent to a child
  integer, parameter :: a3_nb_adj_child(3, 8) = reshape( &
       [1,3,5, 2,3,5, 1,4,5, 2,4,5, 1,3,6, 2,3,6, 1,4,6, 2,4,6], [3,8])
  ! Which children have a low index per dimension
  logical, parameter :: a3_child_low(3, 8) = reshape([ &
       .true., .true., .true., .false., .true., .true., &
       .true., .false., .true., .false., .false., .true., &
       .true., .true., .false., .false., .true., .false., &
       .true., .false., .false., .false., .false., .false.], [3, 8])
  ! Which children have a high index per dimension
  integer, parameter :: a3_child_high_01(3, 8) = reshape([ &
       0, 0, 0, 1, 0, 0, &
       0, 1, 0, 1, 1, 0, &
       0, 0, 1, 1, 0, 1, &
       0, 1, 1, 1, 1, 1], [3, 8])

  ! Neighbor topology information
  integer, parameter :: a3_num_neighbors = 6
  integer, parameter :: a3_neighb_lowx = 1
  integer, parameter :: a3_neighb_highx = 2
  integer, parameter :: a3_neighb_lowy = 3
  integer, parameter :: a3_neighb_highy = 4
  integer, parameter :: a3_neighb_lowz = 5
  integer, parameter :: a3_neighb_highz = 6
  ! Index offsets of neighbors
  integer, parameter :: a3_neighb_dix(3, 6) = reshape( &
       [-1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1], [3,6])
  ! Which neighbors have a lower index
  logical, parameter :: a3_neighb_low(6) = &
       [.true., .false., .true., .false., .true., .false.]
  ! The low neighbors
  integer, parameter :: a3_low_neighbs(3) = [1, 3, 5]
  ! The high neighbors
  integer, parameter :: a3_high_neighbs(3) = [2, 4, 6]
  ! Opposite of nb_low, but now as 0,1 integers
  integer, parameter :: a3_neighb_high_01(6) = [0, 1, 0, 1, 0, 1]
  ! Opposite of nb_low, but now as -1,1 integers
  integer, parameter :: a3_neighb_high_pm(6) = [-1, 1, -1, 1, -1, 1]
  ! Reverse neighbors
  integer, parameter :: a3_neighb_rev(6) = [2, 1, 4, 3, 6, 5]
  ! Direction (dimension) for a neighbor
  integer, parameter :: a3_neighb_dim(6) = [1, 1, 2, 2, 3, 3]

  ! Number of edgse
  integer, parameter :: a3_num_edges = 12
  ! Coordinate parallel to edge
  integer, parameter :: a3_edge_dim(12) = &
       [1,1,1,1, 2,2,2,2, 3,3,3,3]
  ! Direction of edge
  integer, parameter :: a3_edge_dir(3, 12) = reshape( &
       [0,-1,-1, 0,1,-1, 0,-1,1, 0,1,1, &
       -1,0,-1, 1,0,-1, -1,0,1, 1,0,1, &
       -1,-1,0, 1,-1,0, -1,1,0, 1,1,0], [3, 12])
  ! Neighbors adjacent to edges
  integer, parameter :: a3_nb_adj_edge(2, 12) = reshape( &
       [3,5, 4,5, 3,6, 4,6, &
       1,5, 2,5, 1,6, 2,6, &
       1,3, 2,3, 1,4, 2,4], [2,12])
  ! Minimum index of edge (1 indicates n_cell + 1)
  integer, parameter :: a3_edge_min_ix(3, 12) = reshape( &
       [0,0,0, 0,1,0, 0,0,1, 0,1,1, &
       0,0,0, 1,0,0, 0,0,1, 1,0,1, &
       0,0,0, 1,0,0, 0,1,0, 1,1,0], [3,12])
#endif

  !> Collection of methods for a cell-centered variable
  type a$D_cc_methods
     !> Prolongation method
     procedure(a$D_subr_prolong), pointer, nopass :: prolong => null()
     !> Restriction method
     procedure(a$D_subr_restrict), pointer, nopass :: restrict => null()
     !> Boundary condition routine
     procedure(a$D_subr_bc), pointer, nopass :: bc => null()
     !> Refinement boundary routine
     procedure(a$D_subr_rb), pointer, nopass :: rb => null()
  end type a$D_cc_methods

  !> The basic building block of afivo: a box with cell-centered and face
  !> centered data, and information about its position, neighbors, children etc.
  type box$D_t
     integer               :: lvl    !< level of the box
     logical               :: in_use=.false.  !< is the box in use?
     integer               :: tag=af_init_tag !< for the user
     integer               :: ix($D) !< index in the domain
     integer               :: parent !< index of parent in box list
     !> index of children in box list
     integer               :: children(a$D_num_children)
     !> index of neighbors in box list
     integer               :: neighbors(a$D_num_neighbors)
     !> matrix with neighbors (including diagonal ones)
#if $D == 2
     integer               :: neighbor_mat(-1:1, -1:1)
#elif $D == 3
     integer               :: neighbor_mat(-1:1, -1:1, -1:1)
#endif
     integer               :: n_cell    !< number of cells per dimension
     real(dp)              :: dr        !< width/height of a cell
     real(dp)              :: r_min($D) !< min coords. of box
     integer               :: coord_t   !< Coordinate type (e.g. Cartesian)
#if $D == 2
     real(dp), allocatable :: cc(:, :, :) !< cell centered variables
     real(dp), allocatable :: fc(:, :, :, :) !< face centered variables
#elif $D == 3
     real(dp), allocatable :: cc(:, :, :, :) !< cell centered variables
     real(dp), allocatable :: fc(:, :, :, :, :) !< x-face centered variables
#endif
  end type box$D_t

  !> Type which stores all the boxes and levels, as well as some information
  !> about the number of boxes, variables and levels.
  type a$D_t
     logical                    :: ready = .false. !< Is tree ready for use?
     integer                    :: lvl_limit       !< maximum allowed level
     integer                    :: box_limit       !< maximum number of boxes
     integer                    :: highest_lvl     !< highest level present
     integer                    :: highest_id      !< highest box index present
     integer                    :: n_cell     !< number of cells per dimension
     integer                    :: n_var_cell !< number of cell-centered variables
     integer                    :: n_var_face !< number of face-centered variables
     integer                    :: coord_t    !< Type of coordinates
     real(dp)                   :: r_base($D) !< min. coords of box at index (1,1)
     real(dp)                   :: dr_base    !< cell spacing at lvl 1

     !> Names of cell-centered variables
     character(len=af_nlen), allocatable :: cc_names(:)
     !> Names of face-centered variables
     character(len=af_nlen), allocatable :: fc_names(:)

     !> Methods for cell-centered variables
     type(a$D_cc_methods), allocatable :: cc_methods(:)

     !> For which cell-centered variables methods have been set
     logical, allocatable :: has_cc_method(:)
     !> Indices of cell-centered variables with methods
     integer, allocatable :: cc_method_vars(:)

     type(lvl_t), allocatable   :: lvls(:)    !< list storing the tree levels
     type(box$D_t), allocatable :: boxes(:)   !< list of all boxes
  end type a$D_t

  !> Type specifying the location of a cell
  type a$D_loc_t
     integer :: id = -1         !< Id of the box that the cell is in
     integer :: ix($D) = -1     !< Index inside the box
  end type a$D_loc_t

  abstract interface

     !> Subroutine for setting refinement flags
     subroutine a$D_subr_ref(box, cell_flags)
       import
       type(box$D_t), intent(in) :: box !< Box to inspect
       !> Cell refinement flags
#if $D == 2
       integer, intent(out)      :: cell_flags(box%n_cell, box%n_cell)
#elif $D == 3
       integer, intent(out)      :: cell_flags(box%n_cell, &
            box%n_cell, box%n_cell)
#endif
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

     !> To fill ghost cells near refinement boundaries.
     subroutine a$D_subr_rb(boxes, id, nb, iv)
       import
       type(box$D_t), intent(inout) :: boxes(:) !< Array with all boxes
       integer, intent(in)         :: id       !< Id of the box that needs to have ghost cells filled
       integer, intent(in)         :: nb       !< Neighbor direction in which ghost cells need to be filled
       integer, intent(in)         :: iv       !< Variable for which ghost cells are filled
     end subroutine a$D_subr_rb

     !> To fill ghost cells near physical boundaries
     subroutine a$D_subr_bc(box, nb, iv, bc_type)
       import
       type(box$D_t), intent(inout)  :: box     !< Box that needs b.c.
       integer, intent(in)          :: nb      !< Direction
       integer, intent(in)          :: iv      !< Index of variable
       integer, intent(out)         :: bc_type !< Type of b.c.
     end subroutine a$D_subr_bc

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

     !> Subroutine for prolongation
     subroutine a$D_subr_prolong(box_p, box_c, iv, iv_to, add)
       import
       type(box$D_t), intent(in)     :: box_p !< Parent box
       type(box$D_t), intent(inout)  :: box_c !< Child box
       integer, intent(in)           :: iv    !< Variable to fill
       integer, intent(in), optional :: iv_to !< Destination variable
       logical, intent(in), optional :: add   !< Add to old values
     end subroutine a$D_subr_prolong

     !> Subroutine for restriction
     subroutine a$D_subr_restrict(box_c, box_p, iv, iv_to)
       import
       type(box$D_t), intent(in)     :: box_c !< Child box to restrict
       type(box$D_t), intent(inout)  :: box_p !< Parent box to restrict to
       integer, intent(in)           :: iv    !< Variable to restrict
       integer, intent(in), optional :: iv_to !< Destination (if /= iv)
     end subroutine a$D_subr_restrict
  end interface

contains

  !> Get tree info
  subroutine a$D_print_info(tree)
    type(a$D_t), intent(in) :: tree !< The tree

    if (.not. allocated(tree%lvls)) then
       print *, "a$D_init has not been called for this tree"
    else if (.not. tree%ready) then
       print *, "a$D_set_base has not been called for this tree"
    else
       write(*, "(A)") ""
       write(*, "(A)") " *** a$D_print_info(tree) ***"
       write(*, "(A,I10)") " Current maximum level:  ", tree%highest_lvl
       write(*, "(A,I10)") " Maximum allowed level:  ", tree%lvl_limit
       write(*, "(A,I10)") " Number of boxes used:   ", a$D_num_boxes_used(tree)
       write(*, "(A,I10)") " Number of leaves used:  ", a$D_num_leaves_used(tree)
       write(*, "(A,I10)") " Memory limit(boxes):    ", tree%box_limit
       write(*, "(A,E10.2)") " Memory limit(GByte):    ", &
            tree%box_limit * 0.5_dp**30 * &
            a$D_box_bytes(tree%n_cell, tree%n_var_cell, tree%n_var_face)
       write(*, "(A,I10)") " Highest id in box list: ", tree%highest_id
       write(*, "(A,I10)") " Size of box list:       ", size(tree%boxes)
       write(*, "(A,I10)") " Box size (cells):       ", tree%n_cell
       write(*, "(A,I10)") " Number of cc variables: ", tree%n_var_cell
       write(*, "(A,I10)") " Number of fc variables: ", tree%n_var_face
       write(*, "(A,A15)") " Type of coordinates:    ", &
            af_coord_names(tree%coord_t)
       write(*, "(A,$DE12.4)") " min. coords:        ", tree%r_base
       write(*, "(A,2E12.4)")  " dx at min/max level ", tree%dr_base, a$D_min_dr(tree)
       write(*, "(A)") ""
    end if
  end subroutine a$D_print_info

  pure function a$D_box_bytes(n_cell, n_var_cell, n_var_face) result(box_bytes)
    integer, intent(in) :: n_cell     !< number of cells per dimension
    integer, intent(in) :: n_var_cell !< number of cell-centered variables
    integer, intent(in) :: n_var_face !< number of face-centered variables
    integer             :: box_bytes
    type(box$D_t)       :: dummy_box

    box_bytes = 8 * n_var_cell * (n_cell + 2)**$D + &
         8 * n_var_face * (n_cell + 1) * n_cell**($D-1) + &
         int(storage_size(dummy_box) / 8)
  end function a$D_box_bytes

  pure function a$D_num_boxes_used(tree) result(n_boxes)
    type(a$D_t), intent(in) :: tree
    integer                 :: n_boxes, lvl

    n_boxes = 0
    do lvl = lbound(tree%lvls, 1), tree%highest_lvl
       n_boxes = n_boxes + size(tree%lvls(lvl)%ids)
    end do
  end function a$D_num_boxes_used

  pure function a$D_num_leaves_used(tree) result(n_boxes)
    type(a$D_t), intent(in) :: tree
    integer                 :: n_boxes, lvl

    n_boxes = 0
    do lvl = lbound(tree%lvls, 1), tree%highest_lvl
       n_boxes = n_boxes + size(tree%lvls(lvl)%leaves)
    end do
  end function a$D_num_leaves_used

  pure function a$D_num_cells_used(tree) result(n_cells)
    type(a$D_t), intent(in) :: tree
    integer                 :: n_cells

    n_cells = a$D_num_leaves_used(tree) * tree%n_cell**$D
  end function a$D_num_cells_used

  !> Return .true. if a box has children
  elemental logical function a$D_has_children(box)
    type(box$D_t), intent(in) :: box

    ! Boxes are either fully refined or not, so we only need to check one of the
    ! children
    a$D_has_children = (box%children(1) /= af_no_box)
  end function a$D_has_children

  !> Return .true. where there is a physical/periodic boundary. Detecting
  !> physical boundaries is straightforward (simply test whether neighbors(i) <
  !> af_no_box), but periodic boundaries require a comparison of their spatial
  !> index.
  pure function a$D_phys_boundary(boxes, id, nbs) result(boundary)
    type(box$D_t), intent(in) :: boxes(:) !< List of boxes
    integer, intent(in)       :: id       !< Box to inspect
    integer, intent(in)       :: nbs(:)   !< Neighbor directions
    logical                   :: boundary(size(nbs))
    integer                   :: n, nb, p_id, nb_id, dim, dix

    do n = 1, size(nbs)
       nb = nbs(n)
       nb_id = boxes(id)%neighbors(nb)

       if (nb_id < af_no_box) then
          boundary(n) = .true.  ! Physical boundary
       else                     ! There could be a periodic boundary
          dim = a$D_neighb_dim(nb)

          if (nb_id == af_no_box) then
             ! Refinement boundary, compute index difference on parent (which is
             ! guaranteed to be there)
             p_id = boxes(id)%parent
             nb_id = boxes(p_id)%neighbors(nb)
             dix = boxes(nb_id)%ix(dim) - boxes(p_id)%ix(dim)
          else
             ! Normal neighbor, compute index difference
             dix = boxes(nb_id)%ix(dim) - boxes(id)%ix(dim)
          end if

          if (dix /= a$D_neighb_high_pm(nb)) then
             ! The difference in index is not the expected +-1, so a periodic
             ! boundary
             boundary(n) = .true.
          else
             boundary(n) = .false.
          end if
       end if
    end do

  end function a$D_phys_boundary

  !> Get the offset of a box with respect to its parent (e.g. in 2d, there can
  !> be a child at offset 0,0, one at n_cell/2,0, one at 0,n_cell/2 and one at
  !> n_cell/2, n_cell/2)
  pure function a$D_get_child_offset(box, nb) result(ix_offset)
    type(box$D_t), intent(in)           :: box   !< A child box
    integer, intent(in), optional      :: nb     !< Optional: get index on parent neighbor
    integer                            :: ix_offset($D)
    if (box%lvl > 1) then
       ix_offset = iand(box%ix-1, 1) * ishft(box%n_cell, -1) ! * n_cell / 2
       if (present(nb)) ix_offset = ix_offset - a$D_neighb_dix(:, nb) * box%n_cell
    else                        ! In the subtree, parents are half the size
       ix_offset = 0
       if (present(nb)) ix_offset = ix_offset - &
            a$D_neighb_dix(:, nb) * ishft(box%n_cell, -1) ! n_cell / 2
    endif
  end function a$D_get_child_offset

  !> Given a cell index on box, get index of the closest cell at its parent
  pure function a$D_get_ix_on_parent(box, ix) result(p_ix)
    type(box$D_t), intent(in) :: box    !< A child box
    integer, intent(in)       :: ix($D) !< Index on child box
    integer                   :: p_ix($D)
    p_ix = a$D_get_child_offset(box) + ishft(ix+1, -1)
  end function a$D_get_ix_on_parent

  !> Given a cell index on box, get index on a neighbor
  pure function a$D_get_ix_on_neighb(box, ix, nb) result(nb_ix)
    type(box$D_t), intent(in) :: box    !< A box
    integer, intent(in)       :: ix($D) !< Index on box
    integer, intent(in)       :: nb     !< Neighbor identifier
    integer                   :: nb_ix($D), nb_dim

    nb_dim        = a$D_neighb_dim(nb)
    nb_ix         = ix
    nb_ix(nb_dim) = nb_ix(nb_dim) - a$D_neighb_high_pm(nb) * box%n_cell
  end function a$D_get_ix_on_neighb

  !> Get diagonal neighbors. Returns the index of the neighbor if found,
  !> otherwise the result nb_id <= af_no_box.
  pure function a$D_diag_neighb_id(boxes, id, nbs) result(nb_id)
    type(box$D_t), intent(in) :: boxes(:) !< List of all the boxes
    integer, intent(in)       :: id       !< Start index
    integer, intent(in)       :: nbs(:)   ! List of neighbor directions
    integer                   :: i, j, k, nb, nb_id
    integer                   :: nbs_perm(size(nbs))

    if (size(nbs) == 0) then
       nb_id = id
    else
       do i = 1, size(nbs)
          nb_id = id

          ! Check if path exists starting from nbs(i)
          do j = 1, size(nbs)
             ! k starts at i and runs over the neighbors
             k = 1 + mod(i + j - 2, size(nbs))
             nb = nbs(k)

             nb_id = boxes(nb_id)%neighbors(nb)
             if (nb_id <= af_no_box) exit
          end do

          if (nb_id > af_no_box) exit ! Found it
       end do
    end if

    ! For a corner neighbor in 3D, try again using the permuted neighbor list to
    ! covers all paths
    if (size(nbs) == 3 .and. nb_id <= af_no_box) then
       nbs_perm = nbs([2,1,3])

       do i = 1, size(nbs)
          nb_id = id

          do j = 1, size(nbs)
             k = 1 + mod(i + j - 2, size(nbs))
             nb = nbs(k)

             nb_id = boxes(nb_id)%neighbors(nb)
             if (nb_id <= af_no_box) exit
          end do

          if (nb_id > af_no_box) exit ! Found it
       end do
    end if

  end function a$D_diag_neighb_id

  !> Given a list of neighbor directions, compute the index offset
  pure function a$D_neighb_offset(nbs) result(dix)
    integer, intent(in) :: nbs(:) !< List of neighbor directions
    integer             :: n, dim, dix($D)

    dix = 0
    do n = 1, size(nbs)
       dim = a$D_neighb_dim(nbs(n))
       dix(dim) = dix(dim) + a$D_neighb_high_pm(nbs(n))
    end do
  end function a$D_neighb_offset

  !> Compute the 'child index' for a box with spatial index ix. With 'child
  !> index' we mean the index in the children(:) array of its parent.
  pure integer function a$D_ix_to_ichild(ix)
    integer, intent(in) :: ix($D) !< Spatial index of the box
    ! The index can range from 1 (all ix odd) and 2**$D (all ix even)
#if $D == 2
    a$D_ix_to_ichild = 4 - 2 * iand(ix(2), 1) - iand(ix(1), 1)
#elif $D == 3
    a$D_ix_to_ichild = &
         8 - 4 * iand(ix(3), 1) - 2 * iand(ix(2), 1) - iand(ix(1), 1)
#endif
  end function a$D_ix_to_ichild

  !> Get the cell index in which r lies. This routine does not check whether r
  !> is actually located inside the box.
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

  !> Get the location of the face parallel to dim with index fc_ix
  pure function a$D_r_fc(box, dim, fc_ix) result(r)
    type(box$D_t), intent(in) :: box
    integer, intent(in)       :: dim
    integer, intent(in)       :: fc_ix($D)
    real(dp)                  :: r($D)
    r      = box%r_min + (fc_ix-0.5_dp) * box%dr
    r(dim) = r(dim) - 0.5_dp * box%dr
  end function a$D_r_fc

  !> Get a general location with index cc_ix (like a$D_r_cc but using reals)
  pure function a$D_rr_cc(box, cc_ix) result(r)
    type(box$D_t), intent(in) :: box
    real(dp), intent(in)     :: cc_ix($D)
    real(dp)                 :: r($D)
    r = box%r_min + (cc_ix-0.5_dp) * box%dr
  end function a$D_rr_cc

  !> Return the coordinate of the center of a box
  pure function a$D_r_center(box) result(r_center)
    type(box$D_t), intent(in) :: box
    real(dp)                 :: r_center($D)
    r_center = box%r_min + 0.5_dp * box%n_cell * box%dr
  end function a$D_r_center

  !> Return finest dr that is used in the tree
  pure function a$D_min_dr(tree) result(dr)
    type(a$D_t), intent(in) :: tree
    real(dp)               :: dr !< Output: dr at the finest lvl of the tree
    dr = a$D_lvl_dr(tree, tree%highest_lvl)
  end function a$D_min_dr

  !> Return dr at lvl
  pure function a$D_lvl_dr(tree, lvl) result(dr)
    type(a$D_t), intent(in) :: tree
    integer, intent(in)    :: lvl
    real(dp)               :: dr !< Output: dr at the finest lvl of the tree
    dr = tree%dr_base * 0.5_dp**(lvl-1)
  end function a$D_lvl_dr

  subroutine a$D_set_box_gc(box, nb, iv, gc_scalar, gc_array)
    type(box$D_t), intent(inout) :: box      !< Box to operate on
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    real(dp), optional, intent(in) :: gc_scalar !< Scalar value for ghost cells
#if $D == 2
    real(dp), optional, intent(in) :: gc_array(box%n_cell) !< Array for ghost cells
#elif $D == 3
    !> Array for ghost cells
    real(dp), optional, intent(in) :: gc_array(box%n_cell, box%n_cell)
#endif
    integer                     :: nc

    nc = box%n_cell

    if (present(gc_array)) then
       select case (nb)
#if $D == 2
       case (a2_neighb_lowx)
          box%cc(0, 1:nc, iv)    = gc_array
       case (a2_neighb_highx)
          box%cc(nc+1, 1:nc, iv) = gc_array
       case (a2_neighb_lowy)
          box%cc(1:nc, 0, iv)    = gc_array
       case (a2_neighb_highy)
          box%cc(1:nc, nc+1, iv) = gc_array
#elif $D == 3
       case (a3_neighb_lowx)
          box%cc(0, 1:nc, 1:nc, iv)    = gc_array
       case (a3_neighb_highx)
          box%cc(nc+1, 1:nc, 1:nc, iv) = gc_array
       case (a3_neighb_lowy)
          box%cc(1:nc, 0, 1:nc, iv)    = gc_array
       case (a3_neighb_highy)
          box%cc(1:nc, nc+1, 1:nc, iv) = gc_array
       case (a3_neighb_lowz)
          box%cc(1:nc, 1:nc, 0, iv)    = gc_array
       case (a3_neighb_highz)
          box%cc(1:nc, 1:nc, nc+1, iv) = gc_array
#endif
       end select
    else if (present(gc_scalar)) then
       select case (nb)
#if $D == 2
       case (a2_neighb_lowx)
          box%cc(0, 1:nc, iv)    = gc_scalar
       case (a2_neighb_highx)
          box%cc(nc+1, 1:nc, iv) = gc_scalar
       case (a2_neighb_lowy)
          box%cc(1:nc, 0, iv)    = gc_scalar
       case (a2_neighb_highy)
          box%cc(1:nc, nc+1, iv) = gc_scalar
#elif $D == 3
       case (a3_neighb_lowx)
          box%cc(0, 1:nc, 1:nc, iv)    = gc_scalar
       case (a3_neighb_highx)
          box%cc(nc+1, 1:nc, 1:nc, iv) = gc_scalar
       case (a3_neighb_lowy)
          box%cc(1:nc, 0, 1:nc, iv)    = gc_scalar
       case (a3_neighb_highy)
          box%cc(1:nc, nc+1, 1:nc, iv) = gc_scalar
       case (a3_neighb_lowz)
          box%cc(1:nc, 1:nc, 0, iv)    = gc_scalar
       case (a3_neighb_highz)
          box%cc(1:nc, 1:nc, nc+1, iv) = gc_scalar
#endif
       end select
    else
       stop "a$D_set_box_gc: requires gc_scalar or gc_array"
    end if
  end subroutine a$D_set_box_gc

#if $D == 2
  !> Get the radius of the cell center with first index i
  pure function a$D_cyl_radius_cc(box, i) result(r)
    type(box$D_t), intent(in) :: box
    integer, intent(in)      :: i
    real(dp)                 :: r
    r = box%r_min(1) + (i-0.5_dp) * box%dr
  end function a$D_cyl_radius_cc

  !> Get the normalized weights of the 'inner' and 'outer' children of a cell
  !> with index ix. Note that the cell centers of the children are located at
  !> -/+ 0.25 dr compared to the parent.
  subroutine a2_cyl_child_weights(box, i, inner, outer)
    type(box2_t), intent(in) :: box
    integer, intent(in)      :: i
    real(dp), intent(out)    :: inner, outer
    real(dp)                 :: tmp

    tmp = 0.25_dp * box%dr / a2_cyl_radius_cc(box, i)
    inner = 1 - tmp
    outer = 1 + tmp
  end subroutine a2_cyl_child_weights
#endif

end module m_a$D_types
