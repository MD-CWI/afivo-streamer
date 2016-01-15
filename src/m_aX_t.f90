! This module contains the basic types and constants that are used in Afivo.
!
! Author: Jannis Teunissen
! License: GPLv3

module m_a$D_t

  implicit none
  public

  integer, parameter :: dp = kind(0.0d0)

  ! The prefix a5_ means that this constant is the same in the 2D and 3D
  ! version of Afivo, and is inspired by "A-five-o"

  !> Value indicating you want to derefine a box
  integer, parameter :: a5_rm_ref = -1

  !> Value indicating you want to keep a box's refinement
  integer, parameter :: a5_kp_ref = 0

  !> Value indicating you want to refine a box
  integer, parameter :: a5_do_ref = 1

  !> The children of a box are removed (for internal use)
  integer, parameter :: a5_derefine = -2

  !> A box will be refined (for internal use)
  integer, parameter :: a5_refine = 2

  !> Special value indicating there is no box
  integer, parameter :: a5_no_box = 0

  !> Each box contains a tag, for which bits can be set. This is the initial
  !> value, which should not be used by the user
  integer, parameter :: a5_init_tag = -huge(1)

  !> Default coordinate system
  integer, parameter :: a5_xyz = 0

  !> Cylindrical coordinate system
  integer, parameter :: a5_cyl = 1

  !> Value to indicate a Dirichlet boundary condition
  integer, parameter :: a5_bc_dirichlet = 0

  !> Value to indicate a Neumann boundary condition
  integer, parameter :: a5_bc_neumann = 1

#if $D == 2
  ! Numbering of children (same location as **corners**)
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
  ! Numbering of children (same location as **corners**)
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

     !> To fill ghost cells near refinement boundaries.
     subroutine a$D_subr_gc(boxes, id, nb, iv)
       import
       type(box$D_t), intent(inout) :: boxes(:) !< Array with all boxes
       integer, intent(in)         :: id       !< Id of the box that needs to have ghost cells filled
       integer, intent(in)         :: nb       !< Neighbor direction in which ghost cells need to be filled
       integer, intent(in)         :: iv       !< Variable for which ghost cells are filled
     end subroutine a$D_subr_gc

     !> To fill ghost cells near physical boundaries
     subroutine a$D_subr_bc(box, nb, iv, bc_values, bc_type)
       import
       type(box$D_t), intent(in) :: box !< Box that needs b.c.
       integer, intent(in)      :: nb  !< Neighbor direction in which ghost cells need to be filled
       integer, intent(in)      :: iv  !< Variable for which ghost cells are filled
       real(dp), intent(out)    :: bc_values(box%n_cell)
       real(dp), intent(out)    :: bc_type
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
  end interface

end module m_a$D_t
