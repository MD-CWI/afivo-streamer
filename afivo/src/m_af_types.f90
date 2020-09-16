#include "cpp_macros.h"
!> This module contains the basic types and constants that are used in the
!> NDIM-dimensional version of Afivo, together with some basic routines. The
!> dimension-independent types and constant are place in m_afivo_types.
module m_af_types

  implicit none
  public

  ! dp stands for double precision (8 byte reals)
  integer, parameter :: dp = kind(0.0d0)

  !> Highest allowed refinement level
  integer, parameter :: af_max_lvl = 30

  !> Lowest allowed refinement level
  integer, parameter :: af_min_lvl = 1

  !> Maximum number of variables
  integer, parameter :: af_max_num_vars = 100

  !> Value indicating you want to derefine a box
  integer, parameter :: af_rm_ref = -1

  !> Value indicating you want to keep a box's refinement
  integer, parameter :: af_keep_ref = 0

  !> Value indicating you want to refine a box
  integer, parameter :: af_do_ref = 1

  !> The children of a box are removed (for internal use)
  integer, parameter :: af_derefine = -2

  !> A box will be refined (for internal use)
  integer, parameter :: af_refine = 2

  !> Special value indicating there is no box
  integer, parameter :: af_no_box = 0

  !> Special value indicating a physical (non-periodic) boundary
  integer, parameter :: af_phys_boundary = -1

  !> Each box contains a tag, for which bits can be set. This is the initial
  !> value, which should not be used by the user
  integer, parameter :: af_init_tag = -huge(1)

  !> Default coordinate system
  integer, parameter :: af_xyz = 1

  !> Cylindrical coordinate system
  integer, parameter :: af_cyl = 2

  !> Names of coordinate systems
  character(len=*), parameter :: af_coord_names(2) = &
       ["Cartesian  ", "Cylindrical"]

  !> Value to indicate a Dirichlet boundary condition
  integer, parameter :: af_bc_dirichlet = -10

  !> Value to indicate a Neumann boundary condition
  integer, parameter :: af_bc_neumann = -11

  !> Value to indicate a continuous boundary condition
  integer, parameter :: af_bc_continuous = -12

  !> Value to indicate a Dirichlet boundary condition in which a value is copied
  !> to the ghost cells, without any type of extrapolation. This can be useful
  !> for hyperbolic PDEs
  integer, parameter :: af_bc_dirichlet_copy = -13

  !> Maximum length of the names of variables
  integer, parameter :: af_nlen = 20

  !> Type which contains the indices of all boxes at a refinement level, as well
  !> as a list with all the "leaf" boxes and non-leaf (parent) boxes
  type lvl_t
     integer, allocatable :: ids(:)     !< indices of boxes of level
     integer, allocatable :: leaves(:)  !< all ids(:) that are leaves
     integer, allocatable :: parents(:) !< all ids(:) that have children
  end type lvl_t

  !> Type that contains the refinement changes in a level
  type ref_lvl_t
     integer, allocatable :: add(:) !< Id's of newly added boxes
  end type ref_lvl_t

  !> Type that contains the refinement changes in a tree
  type ref_info_t
     integer                      :: n_add = 0 !< Total number of added boxes
     integer                      :: n_rm  = 0 !< Total number removed boxes
     integer, allocatable         :: rm(:)     !< Ids of removed boxes
     type(ref_lvl_t), allocatable :: lvls(:)   !< Ids of added boxes per level
  end type ref_info_t

#if NDIM == 1
    !> Number of children
  integer, parameter :: af_num_children = 2

  !> Index offset for each child
  integer, parameter :: af_child_dix(1, 2) = reshape([0,1], [1,2])
  !> Children adjacent to a neighbor
  integer, parameter :: af_child_adj_nb(1, 2) = reshape([1,2], [1,2])
  !> Which children have a low index per dimension
  logical, parameter :: af_child_low(1, 2) = reshape([.true., .false.], [1, 2])
  !> Whether a child located in the 'upper' direction (1 or 0)

  !> Number of neighbors
  integer, parameter :: af_num_neighbors = 2
  !> Lower-x neighbor
  integer, parameter :: af_neighb_lowx = 1
  !> Upper-x neighbor
  integer, parameter :: af_neighb_highx = 2

  !> Index offsets of neighbors
  integer, parameter :: af_neighb_dix(1, 2) = reshape([-1,1], [1,2])
  !> Which neighbors have a lower index
  logical, parameter :: af_neighb_low(2) = [.true., .false.]
  !> The low neighbors
  integer, parameter :: af_low_neighbs(1) = [1]
  !> The high neighbors
  integer, parameter :: af_high_neighbs(1) = [2]
  !> Opposite of nb_low, but now as -1,1 integers
  integer, parameter :: af_neighb_high_pm(2) = [-1, 1]

  !> Reverse neighbors
  integer, parameter :: af_neighb_rev(2) = [2, 1]
  !> Direction (dimension) for a neighbor
  integer, parameter :: af_neighb_dim(2) = [1, 1]
#elif NDIM == 2
  !> Number of children
  integer, parameter :: af_num_children = 4

  !> Index offset for each child
  integer, parameter :: af_child_dix(2, 4) = reshape([0,0,1,0,0,1,1,1], [2,4])
  !> Children adjacent to a neighbor
  integer, parameter :: af_child_adj_nb(2, 4) = reshape([1,3,2,4,1,2,3,4], [2,4])
  !> Which children have a low index per dimension
  logical, parameter :: af_child_low(2, 4) = reshape([.true., .true., &
       .false., .true., .true., .false., .false., .false.], [2, 4])

  !> Number of neighbors
  integer, parameter :: af_num_neighbors = 4
  !> Lower-x neighbor
  integer, parameter :: af_neighb_lowx = 1
  !> Upper-x neighbor
  integer, parameter :: af_neighb_highx = 2
  !> Lower-y neighbor
  integer, parameter :: af_neighb_lowy = 3
  !> Upper-y neighbor
  integer, parameter :: af_neighb_highy = 4

  !> Index offsets of neighbors
  integer, parameter :: af_neighb_dix(2, 4) = reshape([-1,0,1,0,0,-1,0,1], [2,4])
  !> Which neighbors have a lower index
  logical, parameter :: af_neighb_low(4) = [.true., .false., .true., .false.]
  !> The low neighbors
  integer, parameter :: af_low_neighbs(2) = [1, 3]
  !> The high neighbors
  integer, parameter :: af_high_neighbs(2) = [2, 4]
  !> Opposite of nb_low, but now as -1,1 integers
  integer, parameter :: af_neighb_high_pm(4) = [-1, 1, -1, 1]

  !> Reverse neighbors
  integer, parameter :: af_neighb_rev(4) = [2, 1, 4, 3]
  !> Direction (dimension) for a neighbor
  integer, parameter :: af_neighb_dim(4) = [1, 1, 2, 2]
#elif NDIM == 3
  !> Number of children
  integer, parameter :: af_num_children = 8

  !> Index offset for each child
  integer, parameter :: af_child_dix(3, 8) = reshape( &
       [0,0,0, 1,0,0, 0,1,0, 1,1,0, &
       0,0,1, 1,0,1, 0,1,1, 1,1,1], [3,8])
  !> Children adjacent to a neighbor
  integer, parameter :: af_child_adj_nb(4, 6) = reshape( &
       [1,3,5,7, 2,4,6,8, 1,2,5,6, 3,4,7,8, 1,2,3,4, 5,6,7,8], [4,6])
  !> Which children have a low index per dimension
  logical, parameter :: af_child_low(3, 8) = reshape([ &
       .true., .true., .true., .false., .true., .true., &
       .true., .false., .true., .false., .false., .true., &
       .true., .true., .false., .false., .true., .false., &
       .true., .false., .false., .false., .false., .false.], [3, 8])

  !> Number of neighbors
  integer, parameter :: af_num_neighbors = 6
  !> Lower-x neighbor
  integer, parameter :: af_neighb_lowx = 1
  !> Upper-x neighbor
  integer, parameter :: af_neighb_highx = 2
  !> Lower-y neighbor
  integer, parameter :: af_neighb_lowy = 3
  !> Upper-y neighbor
  integer, parameter :: af_neighb_highy = 4
  !> Lower-z neighbor
  integer, parameter :: af_neighb_lowz = 5
  !> Upper-z neighbor
  integer, parameter :: af_neighb_highz = 6
  !> Index offsets of neighbors
  integer, parameter :: af_neighb_dix(3, 6) = reshape( &
       [-1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1], [3,6])
  !> Which neighbors have a lower index
  logical, parameter :: af_neighb_low(6) = &
       [.true., .false., .true., .false., .true., .false.]
  !> The low neighbors
  integer, parameter :: af_low_neighbs(3) = [1, 3, 5]
  !> The high neighbors
  integer, parameter :: af_high_neighbs(3) = [2, 4, 6]
  !> Opposite of nb_low, but now as -1,1 integers
  integer, parameter :: af_neighb_high_pm(6) = [-1, 1, -1, 1, -1, 1]
  !> Reverse neighbors
  integer, parameter :: af_neighb_rev(6) = [2, 1, 4, 3, 6, 5]
  !> Direction (dimension) for a neighbor
  integer, parameter :: af_neighb_dim(6) = [1, 1, 2, 2, 3, 3]

  !> Number of edgse
  integer, parameter :: af_num_edges = 12
  !> Coordinate parallel to edge
  integer, parameter :: af_edge_dim(12) = &
       [1,1,1,1, 2,2,2,2, 3,3,3,3]
  !> Direction of edge
  integer, parameter :: af_edge_dir(3, 12) = reshape( &
       [0,-1,-1, 0,1,-1, 0,-1,1, 0,1,1, &
       -1,0,-1, 1,0,-1, -1,0,1, 1,0,1, &
       -1,-1,0, 1,-1,0, -1,1,0, 1,1,0], [3, 12])
  !> Neighbors adjacent to edges
  integer, parameter :: af_nb_adj_edge(2, 12) = reshape( &
       [3,5, 4,5, 3,6, 4,6, &
       1,5, 2,5, 1,6, 2,6, &
       1,3, 2,3, 1,4, 2,4], [2,12])
  !> Minimum index of edge (1 indicates n_cell + 1)
  integer, parameter :: af_edge_min_ix(3, 12) = reshape( &
       [0,0,0, 0,1,0, 0,0,1, 0,1,1, &
       0,0,0, 1,0,0, 0,0,1, 1,0,1, &
       0,0,0, 1,0,0, 0,1,0, 1,1,0], [3,12])
#endif

  !> Collection of methods for a cell-centered variable
  type af_cc_methods
     !> Prolongation method
     procedure(af_subr_prolong), pointer, nopass :: prolong => null()
     !> Restriction method
     procedure(af_subr_restrict), pointer, nopass :: restrict => null()
     !> Boundary condition routine
     procedure(af_subr_bc), pointer, nopass :: bc => null()
     !> Refinement boundary routine
     procedure(af_subr_rb), pointer, nopass :: rb => null()
  end type af_cc_methods

  !> The basic building block of afivo: a box with cell-centered and face
  !> centered data, and information about its position, neighbors, children etc.
  type box_t
     integer               :: lvl    !< level of the box
     logical               :: in_use=.false.  !< is the box in use?
     integer               :: tag=af_init_tag !< for the user
     integer               :: ix(NDIM) !< index in the domain
     integer               :: parent !< index of parent in box list
     !> index of children in box list
     integer               :: children(af_num_children)
     !> index of neighbors in box list
     integer               :: neighbors(af_num_neighbors)
     !> matrix with neighbors (including diagonal ones)
     integer               :: neighbor_mat(DTIMES(-1:1))
     integer               :: n_cell    !< number of cells per dimension
     real(dp)              :: dr(NDIM)  !< width/height of a cell
     real(dp)              :: r_min(NDIM) !< min coords. of box
     integer               :: coord_t   !< Coordinate type (e.g. Cartesian)
     real(dp), allocatable :: cc(DTIMES(:), :) !< cell centered variables
     real(dp), allocatable :: fc(DTIMES(:), :, :) !< face centered variables
  end type box_t

  !> Type which stores all the boxes and levels, as well as some information
  !> about the number of boxes, variables and levels.
  type af_t
     logical  :: ready       = .false. !< Is tree ready for use?
     integer  :: box_limit             !< maximum number of boxes
     integer  :: highest_lvl = 0       !< highest level present
     integer  :: highest_id  = 0       !< highest box index present
     integer  :: n_cell                !< number of cells per dimension
     integer  :: n_var_cell  = 0       !< number of cell-centered variables
     integer  :: n_var_face  = 0       !< number of face-centered variables
     integer  :: coord_t               !< Type of coordinates
     integer  :: coarse_grid_size(NDIM) = -1 !< Size of the coarse grid (if rectangular)
     logical  :: periodic(NDIM) = .false.  !< Which dimensions are periodic
     real(dp) :: r_base(NDIM)          !< min. coords of box at index (1,1)
     real(dp) :: dr_base(NDIM)         !< cell spacing at lvl 1

     !> Names of cell-centered variables
     character(len=af_nlen) :: cc_names(af_max_num_vars)

     !> Names of face-centered variables
     character(len=af_nlen) :: fc_names(af_max_num_vars)

     !> Number of copies of the variable (for e.g. time-stepping)
     integer :: cc_num_copies(af_max_num_vars) = 1

     !> Whether to include the variable in the output
     logical :: cc_write_output(af_max_num_vars) = .true.

     !> Whether to include a cell-centered variable in binary output
     logical :: cc_write_binary(af_max_num_vars) = .true.

     !> Whether to include a face-centered variable in binary output
     logical :: fc_write_binary(af_max_num_vars) = .true.

     !> Methods for cell-centered variables
     type(af_cc_methods) :: cc_methods(af_max_num_vars)

     !> For which cell-centered variables methods have been set
     logical :: has_cc_method(af_max_num_vars) = .false.

     !> Indices of cell-centered variables with methods
     integer, allocatable :: cc_auto_vars(:)

     !> List storing the tree levels
     type(lvl_t) :: lvls(af_min_lvl:af_max_lvl)

     !> List of all boxes
     type(box_t), allocatable :: boxes(:)

     !> List of all removed boxes (that can be reused)
     integer, allocatable :: removed_ids(:)

     !> Number of removed boxes
     integer :: n_removed_ids = 0
  end type af_t

  !> Type specifying the location of a cell
  type af_loc_t
     integer :: id = -1         !< Id of the box that the cell is in
     integer :: ix(NDIM) = -1     !< Index inside the box
  end type af_loc_t

  abstract interface

     !> Subroutine for setting refinement flags
     subroutine af_subr_ref(box, cell_flags)
       import
       type(box_t), intent(in) :: box !< Box to inspect
       !> Cell refinement flags
       integer, intent(out)      :: cell_flags(DTIMES(box%n_cell))
     end subroutine af_subr_ref

     !> Subroutine that gets a box
     subroutine af_subr(box)
       import
       type(box_t), intent(inout) :: box
     end subroutine af_subr

     !> Subroutine that gets a box and an array of reals
     subroutine af_subr_arg(box, rarg)
       import
       type(box_t), intent(inout) :: box
       real(dp), intent(in)        :: rarg(:)
     end subroutine af_subr_arg

     !> Subroutine that gets a list of boxes and a box id
     subroutine af_subr_boxes(boxes, id)
       import
       type(box_t), intent(inout) :: boxes(:)
       integer, intent(in)         :: id
     end subroutine af_subr_boxes

     !> Subroutine that gets a list of boxes, an id and an array of reals
     subroutine af_subr_boxes_arg(boxes, id, rarg)
       import
       type(box_t), intent(inout) :: boxes(:)
       integer, intent(in)         :: id
       real(dp), intent(in)        :: rarg(:)
     end subroutine af_subr_boxes_arg

     !> Subroutine that gets a tree and a box id
     subroutine af_subr_tree(tree, id)
       import
       type(af_t), intent(inout) :: tree
       integer, intent(in)       :: id
     end subroutine af_subr_tree

     !> Subroutine that gets a tree, a box id and an array of reals
     subroutine af_subr_tree_arg(tree, id, rarg)
       import
       type(af_t), intent(inout) :: tree
       integer, intent(in)       :: id
       real(dp), intent(in)      :: rarg(:)
     end subroutine af_subr_tree_arg

     !> To fill ghost cells near refinement boundaries.
     subroutine af_subr_rb(boxes, id, nb, iv)
       import
       type(box_t), intent(inout) :: boxes(:) !< Array with all boxes
       integer, intent(in)         :: id       !< Id of the box that needs to have ghost cells filled
       integer, intent(in)         :: nb       !< Neighbor direction in which ghost cells need to be filled
       integer, intent(in)         :: iv       !< Variable for which ghost cells are filled
     end subroutine af_subr_rb

     !> To fill ghost cells near physical boundaries
     subroutine af_subr_bc(box, nb, iv, coords, bc_val, bc_type)
       import
       type(box_t), intent(in) :: box     !< Box that needs b.c.
       integer, intent(in)     :: nb      !< Direction
       integer, intent(in)     :: iv      !< Index of variable
       real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1)) !< Coordinates of boundary
       real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1)) !< Boundary values
       integer, intent(out)    :: bc_type !< Type of b.c.
     end subroutine af_subr_bc

     !> Subroutine for prolongation
     subroutine af_subr_prolong(box_p, box_c, iv, iv_to, add)
       import
       type(box_t), intent(in)     :: box_p !< Parent box
       type(box_t), intent(inout)  :: box_c !< Child box
       integer, intent(in)           :: iv    !< Variable to fill
       integer, intent(in), optional :: iv_to !< Destination variable
       logical, intent(in), optional :: add   !< Add to old values
     end subroutine af_subr_prolong

     !> Subroutine for restriction
     subroutine af_subr_restrict(box_c, box_p, ivs, use_geometry)
       import
       type(box_t), intent(in)       :: box_c        !< Child box to restrict
       type(box_t), intent(inout)    :: box_p        !< Parent box to restrict to
       integer, intent(in)           :: ivs(:)       !< Variables to restrict
       logical, intent(in), optional :: use_geometry !< If set to false, don't use geometry
     end subroutine af_subr_restrict
  end interface

contains

  !> Get number of threads
  integer function af_get_max_threads()
    use omp_lib, only: omp_get_max_threads
    af_get_max_threads = omp_get_max_threads()
  end function af_get_max_threads

  !> Get tree info
  subroutine af_print_info(tree)
    type(af_t), intent(in) :: tree !< The tree
    character(len=15) :: fmt_string

    if (.not. allocated(tree%boxes)) then
       print *, "af_init has not been called for this tree"
    else if (.not. tree%ready) then
       print *, "af_set_base has not been called for this tree"
    else
       write(*, "(A)") ""
       write(*, "(A)") " *** af_print_info(tree) ***"
       write(*, "(A,I10)") " Current maximum level:  ", tree%highest_lvl
       write(*, "(A,I10)") " Number of boxes used:   ", af_num_boxes_used(tree)
       write(*, "(A,I10)") " Number of leaves used:  ", af_num_leaves_used(tree)
       write(*, "(A,I10)") " Memory limit(boxes):    ", tree%box_limit
       write(*, "(A,E10.2)") " Memory limit(GByte):    ", &
            tree%box_limit * 0.5_dp**30 * &
            af_box_bytes(tree%n_cell, tree%n_var_cell, tree%n_var_face)
       write(*, "(A,I10)") " Highest id in box list: ", tree%highest_id
       write(*, "(A,I10)") " Size of box list:       ", size(tree%boxes)
       write(*, "(A,I10)") " Box size (cells):       ", tree%n_cell
       write(*, "(A,I10)") " Number of cc variables: ", tree%n_var_cell
       write(*, "(A,I10)") " Number of fc variables: ", tree%n_var_face
       write(*, "(A,A15)") " Type of coordinates:    ", &
            af_coord_names(tree%coord_t)
       write(fmt_string, "(A,I0,A)") "(A,", NDIM, "E10.2)"
       write(*, fmt_string) " min. coords:            ", tree%r_base
       write(*, fmt_string) " dr at level one:        ", tree%dr_base
       write(*, "(A,E10.2)") " minval(dr):             ", af_min_dr(tree)
       write(*, "(A)") ""
    end if
  end subroutine af_print_info

  pure function af_box_bytes(n_cell, n_var_cell, n_var_face) result(box_bytes)
    integer, intent(in) :: n_cell     !< number of cells per dimension
    integer, intent(in) :: n_var_cell !< number of cell-centered variables
    integer, intent(in) :: n_var_face !< number of face-centered variables
    integer             :: box_bytes
    type(box_t)       :: dummy_box

    box_bytes = 8 * n_var_cell * (n_cell + 2)**NDIM + &
         8 * n_var_face * (n_cell + 1) * n_cell**(NDIM-1) + &
         int(storage_size(dummy_box) / 8)
  end function af_box_bytes

  pure function af_num_boxes_used(tree) result(n_boxes)
    type(af_t), intent(in) :: tree
    integer                 :: n_boxes, lvl

    n_boxes = 0
    do lvl = 1, tree%highest_lvl
       n_boxes = n_boxes + size(tree%lvls(lvl)%ids)
    end do
  end function af_num_boxes_used

  pure function af_num_leaves_used(tree) result(n_boxes)
    type(af_t), intent(in) :: tree
    integer                 :: n_boxes, lvl

    n_boxes = 0
    do lvl = 1, tree%highest_lvl
       n_boxes = n_boxes + size(tree%lvls(lvl)%leaves)
    end do
  end function af_num_leaves_used

  pure function af_num_cells_used(tree) result(n_cells)
    type(af_t), intent(in) :: tree
    integer                 :: n_cells

    n_cells = af_num_leaves_used(tree) * tree%n_cell**NDIM
  end function af_num_cells_used

  pure function af_total_volume(tree) result(vol)
    type(af_t), intent(in) :: tree
    real(dp)                :: vol
    integer                 :: i, id
    real(dp)                :: r0, r1, box_len(NDIM)
    real(dp), parameter     :: pi = acos(-1.0_dp)

    box_len = tree%n_cell * tree%dr_base

    if (NDIM == 2 .and. tree%coord_t == af_cyl) then
       vol     = 0.0_dp
       do i = 1, size(tree%lvls(1)%ids)
          id  = tree%lvls(1)%ids(i)
          r0  = tree%boxes(id)%r_min(1)
          r1  = r0 + box_len(1)
          vol = vol + pi * (r1**2 - r0**2) * box_len(NDIM)
       end do
    else
       vol = size(tree%lvls(1)%ids) * product(box_len)
    end if
  end function af_total_volume

  !> Return .true. if a box has children
  elemental logical function af_has_children(box)
    type(box_t), intent(in) :: box

    ! Boxes are either fully refined or not, so we only need to check one of the
    ! children
    af_has_children = (box%children(1) /= af_no_box)
  end function af_has_children

  !> Return .true. where there is a physical/periodic boundary. Detecting
  !> physical boundaries is straightforward (simply test whether neighbors(i) <
  !> af_no_box), but periodic boundaries require a comparison of their spatial
  !> index.
  pure function af_is_phys_boundary(boxes, id, nbs) result(boundary)
    type(box_t), intent(in) :: boxes(:) !< List of boxes
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
          dim = af_neighb_dim(nb)

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

          if (dix /= af_neighb_high_pm(nb)) then
             ! The difference in index is not the expected +-1, so a periodic
             ! boundary
             boundary(n) = .true.
          else
             boundary(n) = .false.
          end if
       end if
    end do

  end function af_is_phys_boundary

  !> Get the offset of a box with respect to its parent (e.g. in 2d, there can
  !> be a child at offset 0,0, one at n_cell/2,0, one at 0,n_cell/2 and one at
  !> n_cell/2, n_cell/2)
  pure function af_get_child_offset(box, nb) result(ix_offset)
    type(box_t), intent(in)           :: box   !< A child box
    integer, intent(in), optional      :: nb     !< Optional: get index on parent neighbor
    integer                            :: ix_offset(NDIM)

    ix_offset = iand(box%ix-1, 1) * ishft(box%n_cell, -1) ! * n_cell / 2
    if (present(nb)) ix_offset = ix_offset - af_neighb_dix(:, nb) * box%n_cell
  end function af_get_child_offset

  !> Given a cell index on box, get index of the closest cell at its parent
  pure function af_get_ix_on_parent(box, ix) result(p_ix)
    type(box_t), intent(in) :: box    !< A child box
    integer, intent(in)       :: ix(NDIM) !< Index on child box
    integer                   :: p_ix(NDIM)
    p_ix = af_get_child_offset(box) + ishft(ix+1, -1)
  end function af_get_ix_on_parent

  !> Given a cell index on box, get index on a neighbor
  pure function af_get_ix_on_neighb(box, ix, nb) result(nb_ix)
    type(box_t), intent(in) :: box    !< A box
    integer, intent(in)       :: ix(NDIM) !< Index on box
    integer, intent(in)       :: nb     !< Neighbor identifier
    integer                   :: nb_ix(NDIM), nb_dim

    nb_dim        = af_neighb_dim(nb)
    nb_ix         = ix
    nb_ix(nb_dim) = nb_ix(nb_dim) - af_neighb_high_pm(nb) * box%n_cell
  end function af_get_ix_on_neighb

  !> Given a list of neighbor directions, compute the index offset
  pure function af_neighb_offset(nbs) result(dix)
    integer, intent(in) :: nbs(:) !< List of neighbor directions
    integer             :: n, dim, dix(NDIM)

    dix = 0
    do n = 1, size(nbs)
       dim = af_neighb_dim(nbs(n))
       dix(dim) = dix(dim) + af_neighb_high_pm(nbs(n))
    end do
  end function af_neighb_offset

  !> Get index range of boundary cells inside a box facing neighbor nb
  subroutine af_get_index_bc_inside(nb, nc, n_gc, lo, hi)
    integer, intent(in)  :: nb !< Neighbor direction
    integer, intent(in)  :: nc !< box size
    integer, intent(in)  :: n_gc !< Number of ghost cells
    integer, intent(out) :: lo(NDIM)
    integer, intent(out) :: hi(NDIM)
    integer              :: nb_dim

    ! Determine index range next to boundary
    nb_dim     = af_neighb_dim(nb)
    lo(:)      = 1
    hi(:)      = nc
    if (af_neighb_low(nb)) then
       lo(nb_dim) = 1
       hi(nb_dim) = n_gc
    else
       lo(nb_dim) = nc - n_gc + 1
       hi(nb_dim) = nc
    end if
  end subroutine af_get_index_bc_inside

  !> Get index range of ghost cells facing neighbor nb
  subroutine af_get_index_bc_outside(nb, nc, n_gc, lo, hi)
    integer, intent(in)  :: nb  !< Neighbor direction
    integer, intent(in)  :: nc  !< box size
    integer, intent(in)  :: n_gc !< Number of ghost cells
    integer, intent(out) :: lo(NDIM)
    integer, intent(out) :: hi(NDIM)
    integer              :: nb_dim

    ! Determine index range next to boundary
    nb_dim     = af_neighb_dim(nb)
    lo(:)      = 1
    hi(:)      = nc
    if (af_neighb_low(nb)) then
       lo(nb_dim) = 1 - n_gc
       hi(nb_dim) = 0
    else
       lo(nb_dim) = nc + 1
       hi(nb_dim) = nc + n_gc
    end if
  end subroutine af_get_index_bc_outside

  !> Compute the 'child index' for a box with spatial index ix. With 'child
  !> index' we mean the index in the children(:) array of its parent.
  pure integer function af_ix_to_ichild(ix)
    integer, intent(in) :: ix(NDIM) !< Spatial index of the box
    ! The index can range from 1 (all ix odd) and 2**NDIM (all ix even)
#if NDIM == 1
    af_ix_to_ichild = 2 - iand(ix(1), 1)
#elif NDIM == 2
    af_ix_to_ichild = 4 - 2 * iand(ix(2), 1) - iand(ix(1), 1)
#elif NDIM == 3
    af_ix_to_ichild = &
         8 - 4 * iand(ix(3), 1) - 2 * iand(ix(2), 1) - iand(ix(1), 1)
#endif
  end function af_ix_to_ichild

  !> Get the cell index in which r lies. This routine does not check whether r
  !> is actually located inside the box.
  pure function af_cc_ix(box, r) result(cc_ix)
    type(box_t), intent(in) :: box
    real(dp), intent(in)     :: r(NDIM)
    integer                  :: cc_ix(NDIM)
    cc_ix = ceiling((r - box%r_min) / box%dr)
  end function af_cc_ix

  !> Get the location of the cell center with index cc_ix
  pure function af_r_cc(box, cc_ix) result(r)
    type(box_t), intent(in) :: box
    integer, intent(in)      :: cc_ix(NDIM)
    real(dp)                 :: r(NDIM)
    r = box%r_min + (cc_ix-0.5_dp) * box%dr
  end function af_r_cc

  !> Get the location of the face parallel to dim with index fc_ix
  pure function af_r_fc(box, dim, fc_ix) result(r)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: dim
    integer, intent(in)     :: fc_ix(NDIM)
    real(dp)                :: r(NDIM)
    r      = box%r_min + (fc_ix-0.5_dp) * box%dr
    r(dim) = r(dim) - 0.5_dp * box%dr(dim)
  end function af_r_fc

  !> Get a general location with index cc_ix (like af_r_cc but using reals)
  pure function af_rr_cc(box, cc_ix) result(r)
    type(box_t), intent(in) :: box
    real(dp), intent(in)     :: cc_ix(NDIM)
    real(dp)                 :: r(NDIM)
    r = box%r_min + (cc_ix-0.5_dp) * box%dr
  end function af_rr_cc

  !> Return the coordinate of the center of a box
  pure function af_r_center(box) result(r_center)
    type(box_t), intent(in) :: box
    real(dp)                 :: r_center(NDIM)
    r_center = box%r_min + 0.5_dp * box%n_cell * box%dr
  end function af_r_center

  !> Return finest dr that is used in the tree
  pure function af_min_dr(tree) result(dr)
    type(af_t), intent(in) :: tree
    real(dp)               :: dr !< Output: dr at the finest lvl of the tree
    dr = minval(af_lvl_dr(tree, tree%highest_lvl))
  end function af_min_dr

  !> Return dr at lvl
  pure function af_lvl_dr(tree, lvl) result(dr)
    type(af_t), intent(in) :: tree
    integer, intent(in)    :: lvl
    real(dp)               :: dr(NDIM) !< Output: dr at the finest lvl of the tree
    dr = tree%dr_base * 0.5_dp**(lvl-1)
  end function af_lvl_dr

  subroutine af_set_box_gc(box, nb, iv, gc_scalar, gc_array)
    type(box_t), intent(inout) :: box      !< Box to operate on
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    real(dp), optional, intent(in) :: gc_scalar !< Scalar value for ghost cells
#if NDIM == 1
    real(dp), optional, intent(in) :: gc_array(1) !< Array for ghost cells
#elif NDIM == 2
    real(dp), optional, intent(in) :: gc_array(box%n_cell) !< Array for ghost cells
#elif NDIM == 3
    !> Array for ghost cells
    real(dp), optional, intent(in) :: gc_array(box%n_cell, box%n_cell)
#endif
    integer                     :: nc

    nc = box%n_cell

    if (present(gc_array)) then
       select case (nb)
#if NDIM == 1
       case (af_neighb_lowx)
          box%cc(0, iv)    = gc_array(1)
       case (af_neighb_highx)
          box%cc(nc+1, iv) = gc_array(1)
#elif NDIM == 2
       case (af_neighb_lowx)
          box%cc(0, 1:nc, iv)    = gc_array
       case (af_neighb_highx)
          box%cc(nc+1, 1:nc, iv) = gc_array
       case (af_neighb_lowy)
          box%cc(1:nc, 0, iv)    = gc_array
       case (af_neighb_highy)
          box%cc(1:nc, nc+1, iv) = gc_array
#elif NDIM == 3
       case (af_neighb_lowx)
          box%cc(0, 1:nc, 1:nc, iv)    = gc_array
       case (af_neighb_highx)
          box%cc(nc+1, 1:nc, 1:nc, iv) = gc_array
       case (af_neighb_lowy)
          box%cc(1:nc, 0, 1:nc, iv)    = gc_array
       case (af_neighb_highy)
          box%cc(1:nc, nc+1, 1:nc, iv) = gc_array
       case (af_neighb_lowz)
          box%cc(1:nc, 1:nc, 0, iv)    = gc_array
       case (af_neighb_highz)
          box%cc(1:nc, 1:nc, nc+1, iv) = gc_array
#endif
       end select
    else if (present(gc_scalar)) then
       select case (nb)
#if NDIM == 1
       case (af_neighb_lowx)
          box%cc(0, iv)    = gc_scalar
       case (af_neighb_highx)
          box%cc(nc+1, iv) = gc_scalar
#elif NDIM == 2
       case (af_neighb_lowx)
          box%cc(0, 1:nc, iv)    = gc_scalar
       case (af_neighb_highx)
          box%cc(nc+1, 1:nc, iv) = gc_scalar
       case (af_neighb_lowy)
          box%cc(1:nc, 0, iv)    = gc_scalar
       case (af_neighb_highy)
          box%cc(1:nc, nc+1, iv) = gc_scalar
#elif NDIM == 3
       case (af_neighb_lowx)
          box%cc(0, 1:nc, 1:nc, iv)    = gc_scalar
       case (af_neighb_highx)
          box%cc(nc+1, 1:nc, 1:nc, iv) = gc_scalar
       case (af_neighb_lowy)
          box%cc(1:nc, 0, 1:nc, iv)    = gc_scalar
       case (af_neighb_highy)
          box%cc(1:nc, nc+1, 1:nc, iv) = gc_scalar
       case (af_neighb_lowz)
          box%cc(1:nc, 1:nc, 0, iv)    = gc_scalar
       case (af_neighb_highz)
          box%cc(1:nc, 1:nc, nc+1, iv) = gc_scalar
#endif
       end select
    else
       stop "af_set_box_gc: requires gc_scalar or gc_array"
    end if
  end subroutine af_set_box_gc

#if NDIM == 2
  !> Get the radius of the cell center with first index i
  pure function af_cyl_radius_cc(box, i) result(r)
    type(box_t), intent(in) :: box
    integer, intent(in)      :: i
    real(dp)                 :: r
    r = box%r_min(1) + (i-0.5_dp) * box%dr(1)
  end function af_cyl_radius_cc

  !> Get the normalized weights of the 'inner' and 'outer' children of a cell
  !> with index ix. Note that the cell centers of the children are located at
  !> -/+ 0.25 dr compared to the parent.
  subroutine af_cyl_child_weights(box, i, inner, outer)
    type(box_t), intent(in) :: box
    integer, intent(in)      :: i
    real(dp), intent(out)    :: inner, outer
    real(dp)                 :: tmp

    tmp = 0.25_dp * box%dr(1) / af_cyl_radius_cc(box, i)
    inner = 1 - tmp
    outer = 1 + tmp
  end subroutine af_cyl_child_weights

  !> Get the factors for the left and right flux in each cell
  subroutine af_cyl_flux_factors(box, flux_factors)
    type(box_t), intent(in) :: box
    real(dp), intent(out)   :: flux_factors(2, box%n_cell)
    integer                 :: i
    real(dp)                :: r, inv_r

    do i = 1, box%n_cell
       r                  = af_cyl_radius_cc(box, i)
       inv_r              = 1/r
       flux_factors(1, i) = (r - 0.5_dp * box%dr(1)) * inv_r
       flux_factors(2, i) = (r + 0.5_dp * box%dr(1)) * inv_r
    end do
  end subroutine af_cyl_flux_factors
#endif

end module m_af_types
