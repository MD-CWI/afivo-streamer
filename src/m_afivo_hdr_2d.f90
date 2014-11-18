module m_afivo_hdr_2d
  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

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

  type box2_t
     integer               :: lvl          ! level of the box
     integer               :: tag          ! for setting tag bits
     integer               :: ix(2)        ! index in the domain
     integer               :: parent       ! index of parent in box list
     integer               :: children(4)  ! index of children in box list
     integer               :: neighbors(4) ! index of neighbors in box list
     integer               :: n_cell       ! number of cells per dimension
     real(dp)              :: dr           ! width/height of a cell
     real(dp)              :: r_min(2)     ! min coords. of box
     real(dp), allocatable :: cc(:, :, :)  ! cell centered variables
     real(dp), allocatable :: fx(:, :, :)  ! x-face centered variables
     real(dp), allocatable :: fy(:, :, :)  ! y-face centered variables
  end type box2_t

  type lvl2_t
     integer, allocatable      :: ids(:)       ! indices of boxes of level
     integer, allocatable      :: leaves(:)     ! all ids(:) that are leaves
     integer, allocatable      :: parents(:)   ! all ids(:) that have children
  end type lvl2_t

  type a2_t
     integer                   :: max_lvl    ! maximum allowed level
     integer                   :: n_lvls     ! current maximum level
     integer                   :: max_id     ! max index in box list
     integer                   :: n_cell     ! number of cells per dimension
     integer                   :: n_var_cell ! number of cc variables
     integer                   :: n_var_face ! number of fc variables
     real(dp)                  :: r_base(2)  ! coords of box at index (1,1)
     real(dp)                  :: dr_base    ! cell spacing at lvl 1
     type(lvl2_t), allocatable :: lvls(:)    ! list storing the tree levels
     type(box2_t), allocatable :: boxes(:)   ! list of all boxes
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

     subroutine a2_subr_arg(box, rarg)
       import
       type(box2_t), intent(inout) :: box
       real(dp), intent(in)        :: rarg(:)
     end subroutine a2_subr_arg

     subroutine a2_subr_boxes(boxes, id)
       import
       type(box2_t), intent(inout) :: boxes(:)
       integer, intent(in)         :: id
     end subroutine a2_subr_boxes

     subroutine a2_subr_boxes_arg(boxes, id, rarg)
       import
       type(box2_t), intent(inout) :: boxes(:)
       integer, intent(in)         :: id
       real(dp), intent(in)        :: rarg(:)
     end subroutine a2_subr_boxes_arg

     subroutine a2_subr_gc(boxes, id, i, iv)
       import
       type(box2_t), intent(inout) :: boxes(:)
       integer, intent(in)         :: id, i, iv
     end subroutine a2_subr_gc
  end interface
end module m_afivo_hdr_2d