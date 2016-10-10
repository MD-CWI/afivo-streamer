!> This module contains the basic dimension-independent types and constants that
!> are used in Afivo, together with some basic routines. Dimension-dependent
!> types and constant are place in m_a2_types and m_a3_types.
module m_afivo_types

  implicit none
  public

  ! dp stands for double precision (8 byte reals)
  integer, parameter :: dp = kind(0.0d0)

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

  !> Each box contains a tag, for which bits can be set. This is the initial
  !> value, which should not be used by the user
  integer, parameter :: af_init_tag = -huge(1)

  !> Default coordinate system
  integer, parameter :: af_xyz = 0

  !> Cylindrical coordinate system
  integer, parameter :: af_cyl = 1

  !> Value to indicate a Dirichlet boundary condition
  integer, parameter :: af_bc_dirichlet = -1

  !> Value to indicate a Neumann boundary condition
  integer, parameter :: af_bc_neumann = -2

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
     integer, allocatable :: rm(:) !< Id's of removed boxes
  end type ref_lvl_t

  !> Type that contains the refinement changes in a tree
  type ref_info_t
     integer :: n_add = 0                    !< Total number of added boxes
     integer :: n_rm = 0                     !< Total number removed boxes
     type(ref_lvl_t), allocatable :: lvls(:) !< Information per level
  end type ref_info_t

contains

  !> Get number of threads
  integer function af_get_max_threads()
    use omp_lib
    af_get_max_threads = OMP_get_max_threads()
  end function af_get_max_threads

end module m_afivo_types
