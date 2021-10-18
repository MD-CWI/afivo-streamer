#include "../src/cpp_macros.h"
!> Module containing basic types used in the multigrid routines
module m_mg_types
  use m_af_types
  use iso_c_binding, only: c_ptr

  implicit none
  public

  ! The mg module supports different multigrid operators, and uses these tags to
  ! identify boxes / operators
  integer, parameter :: mg_auto_operator = 1 !< Automatic
  integer, parameter :: mg_normal_box = 2 !< Normal box
  integer, parameter :: mg_lsf_box    = 3 !< Box with an internal boundary
  integer, parameter :: mg_ceps_box   = 4 !< Box with constant eps /= 1
  integer, parameter :: mg_veps_box   = 5 !< Box with varying eps (on face)

  integer, parameter :: mg_auto_prolongation = 32
  integer, parameter :: mg_prolong_linear = 33
  integer, parameter :: mg_prolong_sparse = 34
  integer, parameter :: mg_prolong_eps = 35
  integer, parameter :: mg_prolong_lsf = 36

  ! Labels for the different steps of a multigrid cycle
  integer, parameter :: mg_cycle_down = 1
  integer, parameter :: mg_cycle_up   = 3

  !> Generic type for the coarse grid solver
  type coarse_solve_t
     type(c_ptr) :: matrix
     type(c_ptr) :: rhs
     type(c_ptr) :: phi
     type(c_ptr) :: solver
     type(c_ptr) :: grid

     !> Stores coefficient to convert boundary conditions to the right-hand side
     real(dp), allocatable :: bc_to_rhs(:, :, :)

     !> Stores coefficients to use with level set function
     real(dp), allocatable :: lsf_fac(DTIMES(:), :)

     integer  :: symmetric      = 1
     integer  :: solver_type    = -1
     integer  :: max_iterations = 50
     integer  :: n_cycle_down   = 1
     integer  :: n_cycle_up     = 1
     real(dp) :: tolerance      = 1e-6_dp
  end type coarse_solve_t

  !> Type to store multigrid options in
  type :: mg_t
     integer :: i_phi        = -1 !< Variable holding solution
     integer :: i_rhs        = -1 !< Variable holding right-hand side
     integer :: i_tmp        = -1 !< Internal variable (holding prev. solution)

     integer :: i_eps        = -1 !< Optional variable (diel. permittivity)
     integer :: i_lsf        = -1 !< Optional variable for level set function

     integer :: n_cycle_down = -1 !< Number of relaxation cycles in downward sweep
     integer :: n_cycle_up   = -1 !< Number of relaxation cycles in upward sweep

     logical :: initialized = .false. !< Whether the structure has been initialized
     logical :: use_corners = .false. !< Does the smoother use corner ghost cells
     logical :: subtract_mean = .false. !< Whether to subtract mean from solution

     !> Store lambda^2 for Helmholtz equations (L phi - lamda phi = f)
     real(dp) :: helmholtz_lambda = 0.0_dp

     !> Boundary value for level set function
     real(dp) :: lsf_boundary_value = 0.0_dp

     !> Maximal gradient of the level set function
     real(dp) :: lsf_max_gradient = 0.0_dp

     !> Level-set function
     procedure(mg_func_lsf), pointer, nopass :: lsf => null()

     !> Routine to determine distance from level-set function
     procedure(mg_lsf_distf), pointer, nopass :: lsf_dist => null()

     !> Routine to call for filling ghost cells near physical boundaries
     procedure(af_subr_bc), pointer, nopass   :: sides_bc => null()

     !> Routine to call for filling ghost cells near refinement boundaries
     procedure(af_subr_rb), pointer, nopass   :: sides_rb => null()

     !> Subroutine that performs the (non)linear operator
     procedure(mg_box_op), pointer, nopass   :: box_op => null()

     !> Key indicating which stencil is to be used for the operator
     integer :: operator_key = af_stencil_none

     !> Key indicating which stencil is to be used for the operator
     integer :: prolongation_key = af_stencil_none

     !> Subroutine that performs Gauss-Seidel relaxation on a box
     procedure(mg_box_gsrb), pointer, nopass :: box_gsrb => null()

     !> Subroutine that corrects the children of a box
     procedure(mg_box_corr), pointer, nopass :: box_corr => null()

     !> Subroutine for restriction
     procedure(mg_box_rstr), pointer, nopass :: box_rstr => null()

     !> Subroutine for getting the stencil
     procedure(mg_box_stencil), pointer, nopass :: box_stencil => null()

     type(coarse_solve_t) :: csolver
  end type mg_t

  abstract interface
     !> Subroutine that performs A * cc(..., i_in) = cc(..., i_out)
     subroutine mg_box_op(box, i_out, mg)
       import
       type(box_t), intent(inout) :: box   !< The box to operate on
       type(mg_t), intent(in)     :: mg    !< Multigrid options
       integer, intent(in)        :: i_out !< Index of output variable
     end subroutine mg_box_op

     !> Subroutine that performs Gauss-Seidel relaxation
     subroutine mg_box_gsrb(box, redblack_cntr, mg)
       import
       type(box_t), intent(inout) :: box           !< The box to operate on
       type(mg_t), intent(in)     :: mg            !< Multigrid options
       integer, intent(in)        :: redblack_cntr !< Iteration counter
     end subroutine mg_box_gsrb

     subroutine mg_box_corr(box_p, box_c, mg)
       import
       type(box_t), intent(inout) :: box_c
       type(box_t), intent(in)    :: box_p
       type(mg_t), intent(in)     :: mg !< Multigrid options
     end subroutine mg_box_corr

     subroutine mg_box_rstr(box_c, box_p, iv, mg)
       import
       type(box_t), intent(in)    :: box_c !< Child box to restrict
       type(box_t), intent(inout) :: box_p !< Parent box to restrict to
       integer, intent(in)        :: iv    !< Variable to restrict
       type(mg_t), intent(in)     :: mg    !< Multigrid options
     end subroutine mg_box_rstr

     subroutine mg_box_stencil(box, mg, stencil, bc_to_rhs, lsf_fac)
       import
       type(box_t), intent(in) :: box
       type(mg_t), intent(in)  :: mg
       real(dp), intent(inout) :: stencil(2*NDIM+1, DTIMES(box%n_cell))
       real(dp), intent(inout) :: bc_to_rhs(box%n_cell**(NDIM-1), af_num_neighbors)
       real(dp), intent(inout) :: lsf_fac(DTIMES(box%n_cell))
     end subroutine mg_box_stencil

     !> Level set function
     real(dp) function mg_func_lsf(rr)
       import
       real(dp), intent(in) :: rr(NDIM) !< Coordinates
     end function mg_func_lsf

     !> Compute distance to boundary starting at point a going to point b, in
     !> the range from [0, 1], with 1 meaning there is no boundary
     real(dp) function mg_lsf_distf(box_a, IJK_(a), box_b, IJK_(b), mg)
       import
       type(box_t), intent(in) :: box_a   !< Box a (start point)
       integer, intent(in)     :: IJK_(a) !< Cell-centered index in box a
       type(box_t), intent(in) :: box_b   !< Box b (end point)
       integer, intent(in)     :: IJK_(b) !< Cell-centered index in box b
       type(mg_t), intent(in)  :: mg
     end function mg_lsf_distf
  end interface

end module m_mg_types
