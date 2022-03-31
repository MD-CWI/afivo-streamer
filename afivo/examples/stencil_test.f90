#include "../src/cpp_macros.h"
!> \example stencil_test
!>
!> This example shows how to define and apply numerical stencils
program stencil_test
  use m_af_all

  implicit none

  integer, parameter  :: box_size   = 8
  real(dp), parameter :: domain_len = 2 * acos(-1.0_dp)
  integer             :: i_phi
  integer             :: i_phi_old
  type(af_t)          :: tree
  type(ref_info_t)    :: refine_info
  integer             :: it
  character(len=100)  :: fname

  integer :: stencil_key = 1

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "phi_old", ix=i_phi_old)
  call af_set_cc_methods(tree, i_phi, af_bc_neumann_zero)
  call af_set_cc_methods(tree, i_phi_old, af_bc_neumann_zero)

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [DTIMES(domain_len)], &
       [DTIMES(box_size)], &
       periodic=[DTIMES(.true.)])

  do
     call af_loop_box(tree, set_initial_condition)
     call af_adjust_refinement(tree, refine_routine, refine_info)
     if (refine_info%n_add == 0) exit
  end do

  call af_stencil_store(tree, stencil_key, set_stencil)

  do it = 1, 10
     write(fname, "(A,I0)") "output/stencil_test_" // DIMNAME // "_", it

     ! Write the cell centered data of tree to a vtk unstructured file fname.
     ! Only the leaves of the tree are used
     call af_write_silo(tree, trim(fname), it)
     call af_tree_copy_cc(tree, i_phi, i_phi_old)
     call af_stencil_apply(tree, stencil_key, i_phi_old, i_phi)
     call af_gc_tree(tree, [i_phi])
  end do

contains

  !> Return the refinement flag for box
  subroutine refine_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))

    if (maxval(box%dr) > 2e-2_dp * domain_len) then
       cell_flags = af_do_ref
    else
       cell_flags = af_keep_ref
    end if
  end subroutine refine_routine

  !> This routine sets the initial conditions for each box
  subroutine set_initial_condition(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM)

    nc = box%n_cell

    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, i_phi) = product(sin(rr))
    end do; CLOSE_DO
  end subroutine set_initial_condition

  subroutine set_stencil(box, stencil)
    type(box_t), intent(in)        :: box
    type(stencil_t), intent(inout) :: stencil

    stencil%shape = af_stencil_357
    stencil%stype = stencil_constant
    allocate(stencil%c(2*NDIM+1))

    stencil%c(1)  = 0.0_dp
    stencil%c(2:) = 1.0_dp / (2 * NDIM)
  end subroutine set_stencil

end program stencil_test
