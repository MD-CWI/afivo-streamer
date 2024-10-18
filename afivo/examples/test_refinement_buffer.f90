#include "../src/cpp_macros.h"
!> \example test_refinement_buffer.f90
!>
!> This example tests the refinement buffer option, which extends the
!> refinement a given number of cells
program random_refinement
  use m_af_all

  implicit none

  type(af_t)         :: tree
  integer            :: grid_size(NDIM)
  real(dp)           :: domain_size(NDIM)
  logical            :: periodic(NDIM) = .true.
  integer            :: iter
  integer, parameter :: coord_type = af_xyz ! af_xyz or af_cyl

  integer, parameter :: box_size = 8
  integer            :: i_phi
  type(ref_info_t)   :: ref_info
  character(len=100) :: fname

  write(*,'(A,I0,A)') 'program test_refinement_buffer_', NDIM, "d"

  print *, "Number of threads", af_get_max_threads()

  call af_add_cc_variable(tree, "phi", ix=i_phi)

  call af_set_cc_methods(tree, 1, af_bc_dirichlet_zero, &
       prolong=af_prolong_linear)

  ! Initialize tree
  grid_size(:) = box_size
  domain_size(:) = 1.0_dp

  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       domain_size, &
       grid_size, &
       periodic=periodic, &
       coord=coord_type)

  do iter = 1, 10
     call af_adjust_refinement(tree, ref_routine, ref_info, ref_buffer=0)
     if (ref_info%n_add == 0) exit
  end do

  call af_loop_box(tree, set_init_cond)
  call af_gc_tree(tree, [i_phi])

  do iter = 0, box_size
     call af_adjust_refinement(tree, ref_routine, ref_info, ref_buffer=iter)

     call af_loop_box(tree, set_init_cond)
     call af_gc_tree(tree, [i_phi])

     write(fname, "(A,I0)") "output/test_refinement_buffer_" // DIMNAME // "_", iter
     call af_write_silo(tree, trim(fname), n_cycle=iter)
  end do

contains

  ! Return the refinement flag for boxes(id)
  subroutine ref_routine(box, cell_flags)
    type(box_t), intent(in) :: box ! A list of all boxes in the tree
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))
    real(dp)                :: rr(NDIM)
    integer                 :: nc, IJK

    nc = box%n_cell

    do KJI_DO(1,nc)
       rr = af_r_cc(box, [IJK])
       if (norm2(rr - 0.5_dp) < 0.25_dp .and. box%lvl < 5) then
          cell_flags(IJK) = af_do_ref   ! Add refinement
       else
          cell_flags(IJK) = af_keep_ref
       end if
    end do; CLOSE_DO

  end subroutine ref_routine

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box_t), intent(inout) :: box
    integer                     :: IJK, nc
    real(dp)                    :: rr(NDIM)

    nc = box%n_cell
    do KJI_DO(1,nc)
       ! Get the coordinate of the cell center at i,j
       rr = af_r_cc(box, [IJK])

       if (norm2(rr - 0.5_dp) < 0.25_dp) then
          box%cc(IJK, i_phi) = 1
       else
          box%cc(IJK, i_phi) = 0
       end if

    end do; CLOSE_DO
  end subroutine set_init_cond

end program random_refinement
