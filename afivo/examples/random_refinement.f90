#include "../src/cpp_macros.h"
!> \example random_refinement_Xd.f90
!>
!> This example shows how to create an AMR tree, perform random refinement,
!> write output files, and how to fill ghost cells.
program random_refinement_Xd
  use m_af_all

  implicit none

  type(af_t)         :: tree
  integer            :: grid_size(NDIM)
  real(dp)           :: domain_size(NDIM)
  logical            :: periodic(NDIM) = .true.
  integer            :: iter, boxes_used
  integer, parameter :: coord_type     = af_xyz ! af_xyz or af_cyl

  integer, parameter :: box_size   = 8
  integer, parameter :: i_phi      = 1
  type(ref_info_t)   :: ref_info
  real(dp)           :: sum_phi_t0, sum_phi
  character(len=100) :: fname
  integer            :: count_rate,t_start,t_end

  write(*,'(A,I0,A)') 'program random_refinement_', NDIM, "d"

  print *, "Number of threads", af_get_max_threads()

  call af_add_cc_variable(tree, "phi")

  call af_set_cc_methods(tree, 1, af_bc_dirichlet_zero, &
       prolong=af_prolong_linear)

  ! Initialize tree
  grid_size(1) = 2 * box_size
  grid_size(2:) = box_size
  domain_size(1) = 2 * acos(-1.0_dp)
  domain_size(2:) = acos(-1.0_dp)

  if (command_argument_count() == 1) then
     ! Load tree from file
     call get_command_argument(1, fname)
     call af_read_tree(tree, fname, read_other_data=read_string)
  else
     call af_init(tree, & ! Tree to initialize
          box_size, &     ! A box contains box_size**DIM cells
          domain_size, &
          grid_size, &
          periodic=periodic, &
          coord=coord_type)

     call af_print_info(tree)

     ! Set variables on base by using the helper functions af_loop_box(tree, sub)
     ! and af_loop_boxes(tree, sub). These functions call the subroutine sub for
     ! each box in the tree, with a slightly different syntax.
     call af_loop_box(tree, set_init_cond)

     ! Fill ghost cells for phi. The third argument is a subroutine that fills
     ! ghost cells near refinement boundaries, and the fourth argument fill ghost
     ! cells near physical boundaries.
     call af_gc_tree(tree, i_phi, af_gc_interp, af_bc_dirichlet_zero)
  end if

  call af_tree_sum_cc(tree, i_phi, sum_phi_t0)

  call system_clock(t_start, count_rate)
  boxes_used = 1
  do iter = 1, 15
     ! This writes a VTK output file containing the cell-centered values of the
     ! leaves of the tree (the boxes not covered by refinement).
     ! Variables are the names given as the third argument.
     write(fname, "(A,I0)") "output/random_refinement_" // DIMNAME // "_", iter
     call af_write_silo(tree, trim(fname), n_cycle=iter)
     call af_write_tree(tree, fname, write_other_data=write_string)

     ! This updates the refinement of the tree, by at most one level per call.
     ! The second argument is a subroutine that is called for each box that can
     ! be refined or derefined, and it should set refinement flags. Information
     ! about the changes in refinement are returned in the third argument.
     !
     ! Within each (de)refinement step the subroutine af_tidy_up is called. This
     ! means that every now and then whether too much holes appear in the tree
     ! box list. Therefore reordering of tree is not necessary at the end of the
     ! iteration process
     call af_adjust_refinement(tree, ref_routine, ref_info, ref_buffer=0)

     call af_tree_sum_cc(tree, i_phi, sum_phi)
     write(*, "(A,E10.2)") " conservation error: ", sum_phi - sum_phi_t0
     boxes_used = boxes_used + ref_info%n_add - ref_info%n_rm
     write(*,'(4(3x,A,1x,i6))') "# new     boxes", ref_info%n_add, &
          "# removed boxes", ref_info%n_rm,  &
          "# boxes used   ", boxes_used,     &
          " highest level ", tree%highest_lvl
  end do
  call system_clock(t_end, count_rate)

  write(*, '(A,i3,1x,A,f8.2,1x,A,/)') &
       ' Wall-clock time after ',iter, &
       ' iterations: ', (t_end-t_start) / real(count_rate, dp), &
       ' seconds'

  call af_print_info(tree)

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call af_destroy(tree)

contains

  ! Return the refinement flag for boxes(id)
  subroutine ref_routine(box, cell_flags)
    type(box_t), intent(in) :: box ! A list of all boxes in the tree
    integer, intent(out)     :: cell_flags(DTIMES(box%n_cell))
    real(dp)                 :: rr

    ! Draw a [0, 1) random number
    call random_number(rr)

    if (rr < 0.5_dp**NDIM .and. box%lvl < 5) then
       cell_flags = af_do_ref   ! Add refinement
    else
       cell_flags = af_rm_ref   ! Ask to remove this box, which will not always
       ! happen (see documentation)
    end if
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

       ! Set the values at each cell according to some function
       box%cc(IJK, i_phi) = sin(0.5_dp * rr(1))**2 * cos(rr(2))**2
    end do; CLOSE_DO
  end subroutine set_init_cond

  subroutine write_string(my_unit)
    integer, intent(in) :: my_unit
    character(len=10)   :: str

    str = "Hello"
    write(my_unit) str
  end subroutine write_string

  subroutine read_string(my_unit)
    integer, intent(in) :: my_unit
    character(len=10)   :: str
    read(my_unit) str

    if (str /= "Hello") error stop "error reading other data"
  end subroutine read_string

end program random_refinement_Xd
