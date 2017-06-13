#include "../src/cpp_macros_$Dd.h"
!> \example particles_to_grid_$Dd.f90
!>
!> Example showing how to map particles to a grid
program particles_to_grid_$Dd
  use m_a$D_all

  implicit none

  integer, parameter  :: box_size   = 8
  integer, parameter  :: i_phi      = 1
  integer, parameter :: n_particles = 1000*1000
  real(dp), parameter :: domain_len = 2.0_dp
  real(dp), parameter :: r_min($D) = -0.5_dp * domain_len
  real(dp), parameter :: dr         = domain_len / box_size

  type(a$D_t)        :: tree
  type(ref_info_t)   :: refine_info
  integer            :: id, ix_list($D, 1)

  print *, "Running particles_to_grid_$Dd"
  print *, "Number of threads", af_get_max_threads()

  ! Initialize tree
  call a$D_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       n_var_cell=1, & ! Number of cell-centered variables
       n_var_face=0, & ! Number of face-centered variables
       dr=dr, &        ! Distance between cells on base level
       cc_names=["phi"], & ! Variable names
       r_min=r_min)

  ! Set up geometry
  id             = 1
  ix_list(:, id) = 1            ! Set index of box

  ! Create the base mesh, using the box indices and their neighbor information
  call a$D_set_base(tree, 1, ix_list)

  do
    call a$D_adjust_refinement(tree, refine_routine, refine_info, 0)
    if (refine_info%n_add == 0) exit
  end do

  call add_particles(tree)
  call a$D_gc_tree(tree, i_phi, a$D_gc_interp, a$D_bc_neumann_zero)
  call a$D_write_silo(tree, "particles_to_grid", 1, 0.0_dp, dir="output")

contains

  subroutine refine_routine(box, cell_flags)
    type(box$D_t), intent(in) :: box
    integer, intent(out)      :: cell_flags(DTIMES(box%n_cell))

    if (box%r_min(1) < 0.0_dp .and. box%lvl <= 5) then
      cell_flags(DTIMES(:)) = af_do_ref
    else
      cell_flags(DTIMES(:)) = af_keep_ref
    end if
  end subroutine refine_routine

  subroutine add_particles(tree)
    type(a$D_t), intent(inout) :: tree
    integer                      :: n
    real(dp)                     :: rr($D)

    do n = 1, n_particles
      call random_number(rr)
      rr = (rr - 0.5_dp) * domain_len
      if (norm2(rr) < 1.0_dp) then
        call a$D_interp0_to_grid(tree, rr, i_phi, 1.0_dp, .true.)
      end if
    end do
  end subroutine add_particles

end program particles_to_grid_$Dd
