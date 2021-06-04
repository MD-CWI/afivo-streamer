#include "../src/cpp_macros.h"
!> \example particles_to_grid.f90
!>
!> Example showing how to map particles to a grid
program particles_to_grid
  use m_af_all

  implicit none

  integer, parameter    :: box_size    = 8
  integer               :: i_phi
  integer, parameter    :: n_particles = 100*1000
  real(dp), parameter   :: domain_len  = 1.0_dp
  real(dp), parameter   :: r_min(NDIM) = -0.5_dp * domain_len
  real(dp), allocatable :: coordinates(:, :), weights(:)

  type(af_t)      :: tree
  type(ref_info_t) :: refine_info
  integer          :: n
  integer          :: count_rate, t_start, t_end
  real(dp)         :: sum_density

  print *, "Running particles_to_grid_" // DIMNAME // ""
  print *, "Number of threads", af_get_max_threads()

  call af_add_cc_variable(tree, "phi", ix=i_phi)

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [DTIMES(domain_len)], &
       [DTIMES(box_size)], &
       r_min=r_min)

  call af_set_cc_methods(tree, i_phi, af_bc_neumann_zero, af_gc_interp)

  do
     call af_adjust_refinement(tree, refine_routine, refine_info, 0)
     if (refine_info%n_add == 0) exit
  end do

  allocate(coordinates(NDIM, n_particles))
  allocate(weights(n_particles))

  call random_number(coordinates(:, :))
  weights(:) = 1.0_dp

  do n = 1, n_particles
     coordinates(:, n) = (coordinates(:, n) - 0.5_dp) * domain_len
     ! coordinates(:, n) = (sqrt(coordinates(:, n)) - 0.5_dp) * domain_len
  end do

  call system_clock(t_start, count_rate)
  call af_particles_to_grid(tree, i_phi, coordinates, weights, n_particles, 0)
  call system_clock(t_end, count_rate)
  call af_write_silo(tree, "output/particles_to_grid_" // DIMNAME // "_0", 1, 0.0_dp)
  call af_tree_sum_cc(tree, i_phi, sum_density)
  print *, "zeroth order: ", (t_end-t_start) / real(count_rate, dp), " seconds"
  print *, "integrated density: ", sum_density

  call af_tree_clear_cc(tree, i_phi)

  call system_clock(t_start, count_rate)
  call af_particles_to_grid(tree, i_phi, coordinates, weights, n_particles, 1)
  call system_clock(t_end, count_rate)
  call af_write_silo(tree, "output/particles_to_grid_" // DIMNAME // "_1", 1, 0.0_dp)
  call af_tree_sum_cc(tree, i_phi, sum_density)
  print *, "first order: ", (t_end-t_start) / real(count_rate, dp), " seconds"
  print *, "integrated density: ", sum_density

contains

  subroutine refine_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)      :: cell_flags(DTIMES(box%n_cell))

    if (all(box%r_min < 0.0_dp) .and. box%lvl <= 5) then
       cell_flags(DTIMES(:)) = af_do_ref
    else
       cell_flags(DTIMES(:)) = af_keep_ref
    end if
  end subroutine refine_routine

end program particles_to_grid
