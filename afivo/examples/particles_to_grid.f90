#include "../src/cpp_macros.h"
!> \example particles_to_grid.f90
!>
!> Example showing how to map particles to a grid
program particles_to_grid
  use m_af_all

  implicit none

  integer, parameter    :: box_size    = 8
  integer               :: i_phi, i_tmp
  integer, parameter    :: n_particles = 20*1000*1000
  real(dp), parameter   :: r_max(NDIM) = 1.0_dp
  real(dp), allocatable :: coordinates(:, :), weights(:)
  integer               :: coordinate_system = af_xyz
  logical               :: cylindrical = .false.
  logical               :: use_tmp_var = .false.

  type(af_t)      :: tree
  type(ref_info_t) :: refine_info
  integer          :: n
  integer          :: count_rate, t_start, t_end
  real(dp)         :: sum_density
  character(len=10) :: tmp_string

  if (command_argument_count() > 0) then
     call get_command_argument(1, tmp_string)
     read(tmp_string, *) cylindrical
  end if

  if (command_argument_count() > 1) then
     call get_command_argument(2, tmp_string)
     read(tmp_string, *) use_tmp_var
  end if

  if (cylindrical) coordinate_system = af_cyl

  print *, "Running particles_to_grid_" // DIMNAME // ""
  print *, "Number of threads", af_get_max_threads()
  print *, "Cylindrical", cylindrical
  print *, "Number of particles", n_particles

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)

  call af_init(tree, box_size, r_max, [DTIMES(box_size)], &
       coord=coordinate_system)

  call af_set_cc_methods(tree, i_phi, af_bc_neumann_zero, af_gc_interp)

  do
     call af_adjust_refinement(tree, refine_routine, refine_info, 0)
     if (refine_info%n_add == 0) exit
  end do

  allocate(coordinates(NDIM, n_particles))
  allocate(weights(n_particles))

  call random_number(coordinates(:, :))
  weights(:) = 1.0_dp / n_particles

  if (cylindrical) then
     coordinates(1, :) = sqrt(coordinates(1, :))
  end if

  do n = 1, n_particles
     coordinates(:, n) = 0.6_dp * coordinates(:, n) * r_max
  end do

  call system_clock(t_start, count_rate)

  if (use_tmp_var) then
     call af_particles_to_grid(tree, i_phi, n_particles, &
          get_id, get_rw, 0, iv_tmp=i_tmp)
  else
     call af_particles_to_grid(tree, i_phi, n_particles, &
          get_id, get_rw, 0)
  end if

  call system_clock(t_end, count_rate)
  call af_write_silo(tree, "output/particles_to_grid_" // DIMNAME // "_0", 1, 0.0_dp)
  call af_tree_sum_cc(tree, i_phi, sum_density)
  print *, "zeroth order: ", (t_end-t_start) / real(count_rate, dp), " seconds"
  print *, "integrated density: ", sum_density

  call af_tree_clear_cc(tree, i_phi)

  call system_clock(t_start, count_rate)

  if (use_tmp_var) then
     call af_particles_to_grid(tree, i_phi, n_particles, &
          get_id, get_rw, 1, iv_tmp=i_tmp)
  else
     call af_particles_to_grid(tree, i_phi, n_particles, &
          get_id, get_rw, 1)
  end if

  call system_clock(t_end, count_rate)
  call af_write_silo(tree, "output/particles_to_grid_" // DIMNAME // "_1", 1, 0.0_dp)
  call af_tree_sum_cc(tree, i_phi, sum_density)
  print *, "first order: ", (t_end-t_start) / real(count_rate, dp), " seconds"
  print *, "integrated density: ", sum_density

contains

  subroutine refine_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)      :: cell_flags(DTIMES(box%n_cell))

    if (all(box%r_min < 0.5_dp * r_max) .and. box%lvl <= 5) then
       cell_flags(DTIMES(:)) = af_do_ref
    else
       cell_flags(DTIMES(:)) = af_keep_ref
    end if
  end subroutine refine_routine

  subroutine get_id(n, id)
    integer, intent(in)  :: n
    integer, intent(out) :: id

    id = af_get_id_at(tree, coordinates(:, n))
  end subroutine get_id

  subroutine get_rw(n, r, w)
    integer, intent(in)   :: n
    real(dp), intent(out) :: r(NDIM)
    real(dp), intent(out) :: w

    r = coordinates(:, n)
    w = weights(n)
  end subroutine get_rw

end program particles_to_grid
