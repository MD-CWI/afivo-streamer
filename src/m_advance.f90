#include "../afivo/src/cpp_macros.h"
!> Module for advancing solution in time
module m_advance
  use m_af_all
  use m_dt
  use m_advance_base

  implicit none
  private

  !> Maximal allowed time step
  real(dp), public, protected :: advance_max_dt

  ! Public methods
  public :: advance_set_max_dt
  public :: advance
  public :: advance_num_states

contains

  subroutine advance_set_max_dt(dt_max)
    real(dp), intent(in) :: dt_max
    advance_max_dt = dt_max
  end subroutine advance_set_max_dt

  subroutine advance(tree, dt, time)
    use m_chemistry
    use m_fluid_lfa
    use m_field
    use m_streamer
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt
    real(dp), intent(inout)   :: time

    dt_matrix(1:dt_num_cond, :) = dt_max ! Maximum time step

    select case (time_integrator)
    case (forward_euler_t)
       call forward_euler(tree, dt, time, 0, 0, 0, .true., 1)
       call af_restrict_ref_boundary(tree, flux_species)
       time = time + dt
    case (rk2_t)
       call forward_euler(tree, 0.5_dp * dt, time, 0, 0, 1, .false., 1)
       call af_restrict_ref_boundary(tree, flux_species+1)
       call forward_euler(tree, dt, time + 0.5_dp*dt, 1, 0, 0, .true., 2)
       call af_restrict_ref_boundary(tree, flux_species)
       time = time + dt
    case (heuns_method_t)
       call forward_euler(tree, dt, time, 0, 0, 1, .false., 1)
       call af_restrict_ref_boundary(tree, flux_species+1)
       time = time + dt
       call forward_euler(tree, dt, time, 1, 1, 1, .true., 2)
       call combine_substeps(tree, &
            species_itree(n_gas_species+1:n_species), &
            [0, 1], [0.5_dp, 0.5_dp], 0)
       call af_restrict_ref_boundary(tree, flux_species)
    case default
       error stop "Invalid time integrator"
    end select

    ! Determine next time step
    advance_max_dt = min(2 * advance_max_dt, dt_safety_factor * &
         minval(dt_matrix(1:dt_num_cond, :)))

    call field_compute(tree, mg, 0, time, .true.)
  end subroutine advance

  subroutine combine_substeps(tree, ivs, in_steps, coeffs, out_step)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: ivs(:)
    integer, intent(in)       :: in_steps(:)
    real(dp), intent(in)      :: coeffs(:)
    integer, intent(in)       :: out_step
    integer                   :: lvl, i, id, n, nc
    real(dp), allocatable     :: tmp(DTIMES(:), :)

    nc = tree%n_cell

    !$omp parallel private(lvl, i, id, tmp, n)
    allocate(tmp(DTIMES(0:nc+1), size(ivs)))

    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)

          tmp = 0.0_dp
          do n = 1, size(in_steps)
             tmp = tmp + coeffs(n) * &
                  tree%boxes(id)%cc(DTIMES(:), ivs+in_steps(n))
          end do
          tree%boxes(id)%cc(DTIMES(:), ivs+out_step) = tmp
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine combine_substeps

end module m_advance
