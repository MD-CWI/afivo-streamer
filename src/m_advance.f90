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
  public :: advance_initialize
  public :: advance
  public :: advance_num_states

contains

  subroutine advance_initialize()
    advance_max_dt = dt_min
  end subroutine advance_initialize

  subroutine advance(tree, mg, dt, time)
    use m_chemistry
    use m_fluid_lfa
    use m_field
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(in)    :: mg
    real(dp), intent(in)      :: dt
    real(dp), intent(inout)   :: time

    dt_matrix(1:dt_num_cond, :) = dt_max ! Maximum time step

    select case (time_integrator)
    case (forward_euler_t)
       call forward_euler(tree, dt, 0, 0, .true.)
       call restrict_flux_species(tree, 0)
       call field_compute(tree, mg, 0, time, .true.)
       time = time + dt
    case (rk2_t)
       call forward_euler(tree, 0.5_dp * dt, 0, 1, .false.)
       call restrict_flux_species(tree, 1)
       time = time + 0.5_dp * dt
       call field_compute(tree, mg, 1, time, .true.)
       call forward_euler(tree, dt, 1, 0, .true.)
       call restrict_flux_species(tree, 0)
       call field_compute(tree, mg, 0, time, .true.)
    case (heuns_method_t)
       call forward_euler(tree, dt, 0, 1, .false.)
       call restrict_flux_species(tree, 1)
       time = time + dt
       call field_compute(tree, mg, 1, time, .true.)
       call forward_euler(tree, dt, 1, 1, .true.)
       call combine_substeps(tree, species_ix(1:n_species), &
            [0, 1], [0.5_dp, 0.5_dp], 0)
       call restrict_flux_species(tree, 0)
       call field_compute(tree, mg, 0, time, .true.)
    case default
       error stop "Invalid time integrator"
    end select

    ! Determine next time step
    advance_max_dt = min(2 * advance_max_dt, dt_safety_factor * &
         minval(dt_matrix(1:dt_num_cond, :)))
  end subroutine advance

  !> Restrict species for which we compute fluxes near refinement boundaries,
  !> for the coarse grid ghost cells
  subroutine restrict_flux_species(tree, s_out)
    use m_streamer
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: s_out
    integer                   :: lvl, i, id, p_id

    !$omp parallel private(lvl, i, id, p_id)
    do lvl = tree%lowest_lvl, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          p_id = tree%boxes(id)%parent
          if (p_id > af_no_box .and. &
               any(tree%boxes(id)%neighbors == af_no_box)) then
             call af_restrict_box_vars(tree%boxes(id), tree%boxes(p_id), &
                  flux_species + s_out)
          end if
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine restrict_flux_species

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
