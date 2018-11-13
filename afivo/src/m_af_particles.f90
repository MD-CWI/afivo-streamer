#include "cpp_macros.h"
!> This module contains routines related to , which can interpolate
!> 'to' the grid and 'from' the grid (useful for e.g. particle simulations). The
!> interpolation for meshes is called prolongation, see m_aX_prolong.
module m_af_particles
  use m_af_types

  implicit none
  private

  public :: af_particles_to_grid

contains

  !> Map a list of particles to a density. The order can be zero (map particle
  !> to the containing cell) or one (use bi/tri-linear interpolation). Note that
  !> ghost cells are automatically filled by this routine.
  subroutine af_particles_to_grid(tree, iv, coords, weights, n_particles, &
       order, id_guess, density, fill_gc)
    use m_af_restrict, only: af_restrict_tree
    use m_af_ghostcell, only: af_gc_tree
    use m_af_utils, only: af_get_id_at, af_tree_clear_ghostcells
    type(af_t), intent(inout)       :: tree
    integer, intent(in)              :: iv                      !< Variable to store density
    integer, intent(in)              :: n_particles             !< The number of particles
    real(dp), intent(in)             :: coords(NDIM, n_particles) !< The particle coordinates
    real(dp), intent(in)             :: weights(n_particles)    !< Weights for the particles
    integer, intent(in)              :: order                   !< Order of interpolation
    !> Guess for box id containing particle, set to 0 where no guess is available
    integer, intent(inout), optional :: id_guess(n_particles)
    !> Divide by cell area/volume (default: true)
    logical, intent(in), optional    :: density
    !> Fill ghost cells afterwards (default: true)
    logical, intent(in), optional    :: fill_gc

    integer              :: n, m
    integer              :: current_thread, current_work
    integer              :: threads_left, work_left
    integer, allocatable :: ids(:)
    integer, allocatable :: npart_per_box(:)
    integer, allocatable :: box_threads(:)
    integer, allocatable :: threads(:)
    logical              :: as_density
    logical              :: fill_gc_at_end

    allocate(ids(n_particles))
    allocate(npart_per_box(-1:tree%highest_id))
    allocate(box_threads(tree%highest_id))
    allocate(threads(n_particles))

    npart_per_box(:) = 0
    as_density = .true.
    if (present(density)) as_density = density
    fill_gc_at_end = .true.
    if (present(fill_gc)) fill_gc_at_end = fill_gc

    if (present(id_guess)) then
       !$omp parallel do reduction(+:npart_per_box)
       do n = 1, n_particles
          ids(n)                = af_get_id_at(tree, coords(:, n), &
               guess=id_guess(n))
          id_guess(n)           = ids(n)
          npart_per_box(ids(n)) = npart_per_box(ids(n)) + 1
       end do
       !$omp end parallel do
    else
       !$omp parallel do reduction(+:npart_per_box)
       do n = 1, n_particles
          ids(n)                = af_get_id_at(tree, coords(:, n))
          npart_per_box(ids(n)) = npart_per_box(ids(n)) + 1
       end do
       !$omp end parallel do
    end if

    if (sum(npart_per_box(-1:0)) > 0) then
       print *, "af_particles_to_grid: some are outside domain"
       m = 0
       do n = 1, n_particles
          if (ids(n) <= af_no_box) then
             print *, n, coords(:, n)
             m = m + 1
          end if
          if (m > 10) then
             print *, "..."
             exit
          end if
       end do
       error stop "af_particles_to_grid: some are outside domain"
    end if

    threads_left   = af_get_max_threads()
    current_thread = 0
    current_work   = 0
    work_left      = n_particles

    do m = 1, tree%highest_id
       box_threads(m) = current_thread
       current_work  = current_work + npart_per_box(m)

       if (current_work > work_left/threads_left) then
          current_thread = current_thread + 1
          threads_left   = threads_left - 1
          work_left      = work_left - current_work
          current_work   = 0
       end if
    end do

    !$omp parallel do
    do n = 1, n_particles
       threads(n) = box_threads(ids(n))
    end do
    !$omp end parallel do

    ! Set density to zero
    call af_tree_clear_ghostcells(tree, iv)

    select case (order)
    case (0)
       call particles_to_grid_0(tree, iv, coords, weights, ids, &
            threads, n_particles, as_density)
    case (1)
       call particles_to_grid_1(tree, iv, coords, weights, ids, &
            threads, n_particles, as_density)
       call tree_add_from_ghostcells(tree, iv)
    case default
       error stop "af_particles_to_grid: Invalid interpolation order"
    end select

    call af_restrict_tree(tree, iv)

    if (fill_gc_at_end) then
       if (.not. tree%has_cc_method(iv)) then
          print *, "Variable with index ", iv, "has no cc_method"
          print *, "do this with call af_set_cc_methods(tree, iv, ...)"
          error stop "af_particles_to_grid: no ghost cell method defined"
       else
          call af_gc_tree(tree, iv)
       end if
    end if
  end subroutine af_particles_to_grid

  subroutine particles_to_grid_0(tree, iv, coords, weights, ids, &
       threads, n_particles, density)
    use omp_lib
    type(af_t), intent(inout) :: tree
    integer, intent(in)        :: iv !< Variable to store particle density
    integer, intent(in)        :: n_particles
    real(dp), intent(in)       :: coords(NDIM, n_particles)
    real(dp), intent(in)       :: weights(n_particles)
    integer, intent(in)        :: ids(n_particles)
    integer, intent(in)        :: threads(n_particles)
    logical, intent(in)        :: density
    integer                    :: n, thread_id, ix(NDIM)

    !$omp parallel private(n, thread_id, ix)
    thread_id = omp_get_thread_num()

    do n = 1, n_particles
       if (threads(n) /= thread_id) cycle
       ! Handle this particle
       ix = af_cc_ix(tree%boxes(ids(n)), coords(:, n))

       ! Fix indices for points exactly on the boundaries of a box (which could
       ! get a ghost cell index)
       where (ix < 1) ix = 1
       where (ix > tree%n_cell) ix = tree%n_cell

       if (density) then
#if NDIM == 2
          tree%boxes(ids(n))%cc(ix(1), ix(2), iv) = &
               tree%boxes(ids(n))%cc(ix(1), ix(2), iv) + &
               weights(n) / tree%boxes(ids(n))%dr**NDIM
#elif NDIM == 3
          tree%boxes(ids(n))%cc(ix(1), ix(2), ix(3), iv) = &
               tree%boxes(ids(n))%cc(ix(1), ix(2), ix(3), iv) + &
               weights(n) / tree%boxes(ids(n))%dr**NDIM
#endif
       else
#if NDIM == 2
          tree%boxes(ids(n))%cc(ix(1), ix(2), iv) = &
               tree%boxes(ids(n))%cc(ix(1), ix(2), iv) + &
               weights(n)
#elif NDIM == 3
          tree%boxes(ids(n))%cc(ix(1), ix(2), ix(3), iv) = &
               tree%boxes(ids(n))%cc(ix(1), ix(2), ix(3), iv) + &
               weights(n)
#endif
       end if
    end do
    !$omp end parallel
  end subroutine particles_to_grid_0

  !> Add weights to the cell centers using linear interpolation @todo Support
  !> cylindrical coordinates
  subroutine particles_to_grid_1(tree, iv, coords, weights, ids, &
       threads, n_particles, density)
    use omp_lib
    type(af_t), intent(inout) :: tree
    integer, intent(in)        :: iv !< Variable to store particle density
    integer, intent(in)        :: n_particles
    real(dp), intent(in)       :: coords(NDIM, n_particles)
    real(dp), intent(in)       :: weights(n_particles)
    integer, intent(in)        :: ids(n_particles)
    integer, intent(in)        :: threads(n_particles)
    logical, intent(in)        :: density
    real(dp)                   :: tmp(NDIM), inv_dr
    real(dp)                   :: wu(NDIM), wl(NDIM), w(DTIMES(2))
    integer                    :: id, ix(NDIM), n, thread_id

    !$omp parallel private(n, inv_dr, tmp, thread_id, ix, id, wu, wl, w)
    thread_id = omp_get_thread_num()

    do n = 1, n_particles
       if (threads(n) /= thread_id) cycle

       id     = ids(n)
       inv_dr = 1.0_dp/tree%boxes(id)%dr
       tmp    = (coords(:, n) - tree%boxes(id)%r_min) * inv_dr + 0.5_dp
       ix     = floor(tmp)
       wu     = tmp - ix
       wl     = 1 - wu

#if NDIM == 2
       w(:, 1) = [wl(1), wu(1)] * wl(2)
       w(:, 2) = [wl(1), wu(1)] * wu(2)
#elif NDIM == 3
       w(:, 1, 1) = [wl(1), wu(1)] * wl(2) * wl(3)
       w(:, 2, 1) = [wl(1), wu(1)] * wu(2) * wl(3)
       w(:, 1, 2) = [wl(1), wu(1)] * wl(2) * wu(3)
       w(:, 2, 2) = [wl(1), wu(1)] * wu(2) * wu(3)
#endif

       ! Linear interpolation
       if (density) then
#if NDIM == 2
          tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, iv) = &
               tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, iv) + &
               w * weights(n) / tree%boxes(id)%dr**NDIM
#elif NDIM == 3
          tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, ix(3):ix(3)+1, iv) = &
               tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, ix(3):ix(3)+1, iv) + &
               w * weights(n) / tree%boxes(id)%dr**NDIM
#endif
       else
#if NDIM == 2
          tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, iv) = &
               tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, iv) + &
               w * weights(n)
#elif NDIM == 3
          tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, ix(3):ix(3)+1, iv) = &
               tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, ix(3):ix(3)+1, iv) + &
               w * weights(n)
#endif
       end if
    end do
    !$omp end parallel

  end subroutine particles_to_grid_1

  subroutine tree_add_from_ghostcells(tree, iv)
    type(af_t), intent(inout) :: tree
    integer, intent(in)        :: iv !< Index of variable
    integer                    :: lvl, i, id

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call add_from_ghostcells(tree%boxes, id, iv)
       end do
       !$omp end do nowait
    end do
    !$omp end parallel
  end subroutine tree_add_from_ghostcells

  subroutine add_from_ghostcells(boxes, id, iv)
    type(box_t), intent(inout) :: boxes(:)
    integer, intent(in)          :: id !< Index of box
    integer, intent(in)          :: iv !< Index of variable
    integer                      :: IJK, i0(NDIM), i1(NDIM)
    integer                      :: n0(NDIM), n1(NDIM), nb_id, nc
    logical                      :: copy_own

    nc = boxes(id)%n_cell

    do KJI_DO(-1,1)
       if (all([IJK] == 0)) cycle

       nb_id = boxes(id)%neighbor_mat(IJK)

       copy_own = .false.

       if (nb_id <= af_no_box) then
          copy_own = .true.
       else if (nb_id > af_no_box) then
          if (af_has_children(boxes(nb_id))) copy_own = .true.
       end if

       if (copy_own) then
          ! Physical boundary
          i0 = 1
          i1 = nc
          n0 = 1
          n1 = nc

          where ([IJK] == 1)
             i0 = nc
             n0 = nc+1
             n1 = nc+1
          elsewhere ([IJK] == -1)
             i1 = 1
             n0 = 0
             n1 = 0
          end where

#if NDIM == 2
          boxes(id)%cc(i0(1):i1(1), i0(2):i1(2), iv) = &
               boxes(id)%cc(i0(1):i1(1), i0(2):i1(2), iv) + &
               boxes(id)%cc(n0(1):n1(1), n0(2):n1(2), iv)
#elif NDIM == 3
          boxes(id)%cc(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), iv) = &
               boxes(id)%cc(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), iv) + &
               boxes(id)%cc(n0(1):n1(1), n0(2):n1(2), n0(3):n1(3), iv)
#endif
       else
          i0 = 1
          i1 = nc
          n0 = 1
          n1 = nc
          where ([IJK] == 1)
             i0 = nc
             n0 = 0
             n1 = 0
          elsewhere ([IJK] == -1)
             i1 = 1
             n0 = nc+1
             n1 = nc+1
          end where

#if NDIM == 2
          boxes(id)%cc(i0(1):i1(1), i0(2):i1(2), iv) = &
               boxes(id)%cc(i0(1):i1(1), i0(2):i1(2), iv) + &
               boxes(nb_id)%cc(n0(1):n1(1), n0(2):n1(2), iv)
#elif NDIM == 3
          boxes(id)%cc(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), iv) = &
               boxes(id)%cc(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), iv) + &
               boxes(nb_id)%cc(n0(1):n1(1), n0(2):n1(2), n0(3):n1(3), iv)
#endif
       end if
    end do; CLOSE_DO
  end subroutine add_from_ghostcells

end module m_af_particles
