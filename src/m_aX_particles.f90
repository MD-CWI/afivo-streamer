#include "cpp_macros_$Dd.h"
!> This module contains routines related to , which can interpolate
!> 'to' the grid and 'from' the grid (useful for e.g. particle simulations). The
!> interpolation for meshes is called prolongation, see m_aX_prolong.
module m_a$D_particles
  use m_a$D_types

  implicit none
  private

  public :: a$D_particles_to_grid

contains

  !> Map a list of particles to a density. The order can be zero (map particle
  !> to the containing cell) or one (use bi/tri-linear interpolation). Note that
  !> ghost cells are automatically filled by this routine.
  subroutine a$D_particles_to_grid(tree, iv, coords, weights, n_particles, order)
    use m_a$D_restrict, only: a$D_restrict_tree
    use m_a$D_ghostcell, only: a$D_gc_tree
    use m_a$D_utils, only: a$D_get_id_at, a$D_tree_clear_ghostcells
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)        :: iv !< Variable to store particle density
    integer, intent(in)        :: n_particles !< The number of particles
    real(dp), intent(in)       :: coords($D, n_particles) !< The particle coordinates
    real(dp), intent(in)       :: weights(n_particles) !< Weights for the particles
    integer, intent(in)        :: order !< Which order of interpolation to use

    integer              :: n, m
    integer              :: current_thread, current_work
    integer              :: threads_left, work_left
    integer, allocatable :: ids(:)
    integer, allocatable :: npart_per_box(:)
    integer, allocatable :: box_threads(:)
    integer, allocatable :: threads(:)

    allocate(ids(n_particles))
    allocate(npart_per_box(-1:tree%highest_id))
    allocate(box_threads(tree%highest_id))
    allocate(threads(n_particles))

    npart_per_box(:) = 0

    !$omp parallel do reduction(+:npart_per_box)
    do n = 1, n_particles
       ids(n) = a$D_get_id_at(tree, coords(:, n))
       npart_per_box(ids(n)) = npart_per_box(ids(n)) + 1
    end do
    !$omp end parallel do

    if (sum(npart_per_box(-1:0)) > 0) then
       error stop "Particles outside computational domain"
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
    call a$D_tree_clear_ghostcells(tree, iv)

    select case (order)
    case (0)
       call particles_to_grid_0(tree, iv, coords, weights, ids, &
            threads, n_particles)
    case (1)
       call particles_to_grid_1(tree, iv, coords, weights, ids, &
            threads, n_particles)
       call tree_add_from_ghostcells(tree, iv)
    case default
       error stop "a$D_particles_to_grid: Invalid interpolation order"
    end select

    call a$D_restrict_tree(tree, iv)
    call a$D_gc_tree(tree, iv)
  end subroutine a$D_particles_to_grid

  subroutine particles_to_grid_0(tree, iv, coords, weights, ids, &
       threads, n_particles)
    use omp_lib
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)        :: iv !< Variable to store particle density
    integer, intent(in)        :: n_particles
    real(dp), intent(in)       :: coords($D, n_particles)
    real(dp), intent(in)       :: weights(n_particles)
    integer, intent(in)        :: ids(n_particles)
    integer, intent(in)        :: threads(n_particles)
    integer                    :: n, thread_id, ix($D)

    !$omp parallel private(n, thread_id, ix)
    thread_id = omp_get_thread_num()

    do n = 1, n_particles
       if (threads(n) /= thread_id) cycle
       ! Handle this particle
       ix = a$D_cc_ix(tree%boxes(ids(n)), coords(:, n))

       ! Fix indices for points exactly on the boundaries of a box (which could
       ! get a ghost cell index)
       where (ix < 1) ix = 1
       where (ix > tree%n_cell) ix = tree%n_cell

#if $D == 2
       tree%boxes(ids(n))%cc(ix(1), ix(2), iv) = &
            tree%boxes(ids(n))%cc(ix(1), ix(2), iv) + &
            weights(n) / tree%boxes(ids(n))%dr**$D
#elif $D == 3
       tree%boxes(ids(n))%cc(ix(1), ix(2), ix(3), iv) = &
            tree%boxes(ids(n))%cc(ix(1), ix(2), ix(3), iv) + &
            weights(n) / tree%boxes(ids(n))%dr**$D
#endif
    end do
    !$omp end parallel
  end subroutine particles_to_grid_0

  !> Add weights to the cell centers using linear interpolation @todo Support
  !> cylindrical coordinates
  subroutine particles_to_grid_1(tree, iv, coords, weights, ids, &
       threads, n_particles)
    use omp_lib
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)        :: iv !< Variable to store particle density
    integer, intent(in)        :: n_particles
    real(dp), intent(in)       :: coords($D, n_particles)
    real(dp), intent(in)       :: weights(n_particles)
    integer, intent(in)        :: ids(n_particles)
    integer, intent(in)        :: threads(n_particles)
    real(dp)                   :: tmp($D), inv_dr
    real(dp)                   :: wu($D), wl($D), w(DTIMES(2))
    integer                    :: id, ix($D), n, thread_id

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

#if $D == 2
       w(:, 1) = [wl(1), wu(1)] * wl(2)
       w(:, 2) = [wl(1), wu(1)] * wu(2)
#elif $D == 3
       w(:, 1, 1) = [wl(1), wu(1)] * wl(2) * wl(3)
       w(:, 2, 1) = [wl(1), wu(1)] * wu(2) * wl(3)
       w(:, 1, 2) = [wl(1), wu(1)] * wl(2) * wu(3)
       w(:, 2, 2) = [wl(1), wu(1)] * wu(2) * wu(3)
#endif

       ! Linear interpolation
#if $D == 2
       tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, iv) = &
            tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, iv) + &
            w * weights(n) / tree%boxes(id)%dr**$D
#elif $D == 3
       tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, ix(3):ix(3)+1, iv) = &
            tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, ix(3):ix(3)+1, iv) + &
            w * weights(n) / tree%boxes(id)%dr**$D
#endif
    end do
    !$omp end parallel

  end subroutine particles_to_grid_1

  subroutine tree_add_from_ghostcells(tree, iv)
    type(a$D_t), intent(inout) :: tree
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
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)          :: id !< Index of box
    integer, intent(in)          :: iv !< Index of variable
    integer                      :: IJK, i0($D), i1($D)
    integer                      :: n0($D), n1($D), nb_id, nc
    logical                      :: copy_own

    nc = boxes(id)%n_cell

    do KJI_DO(-1,1)
       if (all([IJK] == 0)) cycle

       nb_id = boxes(id)%neighbor_mat(IJK)

       copy_own = .false.

       if (nb_id <= af_no_box) then
          copy_own = .true.
       else if (nb_id > af_no_box) then
          if (a$D_has_children(boxes(nb_id))) copy_own = .true.
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

#if $D == 2
          boxes(id)%cc(i0(1):i1(1), i0(2):i1(2), iv) = &
               boxes(id)%cc(i0(1):i1(1), i0(2):i1(2), iv) + &
               boxes(id)%cc(n0(1):n1(1), n0(2):n1(2), iv)
#elif $D == 3
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

#if $D == 2
          boxes(id)%cc(i0(1):i1(1), i0(2):i1(2), iv) = &
               boxes(id)%cc(i0(1):i1(1), i0(2):i1(2), iv) + &
               boxes(nb_id)%cc(n0(1):n1(1), n0(2):n1(2), iv)
#elif $D == 3
          boxes(id)%cc(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), iv) = &
               boxes(id)%cc(i0(1):i1(1), i0(2):i1(2), i0(3):i1(3), iv) + &
               boxes(nb_id)%cc(n0(1):n1(1), n0(2):n1(2), n0(3):n1(3), iv)
#endif
       end if
    end do; CLOSE_DO
  end subroutine add_from_ghostcells

end module m_a$D_particles
