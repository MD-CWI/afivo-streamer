#include "cpp_macros_$Dd.h"
!> This module contains routines related to , which can interpolate
!> 'to' the grid and 'from' the grid (useful for e.g. particle simulations). The
!> interpolation for meshes is called prolongation, see m_aX_prolong.
module m_a$D_interp
  use m_a$D_types

  implicit none
  private

  public :: a$D_particles_to_grid

contains

  subroutine a$D_particles_to_grid(tree, iv, coords, weights, n_particles, order)
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)        :: iv !< Variable to store particle density
    integer, intent(in)        :: n_particles
    real(dp), intent(in)       :: coords($D, n_particles)
    real(dp), intent(in)       :: weights(n_particles)
    integer, intent(in)        :: order

    integer :: n
    integer, allocatable :: ids(:), npart_per_box(:)

    allocate(ids(n_particles))
    allocate(npart_per_box(-1:tree%highest_id))

    npart_per_box(:) = 0

    !$omp do
    do n = 1, n_particles
       ids(n) = a$D_get_id_at(tree, coords(:, n))
    end do
    !$omp end do

    do n = 1, n_particles
       npart_per_box(ids(n)) = npart_per_box(ids(n)) + 1
    end do

    select case (order)
    case (0)
       call particles_to_grid_0(tree, iv, coords, weights, n_particles)
    case (1)
       call particles_to_grid_1(tree, iv, coords, weights, n_particles)
    case default
       error stop "a$D_particles_to_grid: Invalid interpolation order"
    end select

    call a$D_restrict_tree(tree, iv)
    call a$D_gc_tree(tree, iv)
  end subroutine a$D_particles_to_grid

  subroutine particles_to_grid_0(tree, iv, coords, weights, n_particles)
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)        :: iv !< Variable to store particle density
    integer, intent(in)        :: n_particles
    real(dp), intent(in)       :: coords($D, n_particles)
    real(dp), intent(in)       :: weights(n_particles)
    integer :: n, loc

    !$omp do private(loc)
    do n = 1, n_particles
       
    end do
    !$omp end do
  end subroutine particles_to_grid_0

  !> Add 'amount' to the cell centers around rr, using linear interpolation
  subroutine a$D_interp1_to_grid(tree, rr, iv, amount, to_density)
    use m_a$D_utils, only: a$D_get_id_at
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)        :: iv         !< Index of variable
    real(dp), intent(in)       :: rr($D)     !< Location
    real(dp), intent(in)       :: amount     !< How much to add
    logical, intent(in)        :: to_density !< If true, divide by cell volume
    real(dp)                   :: actual_amount, inv_dr, tmp($D)
    real(dp)                   :: wu($D), wl($D), w(DTIMES(2))
    integer                    :: id, ix($D), nc
    logical                    :: rb_low($D), rb_high($D)

    id = a$D_get_id_at(tree, rr)

    if (id == -1) then
       print *, "a$D_interp1_to_grid error, no box at ", rr
       stop
    end if

    nc     = tree%boxes(id)%n_cell
    inv_dr = 1.0_dp/tree%boxes(id)%dr
    tmp    = (rr - tree%boxes(id)%r_min) * inv_dr + 0.5_dp
    ix     = floor(tmp)
    wu     = tmp - ix
    wl     = 1 - wu

#if $D == 2
    w(:, 1) = [wl(1), wu(1)] * wl(2)
    w(:, 2) = [wl(1), wu(1)] * wu(2)
#elif $D == 3
#endif

    !> @todo Support cylindrical coordinates
    if (to_density) then
       actual_amount = amount / tree%boxes(id)%dr**$D
    else
       actual_amount = amount
    end if

    !     if (any(ix == 0) .or. any(ix == nc)) then
    !        ix = ceiling((rr - tree%boxes(id)%r_min) * inv_dr)

    !        ! Apply special method near refinement boundaries
    ! #if $D == 2
    !        tree%boxes(id)%cc(ix(1), ix(2), iv) = &
    !             tree%boxes(id)%cc(ix(1), ix(2), iv) + &
    !             actual_amount
    ! #elif $D == 3
    !        tree%boxes(id)%cc(ix(1), ix(2), ix(3), iv) = &
    !             tree%boxes(id)%cc(ix(1), ix(2), ix(3), iv) + &
    !             actual_amount
    ! #endif
    !     else
    ! Linear interpolation
#if $D == 2
    tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, iv) = &
         tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, iv) + &
         w * actual_amount
#elif $D == 3
    ! tree%boxes(id)%cc(ix(1), ix(2), ix(3), iv) = &
    !      tree%boxes(id)%cc(ix(1), ix(2), ix(3), iv) + &
    !      actual_amount
#endif
    ! end if
  end subroutine a$D_interp1_to_grid

  subroutine a$D_interp1_copy_gc(tree, iv)
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)        :: iv !< Index of variable
    integer                    :: lvl, i, id

    do lvl = 1, tree%highest_lvl
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call a$D_copy_gc(tree%boxes, id, iv)
       end do
    end do
  end subroutine a$D_interp1_copy_gc

  subroutine a$D_copy_gc(boxes, id, iv)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)          :: id !< Index of box
    integer, intent(in)          :: iv !< Index of variable
    integer                      :: IJK, ix0($D), ix1($D)
    integer :: ni0($D), ni1($D), nb_id, nc
    logical :: copy_own

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
          ix0 = 1
          ix1 = nc
          ni0 = 1
          ni1 = nc

          where ([IJK] == 1)
             ix0 = nc
             ni0 = nc+1
             ni1 = nc+1
          elsewhere ([IJK] == -1)
             ix1 = 1
             ni0 = 0
             ni1 = 0
          end where

#if $D == 2
          boxes(id)%cc(ix0(1):ix1(1), ix0(2):ix1(2), iv) = &
               boxes(id)%cc(ix0(1):ix1(1), ix0(2):ix1(2), iv) + &
               boxes(id)%cc(ni0(1):ni1(1), ni0(2):ni1(2), iv)
#elif $D == 3
          error stop
#endif
       else
          ix0 = 1
          ix1 = nc
          ni0 = 1
          ni1 = nc
          where ([IJK] == 1)
             ix0 = nc
             ni0 = 0
             ni1 = 0
          elsewhere ([IJK] == -1)
             ix1 = 1
             ni0 = nc+1
             ni1 = nc+1
          end where

#if $D == 2
          boxes(id)%cc(ix0(1):ix1(1), ix0(2):ix1(2), iv) = &
               boxes(id)%cc(ix0(1):ix1(1), ix0(2):ix1(2), iv) + &
               boxes(nb_id)%cc(ni0(1):ni1(1), ni0(2):ni1(2), iv)
#elif $D == 3
          error stop
#endif
       end if
    end do; CLOSE_DO
  end subroutine a$D_copy_gc

end module m_a$D_interp
