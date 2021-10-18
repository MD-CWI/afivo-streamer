#include "cpp_macros.h"
!> This module contains routines related to , which can interpolate
!> 'to' the grid and 'from' the grid (useful for e.g. particle simulations). The
!> interpolation for meshes is called prolongation, see m_aX_prolong.
!>
!> Note that the particle coordinates are transferred via subroutines. This has
!> two advantages: first, data does not need to be copied, saving memory.
!> Second, the particle id can be cached in the particle's data structured.
module m_af_particles
  use m_af_types

  implicit none
  private

  abstract interface
     !> To get a particle id
     subroutine subr_particle_id(ix, id)
       import
       integer, intent(in)  :: ix !< Particle index
       integer, intent(out) :: id !< Particle id
     end subroutine subr_particle_id

     !> To get a particle's coordinates and weight
     subroutine subr_particle_rw(ix, r, w)
       import
       integer, intent(in)   :: ix      !< Particle index
       real(dp), intent(out) :: r(NDIM) !< Particle coordinates
       real(dp), intent(out) :: w       !< Particle weight
     end subroutine subr_particle_rw
  end interface

  public :: af_particles_to_grid

contains

  !> Map a list of particles to a density. The order can be zero (map particle
  !> to the containing cell) or one (use bi/tri-linear interpolation). Note that
  !> ghost cells are automatically filled by this routine.
  subroutine af_particles_to_grid(tree, iv, n_particles, get_id, get_rw, &
       order, density, fill_gc, iv_tmp, offset_particles)
    use m_af_restrict, only: af_restrict_tree
    use m_af_ghostcell, only: af_gc_tree
    use m_af_utils, only: af_get_id_at, af_tree_clear_cc, af_tree_clear_ghostcells
    type(af_t), intent(inout)     :: tree
    integer, intent(in)           :: iv          !< Variable to store density
    integer, intent(in)           :: n_particles !< The number of particles
    procedure(subr_particle_id)   :: get_id      !< To get the particle id
    procedure(subr_particle_rw)   :: get_rw      !< To get the particle position and weight
    integer, intent(in)           :: order       !< Order of interpolation
    !> Divide by cell area/volume (default: true)
    logical, intent(in), optional :: density
    !> Fill ghost cells afterwards (default: true)
    logical, intent(in), optional :: fill_gc
    !> Use temporary variable to convert to density. This can be faster, and is
    !> slightly more accurate for cylindrical coordinate systems, due to the way
    !> ghost cells are exchanged near refinement boundaries
    integer, intent(in), optional :: iv_tmp
    !> Start offset for indexing the particles (if zero, start at index 1)
    integer, intent(in), optional :: offset_particles

    integer              :: n, m, p_offset
    integer              :: current_thread, current_work
    integer              :: threads_left, work_left
    integer, allocatable :: ids(:)
    integer, allocatable :: npart_per_box(:)
    integer, allocatable :: box_threads(:)
    integer, allocatable :: threads(:)
    real(dp)             :: r(NDIM), weight
    logical              :: as_density
    logical              :: fill_gc_at_end
    logical              :: use_tmp_var

    allocate(ids(n_particles))
    allocate(npart_per_box(-1:tree%highest_id))
    allocate(box_threads(tree%highest_id))
    allocate(threads(n_particles))

    npart_per_box(:) = 0
    as_density = .true.
    if (present(density)) as_density = density
    fill_gc_at_end = .true.
    if (present(fill_gc)) fill_gc_at_end = fill_gc
    use_tmp_var = .false.
    if (present(iv_tmp)) then
       if (iv_tmp > 0) use_tmp_var = .true.
    end if
    p_offset = 0
    if (present(offset_particles)) p_offset = offset_particles

    if (use_tmp_var .and. .not. as_density) &
         error stop "Use iv_tmp only for density = .true."

    !$omp parallel do reduction(+:npart_per_box)
    do n = 1, n_particles
       call get_id(p_offset + n, ids(n))
       npart_per_box(ids(n)) = npart_per_box(ids(n)) + 1
    end do
    !$omp end parallel do

    if (sum(npart_per_box(-1:0)) > 0) then
       print *, "af_particles_to_grid: some are outside domain"
       m = 0
       do n = 1, n_particles
          if (ids(n) <= af_no_box) then
             call get_rw(p_offset + n, r, weight)
             print *, n, r
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

    ! Set density to zero in ghost cells
    call af_tree_clear_ghostcells(tree, iv)

    ! Set density to zero in all cells
    if (use_tmp_var) call af_tree_clear_cc(tree, iv_tmp)

    select case (order)
    case (0)
       if (use_tmp_var) then
          call particles_to_grid_0(tree, iv_tmp, get_rw, ids, &
               threads, n_particles, .false., p_offset)
          call add_as_density(tree, iv_tmp, iv)
       else
          call particles_to_grid_0(tree, iv, get_rw, ids, &
               threads, n_particles, as_density, p_offset)
       end if
    case (1)
       if (use_tmp_var) then
          call particles_to_grid_1(tree, iv_tmp, get_rw, ids, &
               threads, n_particles, .false., p_offset)
          call tree_add_from_ghostcells(tree, iv_tmp)
          call add_as_density(tree, iv_tmp, iv)
       else
          call particles_to_grid_1(tree, iv, get_rw, ids, &
               threads, n_particles, as_density, p_offset)
          call tree_add_from_ghostcells(tree, iv)
       end if
    case default
       error stop "af_particles_to_grid: Invalid interpolation order"
    end select

    call af_restrict_tree(tree, [iv])

    if (fill_gc_at_end) then
       if (.not. tree%has_cc_method(iv)) then
          print *, "Variable with index ", iv, "has no cc_method"
          print *, "do this with call af_set_cc_methods(tree, iv, ...)"
          error stop "af_particles_to_grid: no ghost cell method defined"
       else
          call af_gc_tree(tree, [iv])
       end if
    end if
  end subroutine af_particles_to_grid

  subroutine particles_to_grid_0(tree, iv, get_rw, ids, &
       threads, n_particles, density, p_offset)
    use omp_lib
    type(af_t), intent(inout) :: tree
    integer, intent(in)        :: iv !< Variable to store particle density
    integer, intent(in)        :: n_particles
    procedure(subr_particle_rw) :: get_rw
    integer, intent(in)        :: ids(n_particles)
    integer, intent(in)        :: threads(n_particles)
    logical, intent(in)        :: density
    integer, intent(in)        :: p_offset !< Offset for particle indexing
    integer                    :: n, thread_id, ix(NDIM)
    real(dp)                   :: r(NDIM), weight, inv_volume

    !$omp parallel private(n, thread_id, ix, inv_volume, r, weight)
    thread_id = omp_get_thread_num()

    do n = 1, n_particles
       if (threads(n) /= thread_id) cycle
       ! Handle this particle
       call get_rw(p_offset + n, r, weight)

       ix = af_cc_ix(tree%boxes(ids(n)), r)

       ! Fix indices for points exactly on the boundaries of a box (which could
       ! get a ghost cell index)
       where (ix < 1) ix = 1
       where (ix > tree%n_cell) ix = tree%n_cell

       if (density) then
#if NDIM == 2
          if (tree%coord_t == af_cyl) then
             inv_volume = 1 / af_cyl_volume_cc(tree%boxes(ids(n)), ix(1))
          else
             ! Cartesian
             inv_volume = 1 / product(tree%boxes(ids(n))%dr)
          end if
#else
          inv_volume = 1 / product(tree%boxes(ids(n))%dr)
#endif

          tree%boxes(ids(n))%cc(DINDEX(ix), iv) = &
               tree%boxes(ids(n))%cc(DINDEX(ix), iv) + &
               weight * inv_volume
       else
          tree%boxes(ids(n))%cc(DINDEX(ix), iv) = &
               tree%boxes(ids(n))%cc(DINDEX(ix), iv) + &
               weight
       end if
    end do
    !$omp end parallel
  end subroutine particles_to_grid_0

  !> Add weights to the cell centers using linear interpolation @todo Support
  !> cylindrical coordinates
  subroutine particles_to_grid_1(tree, iv, get_rw, ids, &
       threads, n_particles, density, p_offset)
    use omp_lib
    type(af_t), intent(inout) :: tree
    integer, intent(in)        :: iv !< Variable to store particle density
    integer, intent(in)        :: n_particles
    procedure(subr_particle_rw) :: get_rw
    integer, intent(in)        :: ids(n_particles)
    integer, intent(in)        :: threads(n_particles)
    logical, intent(in)        :: density !< Add particle as a density
    integer, intent(in)        :: p_offset !< Offset for particle indexing
    real(dp)                   :: tmp(NDIM), inv_dr(NDIM)
    real(dp)                   :: wu(NDIM), wl(NDIM), w(DTIMES(2))
    real(dp)                   :: inv_volume, r(NDIM), weight
    integer                    :: id, ix(NDIM), n, thread_id

    if (tree%coord_t == af_cyl .and. density) &
         error stop "For cylindrical coordinates, use iv_tmp"

    !$omp parallel private(n, inv_dr, tmp, thread_id, ix, id, wu, wl, &
    !$omp& w, inv_volume, r, weight)
    thread_id = omp_get_thread_num()

    do n = 1, n_particles
       if (threads(n) /= thread_id) cycle

       call get_rw(p_offset + n, r, weight)

       id     = ids(n)
       inv_dr = 1.0_dp/tree%boxes(id)%dr
       tmp    = (r - tree%boxes(id)%r_min) * inv_dr + 0.5_dp
       ix     = floor(tmp)
       wu     = tmp - ix
       wl     = 1 - wu

#if NDIM == 1
       w(1) = wl(1)
       w(2) = wu(1)
#elif NDIM == 2
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
          inv_volume = 1 / product(tree%boxes(ids(n))%dr)
#if NDIM == 1
          tree%boxes(id)%cc(ix(1):ix(1)+1, iv) = &
               tree%boxes(id)%cc(ix(1):ix(1)+1, iv) + &
               w * inv_volume * weight
#elif NDIM == 2
          tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, iv) = &
               tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, iv) + &
               w * inv_volume * weight
#elif NDIM == 3
          tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, ix(3):ix(3)+1, iv) = &
               tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, ix(3):ix(3)+1, iv) + &
               w * inv_volume * weight
#endif
       else
#if NDIM == 1
          tree%boxes(id)%cc(ix(1):ix(1)+1, iv) = &
               tree%boxes(id)%cc(ix(1):ix(1)+1, iv) + w * weight
#elif NDIM == 2
          tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, iv) = &
               tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, iv) + &
               w * weight
#elif NDIM == 3
          tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, ix(3):ix(3)+1, iv) = &
               tree%boxes(id)%cc(ix(1):ix(1)+1, ix(2):ix(2)+1, ix(3):ix(3)+1, iv) + &
               w * weight
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

          boxes(id)%cc(DSLICE(i0, i1), iv) = &
               boxes(id)%cc(DSLICE(i0, i1), iv) + &
               boxes(id)%cc(DSLICE(n0, n1), iv)
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

          boxes(id)%cc(DSLICE(i0, i1), iv) = &
               boxes(id)%cc(DSLICE(i0, i1), iv) + &
               boxes(nb_id)%cc(DSLICE(n0, n1), iv)
       end if
    end do; CLOSE_DO
  end subroutine add_from_ghostcells

  !> Convert particle weights to densities and add to another variable
  subroutine add_as_density(tree, iv_from, iv_to)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: iv_from
    integer, intent(in)       :: iv_to
    integer                   :: lvl, i, id

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call add_as_density_box(tree%boxes(id), iv_from, iv_to)
       end do
       !$omp end do nowait
    end do
    !$omp end parallel
  end subroutine add_as_density

  !> Convert particle weights to densities and add to another variable
  subroutine add_as_density_box(box, iv_from, iv_to)
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: iv_from
    integer, intent(in)        :: iv_to
    real(dp)                   :: inv_volume
#if NDIM == 2
    integer                    :: i
    real(dp), parameter        :: twopi = 2 * acos(-1.0_dp)
    real(dp)                   :: radius, inv_cyl
#endif

    inv_volume = 1.0_dp / product(box%dr)

#if NDIM == 2
    if (box%coord_t == af_cyl) then
       do i = 0, box%n_cell+1
          ! abs() accounts for ghost cell on other side of r = 0
          radius = abs(af_cyl_radius_cc(box, i))
          inv_cyl = inv_volume / (twopi * radius)
          box%cc(i, :, iv_to) = box%cc(i, :, iv_to) + &
               box%cc(i, :, iv_from) * inv_cyl
       end do
    else
       box%cc(DTIMES(:), iv_to) = box%cc(DTIMES(:), iv_to) + &
         box%cc(DTIMES(:), iv_from) * inv_volume
    end if
#else
    box%cc(DTIMES(:), iv_to) = box%cc(DTIMES(:), iv_to) + &
         box%cc(DTIMES(:), iv_from) * inv_volume
#endif
  end subroutine add_as_density_box

end module m_af_particles
