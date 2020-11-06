#include "cpp_macros.h"
!> This module contains routines related to the filling of ghost cells. Note that
!> corner ghost cells are not used in Afivo.
module m_af_ghostcell
  use m_af_types

  implicit none
  private

  public :: af_gc_tree
  public :: af_gc_ids
  public :: af_gc_box
  public :: af_gc_get_boundary_coords
  public :: af_bc_dirichlet_zero
  public :: af_bc_neumann_zero
  public :: af_bc_set_continuous
  public :: af_gc_interp
  public :: af_gc_prolong_copy
  public :: af_gc_interp_lim
  public :: af_gc2_box

contains

  !> Fill ghost cells for variables ivs on the sides of all boxes, using
  !> subr_rb on refinement boundaries and subr_bc on physical boundaries
  subroutine af_gc_tree(tree, ivs, corners, leaves_only)
    type(af_t), intent(inout)       :: tree        !< Tree to fill ghost cells on
    integer, intent(in)              :: ivs(:)      !< Variables for which ghost cells are set
    logical, intent(in), optional    :: corners     !< Fill corner ghost cells (default: yes)
    logical, intent(in), optional    :: leaves_only !< Fill only leaves' ghost cells (default: false)
    integer                          :: lvl
    logical                          :: all_ids

    if (.not. tree%ready) error stop "af_gc_tree: tree not ready"
    all_ids = .true.
    if (present(leaves_only)) all_ids = .not. leaves_only

    if (all_ids) then
       do lvl = 1, tree%highest_lvl
          call af_gc_ids(tree, tree%lvls(lvl)%ids, ivs, corners)
       end do
    else
       do lvl = 1, tree%highest_lvl
          call af_gc_ids(tree, tree%lvls(lvl)%leaves, ivs, corners)
       end do
    end if
  end subroutine af_gc_tree

  !> Fill ghost cells for variables ivs on the sides of all boxes, using subr_rb
  !> on refinement boundaries and subr_bc on physical boundaries. This routine
  !> assumes that ghost cells on other ids have been set already.
  subroutine af_gc_ids(tree, ids, ivs, corners)
    type(af_t), intent(inout)       :: tree    !< Tree to fill ghost cells on
    integer, intent(in)              :: ids(:)  !< Ids of boxes for which we set ghost cells
    integer, intent(in)              :: ivs(:)  !< Variables for which ghost cells are set
    logical, intent(in), optional    :: corners !< Fill corner ghost cells (default: yes)
    integer                          :: i

    !$omp parallel do
    do i = 1, size(ids)
       call af_gc_box(tree, ids(i), ivs, corners)
    end do
    !$omp end parallel do
  end subroutine af_gc_ids

  !> Fill ghost cells for variables ivs
  subroutine af_gc_box(tree, id, ivs, corners)
    type(af_t), intent(inout)     :: tree    !< Tree to fill ghost cells on
    integer, intent(in)           :: id      !< Id of box for which we set ghost cells
    integer, intent(in)           :: ivs(:)  !< Variables for which ghost cells are set
    logical, intent(in), optional :: corners !< Fill corner ghost cells (default: yes)
    logical                       :: do_corners
    integer                       :: i, iv
    integer                       :: nb, nb_id, bc_type
    integer                       :: lo(NDIM), hi(NDIM), dnb(NDIM)
    real(dp)                      :: coords(NDIM, tree%n_cell**(NDIM-1))
    real(dp)                      :: bc_val(tree%n_cell**(NDIM-1))

    do_corners = .true.
    if (present(corners)) do_corners = corners

    do i = 1, size(ivs)
       iv = ivs(i)
       if (.not. tree%has_cc_method(iv)) then
          print *, "For variable ", trim(tree%cc_names(iv))
          error stop "af_gc_box: no methods stored"
       end if
    end do

    do nb = 1, af_num_neighbors
       nb_id = tree%boxes(id)%neighbors(nb)

       if (nb_id > af_no_box) then
          ! There is a neighbor
          call af_get_index_bc_outside(nb, tree%boxes(id)%n_cell, 1, lo, hi)
          dnb = af_neighb_offset([nb])
          call copy_from_nb(tree%boxes(id), tree%boxes(nb_id), dnb, lo, hi, ivs)
       else if (nb_id == af_no_box) then
          ! Refinement boundary
          do i = 1, size(ivs)
             iv = ivs(i)
             call tree%cc_methods(iv)%rb(tree%boxes, id, nb, iv)
          end do
       else
          ! Physical boundary
          call af_gc_get_boundary_coords(tree%boxes(id), nb, coords)

          do i = 1, size(ivs)
             iv = ivs(i)
             if (associated(tree%cc_methods(iv)%bc_custom)) then
                call tree%cc_methods(iv)%bc_custom(tree%boxes(id), &
                     nb, iv, 1)
             else
                call tree%cc_methods(iv)%bc(tree%boxes(id), nb, iv, &
                     coords, bc_val, bc_type)
                call bc_to_gc(tree%boxes(id), nb, iv, bc_val, bc_type)
             end if
          end do
       end if
    end do

    if (do_corners) call af_gc_box_corner(tree%boxes, id, ivs)
  end subroutine af_gc_box

  !> Get coordinates at the faces of a box boundary
  subroutine af_gc_get_boundary_coords(box, nb, coords)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    real(dp), intent(out)   :: coords(NDIM, box%n_cell**(NDIM-1))
    integer                 :: i, nb_dim, bc_dim(NDIM-1)
    integer, parameter      :: all_dims(NDIM) = [(i, i = 1, NDIM)]
    real(dp)                :: rmin(NDIM)
#if NDIM == 3
    integer                 :: j, ix
#endif

    nb_dim       = af_neighb_dim(nb)
    bc_dim       = pack(all_dims, all_dims /= nb_dim)
    rmin(bc_dim) = box%r_min(bc_dim) + 0.5_dp * box%dr(bc_dim)

    if (af_neighb_low(nb)) then
       rmin(nb_dim) = box%r_min(nb_dim)
    else
       rmin(nb_dim) = box%r_min(nb_dim) + box%n_cell * box%dr(nb_dim)
    end if

#if NDIM == 1
    coords(nb_dim, 1) = rmin(nb_dim)
#elif NDIM == 2
    do i = 1, box%n_cell
       coords(bc_dim, i) = rmin(bc_dim) + (i-1) * box%dr(bc_dim)
       coords(nb_dim, i) = rmin(nb_dim)
    end do
#elif NDIM == 3
    ix = 0
    do j = 1, box%n_cell
       do i = 1, box%n_cell
          ix = ix + 1
          coords(bc_dim, ix) = rmin(bc_dim) + [i-1, j-1] * box%dr(bc_dim)
          coords(nb_dim, ix) = rmin(nb_dim)
       end do
    end do
#endif
  end subroutine af_gc_get_boundary_coords

  !> Fill corner ghost cells for variable iv on corners/edges of a box. If there
  !> is no box to copy the data from, use linear extrapolation. This routine
  !> assumes ghost cells on the sides of the box are available.
  subroutine af_gc_box_corner(boxes, id, ivs)
    type(box_t), intent(inout)  :: boxes(:)              !< List of all the boxes
    integer, intent(in)          :: id                    !< Id of box for which we set ghost cells
    integer, intent(in)          :: ivs(:)                !< Variable for which ghost cells are set
    integer                      :: n, nb_id, dnb(NDIM), lo(NDIM)
#if NDIM == 3
    integer                      :: hi(NDIM), dim
#endif

#if NDIM == 3
    ! Have to do edges before corners (since extrapolation can depend on edge values)
    do n = 1, af_num_edges
       dim = af_edge_dim(n)

       ! Check whether there is a neighbor, and find its index
       nb_id = boxes(id)%neighbor_mat(af_edge_dir(1, n), &
            af_edge_dir(2, n), af_edge_dir(3, n))

       lo = af_edge_min_ix(:, n) * (boxes(id)%n_cell + 1)
       lo(dim) = 1

       if (nb_id > af_no_box) then
          hi      = lo
          hi(dim) = boxes(id)%n_cell
          dnb   = af_neighb_offset(af_nb_adj_edge(:, n))
          call copy_from_nb(boxes(id), boxes(nb_id), dnb, lo, hi, ivs)
       else
          call af_edge_gc_extrap(boxes(id), lo, dim, ivs)
       end if
    end do
#endif

    do n = 1, af_num_children
       ! Check whether there is a neighbor, and find its index
       dnb   = 2 * af_child_dix(:, n) - 1

       nb_id = boxes(id)%neighbor_mat(DINDEX(dnb))
       lo    = af_child_dix(:, n) * (boxes(id)%n_cell + 1)

       if (nb_id > af_no_box) then
          call copy_from_nb(boxes(id), boxes(nb_id), dnb, lo, lo, ivs)
       else
          call af_corner_gc_extrap(boxes(id), lo, ivs)
       end if
    end do
  end subroutine af_gc_box_corner

  !> Convert a boundary condition to ghost cell data
  subroutine bc_to_gc(box, nb, iv, bc_val, bc_type)
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: iv      !< Variable to fill
    integer, intent(in)        :: nb      !< Neighbor direction
    real(dp), intent(in)       :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(in)        :: bc_type !< Type of b.c.
    real(dp)                   :: c0, c1, c2
    integer                    :: nc

    nc = box%n_cell

    ! If we call the interior point x1, x2 and the ghost point x0, then a
    ! Dirichlet boundary value b can be imposed as:
    ! x0 = -x1 + 2*b
    ! A Neumann b.c. can be imposed as:
    ! x0 = x1 +/- dx * b
    ! A continuous boundary (same slope) as:
    ! x0 = 2 * x1 - x2
    ! Below, we set coefficients to handle these cases
    select case (bc_type)
    case (af_bc_dirichlet)
       c0 = 2
       c1 = -1
       c2 = 0
    case (af_bc_neumann)
       c0 = box%dr(af_neighb_dim(nb)) * af_neighb_high_pm(nb) ! This gives a + or - sign
       c1 = 1
       c2 = 0
    case (af_bc_continuous)
       c0 = 0
       c1 = 2
       c2 = -1
    case (af_bc_dirichlet_copy)
       c0 = 1
       c1 = 0
       c2 = 0
    case default
       stop "fill_bc: unknown boundary condition"
    end select

    select case (nb)
#if NDIM == 1
    case (af_neighb_lowx)
       box%cc(0, iv) = &
            c0 * bc_val(1) + &
            c1 * box%cc(1, iv) + &
            c2 * box%cc(2, iv)
    case (af_neighb_highx)
       box%cc(nc+1, iv) = &
            c0 * bc_val(1) + &
            c1 * box%cc(nc, iv) + &
            c2 * box%cc(nc-1, iv)
#elif NDIM == 2
    case (af_neighb_lowx)
       box%cc(0, 1:nc, iv) = &
            c0 * bc_val + &
            c1 * box%cc(1, 1:nc, iv) + &
            c2 * box%cc(2, 1:nc, iv)
    case (af_neighb_highx)
       box%cc(nc+1, 1:nc, iv) = &
            c0 * bc_val + &
            c1 * box%cc(nc, 1:nc, iv) + &
            c2 * box%cc(nc-1, 1:nc, iv)
    case (af_neighb_lowy)
       box%cc(1:nc, 0, iv) = &
            c0 * bc_val + &
            c1 * box%cc(1:nc, 1, iv) + &
            c2 * box%cc(1:nc, 2, iv)
    case (af_neighb_highy)
       box%cc(1:nc, nc+1, iv) = &
            c0 * bc_val + &
            c1 * box%cc(1:nc, nc, iv) + &
            c2 * box%cc(1:nc, nc-1, iv)
#elif NDIM == 3
    case (af_neighb_lowx)
       box%cc(0, 1:nc, 1:nc, iv) = &
            c0 * reshape(bc_val, [nc, nc]) + &
            c1 * box%cc(1, 1:nc, 1:nc, iv) + &
            c2 * box%cc(2, 1:nc, 1:nc, iv)
    case (af_neighb_highx)
       box%cc(nc+1, 1:nc, 1:nc, iv) = &
            c0 * reshape(bc_val, [nc, nc]) + &
            c1 * box%cc(nc, 1:nc, 1:nc, iv) + &
            c2 * box%cc(nc-1, 1:nc, 1:nc, iv)
    case (af_neighb_lowy)
       box%cc(1:nc, 0, 1:nc, iv) = &
            c0 * reshape(bc_val, [nc, nc]) + &
            c1 * box%cc(1:nc, 1, 1:nc, iv) + &
            c2 * box%cc(1:nc, 2, 1:nc, iv)
    case (af_neighb_highy)
       box%cc(1:nc, nc+1, 1:nc, iv) = &
            c0 * reshape(bc_val, [nc, nc]) + &
            c1 * box%cc(1:nc, nc, 1:nc, iv) + &
            c2 * box%cc(1:nc, nc-1, 1:nc, iv)
    case (af_neighb_lowz)
       box%cc(1:nc, 1:nc, 0, iv) = &
            c0 * reshape(bc_val, [nc, nc]) + &
            c1 * box%cc(1:nc, 1:nc, 1, iv) + &
            c2 * box%cc(1:nc, 1:nc, 2, iv)
    case (af_neighb_highz)
       box%cc(1:nc, 1:nc, nc+1, iv) = &
            c0 * reshape(bc_val, [nc, nc]) + &
            c1 * box%cc(1:nc, 1:nc, nc, iv) + &
            c2 * box%cc(1:nc, 1:nc, nc-1, iv)
#endif
    end select
  end subroutine bc_to_gc

  !> Convert a boundary condition to two layers of ghost cell data
  subroutine bc_to_gc2(nc, cc, nb, bc_val, bc_type, dr)
    integer, intent(in)     :: nc                   !< Number of cells
    real(dp), intent(inout) :: cc(DTIMES(-1:nc+2))  !< Cell-centered data
    integer, intent(in)     :: nb                   !< Neighbor direction
    real(dp), intent(in)    :: bc_val(nc**(NDIM-1)) !< Boundary condition
    integer, intent(in)     :: bc_type              !< Type of b.c.
    real(dp), intent(in)    :: dr(NDIM)             !< Grid spacing
    real(dp)                :: c0, c1, c2

    ! If we call the interior point x1, x2 and the ghost point x0, then a
    ! Dirichlet boundary value b can be imposed as:
    ! x0 = -x1 + 2*b = c0 * b + c1 * x1
    ! x-1 = -x2 + 2*b = c2 * b + c1 * x2
    ! A Neumann b.c. can be imposed as:
    ! x0 = x1 +/- dx * b = c0 * b + c1 * x1
    ! x-1 = x2 +/- 3 * dx * b = c2 * b + c1 * x2
    !
    ! The second ghost cell here is a copy of the first one, this might not
    ! always be ideal, but it ensures the af_bc_dirichlet_copy variant does not
    ! introduce negative values
    !
    ! Below, we set coefficients to handle these cases
    select case (bc_type)
    case (af_bc_dirichlet)
       c0 = 2
       c1 = -1
       c2 = c0
    case (af_bc_neumann)
       c0 = dr(af_neighb_dim(nb)) * af_neighb_high_pm(nb) ! This gives a + or - sign
       c1 = 1
       c2 = 3 * c0
    case (af_bc_dirichlet_copy)
       c0 = 1
       c1 = 0
       c2 = c0
    case default
       stop "fill_bc: unknown boundary condition"
    end select

    select case (nb)
#if NDIM == 1
    case (af_neighb_lowx)
       cc(0) = c0 * bc_val(1) + c1 * cc(1)
       cc(-1) = c2 * bc_val(1) + c1 * cc(2)
    case (af_neighb_highx)
       cc(nc+1) = c0 * bc_val(1) + c1 * cc(nc)
       cc(nc+2) = c2 * bc_val(1) + c1 * cc(nc-1)
#elif NDIM == 2
    case (af_neighb_lowx)
       cc(0, 1:nc) = c0 * bc_val + c1 * cc(1, 1:nc)
       cc(-1, 1:nc) = c2 * bc_val + c1 * cc(2, 1:nc)
    case (af_neighb_highx)
       cc(nc+1, 1:nc) = c0 * bc_val + c1 * cc(nc, 1:nc)
       cc(nc+2, 1:nc) = c2 * bc_val + c1 * cc(nc-1, 1:nc)
    case (af_neighb_lowy)
       cc(1:nc, 0) = c0 * bc_val + c1 * cc(1:nc, 1)
       cc(1:nc, -1) = c2 * bc_val + c1 * cc(1:nc, 2)
    case (af_neighb_highy)
       cc(1:nc, nc+1) = c0 * bc_val + c1 * cc(1:nc, nc)
       cc(1:nc, nc+2) = c2 * bc_val + c1 * cc(1:nc, nc-1)
#elif NDIM == 3
    case (af_neighb_lowx)
       cc(0, 1:nc, 1:nc) = c0 * reshape(bc_val, [nc, nc]) + &
            c1 * cc(1, 1:nc, 1:nc)
       cc(-1, 1:nc, 1:nc) = c2 * reshape(bc_val, [nc, nc]) + &
            c1 * cc(2, 1:nc, 1:nc)
    case (af_neighb_highx)
       cc(nc+1, 1:nc, 1:nc) = c0 * reshape(bc_val, [nc, nc]) + &
            c1 * cc(nc, 1:nc, 1:nc)
       cc(nc+2, 1:nc, 1:nc) = c2 * reshape(bc_val, [nc, nc]) + &
            c1 * cc(nc-1, 1:nc, 1:nc)
    case (af_neighb_lowy)
       cc(1:nc, 0, 1:nc) = c0 * reshape(bc_val, [nc, nc]) + &
            c1 * cc(1:nc, 1, 1:nc)
       cc(1:nc, -1, 1:nc) = c2 * reshape(bc_val, [nc, nc]) + &
            c1 * cc(1:nc, 2, 1:nc)
    case (af_neighb_highy)
       cc(1:nc, nc+1, 1:nc) = c0 * reshape(bc_val, [nc, nc]) + &
            c1 * cc(1:nc, nc, 1:nc)
       cc(1:nc, nc+2, 1:nc) = c0 * reshape(bc_val, [nc, nc]) + &
            c1 * cc(1:nc, nc-1, 1:nc)
    case (af_neighb_lowz)
       cc(1:nc, 1:nc, 0) = c0 * reshape(bc_val, [nc, nc]) + &
            c1 * cc(1:nc, 1:nc, 1)
       cc(1:nc, 1:nc, -1) = c2 * reshape(bc_val, [nc, nc]) + &
            c1 * cc(1:nc, 1:nc, 2)
    case (af_neighb_highz)
       cc(1:nc, 1:nc, nc+1) = c0 * reshape(bc_val, [nc, nc]) + &
            c1 * cc(1:nc, 1:nc, nc)
       cc(1:nc, 1:nc, nc+2) = c2 * reshape(bc_val, [nc, nc]) + &
            c1 * cc(1:nc, 1:nc, nc-1)
#endif
    end select
  end subroutine bc_to_gc2

  !> Partial prolongation to the ghost cells of box id from parent
  subroutine af_gc_prolong_copy(boxes, id, nb, iv)
    use m_af_prolong, only: af_prolong_copy
    type(box_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Id of child
    integer, intent(in)           :: iv       !< Variable to fill
    integer, intent(in)           :: nb       !< Neighbor to get data from
    integer                       :: p_id, lo(NDIM), hi(NDIM)

    p_id = boxes(id)%parent
    call af_get_index_bc_outside(nb, boxes(id)%n_cell, 1, lo, hi)
    call af_prolong_copy(boxes(p_id), boxes(id), iv, low=lo, high=hi)
  end subroutine af_gc_prolong_copy

  !> Interpolate between fine points and coarse neighbors to fill ghost cells
  !> near refinement boundaries
  subroutine af_gc_interp(boxes, id, nb, iv)
    type(box_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc, ix, ix_c, ix_f
    integer                     :: p_nb_id
    integer                     :: p_id, ix_offset(NDIM)
    real(dp), parameter         :: third=1/3.0_dp
#if NDIM > 1
    real(dp), parameter         :: sixth=1/6.0_dp
    integer                     :: i_c1, i_c2, j_c1, j_c2, i, j
#endif
#if NDIM > 2
    integer                     :: k_c1, k_c2, k
#endif

    nc        = boxes(id)%n_cell
    p_id      = boxes(id)%parent
    p_nb_id   = boxes(p_id)%neighbors(nb)
    ix_offset = af_get_child_offset(boxes(id), nb)

    if (af_neighb_low(nb)) then
       ix = 0
       ix_f = 1
       ix_c = nc
    else
       ix = nc+1
       ix_f = nc
       ix_c = 1
    end if

    select case (af_neighb_dim(nb))
#if NDIM == 1
    case (1)
          boxes(id)%cc(ix, iv) = &
               2 * third * boxes(p_nb_id)%cc(ix_c, iv) + &
               third * boxes(id)%cc(ix_f, iv)
#elif NDIM == 2
    case (1)
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          boxes(id)%cc(ix, j, iv) = &
               0.5_dp * boxes(p_nb_id)%cc(ix_c, j_c1, iv) + &
               sixth * boxes(p_nb_id)%cc(ix_c, j_c2, iv) + &
               third * boxes(id)%cc(ix_f, j, iv)
       end do
    case (2)
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
          boxes(id)%cc(i, ix, iv) = &
               0.5_dp * boxes(p_nb_id)%cc(i_c1, ix_c, iv) + &
               sixth * boxes(p_nb_id)%cc(i_c2, ix_c, iv) + &
               third * boxes(id)%cc(i, ix_f, iv)
       end do
#elif NDIM==3
    case (1)
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)          ! even: +1, odd: -1
          do j = 1, nc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
             boxes(id)%cc(ix, j, k, iv) = &
                  third * boxes(p_nb_id)%cc(ix_c, j_c1, k_c1, iv) + &
                  sixth * boxes(p_nb_id)%cc(ix_c, j_c2, k_c1, iv) + &
                  sixth * boxes(p_nb_id)%cc(ix_c, j_c1, k_c2, iv) + &
                  third * boxes(id)%cc(ix_f, j, k, iv)
          end do
       end do
    case (2)
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             boxes(id)%cc(i, ix, k, iv) = &
                  third * boxes(p_nb_id)%cc(i_c1, ix_c, k_c1, iv) + &
                  sixth * boxes(p_nb_id)%cc(i_c2, ix_c, k_c1, iv) + &
                  sixth * boxes(p_nb_id)%cc(i_c1, ix_c, k_c2, iv) + &
                  third * boxes(id)%cc(i, ix_f, k, iv)
          end do
       end do
    case (3)
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             boxes(id)%cc(i, j, ix, iv) = &
                  third * boxes(p_nb_id)%cc(i_c1, j_c1, ix_c, iv) + &
                  sixth * boxes(p_nb_id)%cc(i_c1, j_c2, ix_c, iv) + &
                  sixth * boxes(p_nb_id)%cc(i_c2, j_c1, ix_c, iv) + &
                  third * boxes(id)%cc(i, j, ix_f, iv)
          end do
       end do
#endif
    end select

  end subroutine af_gc_interp

  !> Interpolate between fine points and coarse neighbors to fill ghost cells
  !> near refinement boundaries. The ghost values are less than twice the coarse
  !> values.
  subroutine af_gc_interp_lim(boxes, id, nb, iv)
    type(box_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)        :: id       !< Id of box
    integer, intent(in)        :: nb       !< Ghost cell direction
    integer, intent(in)        :: iv       !< Ghost cell variable
    integer                    :: nc, ix, ix_c, ix_f
    integer                    :: p_id, ix_offset(NDIM), p_nb_id
    real(dp)                   :: c(NDIM)
    real(dp), parameter        :: third = 1/3.0_dp
#if NDIM > 1
    integer                    :: IJK, IJK_(c1), IJK_(c2)
    real(dp), parameter        :: sixth = 1/6.0_dp
#endif

    nc        = boxes(id)%n_cell
    p_id      = boxes(id)%parent
    p_nb_id   = boxes(p_id)%neighbors(nb)
    ix_offset = af_get_child_offset(boxes(id), nb)

    if (af_neighb_low(nb)) then
       ix = 0
       ix_f = 1
       ix_c = nc
    else
       ix = nc+1
       ix_f = nc
       ix_c = 1
    end if

    select case (af_neighb_dim(nb))
#if NDIM == 1
    case (1)
       c(1) = boxes(p_nb_id)%cc(ix_c, iv)
       boxes(id)%cc(ix, iv) = (2 * c(1) + boxes(id)%cc(ix_f, iv)) * third
       if (boxes(id)%cc(ix, iv) > 2 * c(1)) boxes(id)%cc(ix, iv) = 2 * c(1)
#elif NDIM == 2
    case (1)
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          c(1) = boxes(p_nb_id)%cc(ix_c, j_c1, iv)
          c(2) = boxes(p_nb_id)%cc(ix_c, j_c2, iv)
          boxes(id)%cc(ix, j, iv) = 0.5_dp * c(1) + sixth * c(2) + &
               third * boxes(id)%cc(ix_f, j, iv)
          if (boxes(id)%cc(ix, j, iv) > 2 * c(1)) boxes(id)%cc(ix, j, iv) = 2 * c(1)
       end do
    case (2)
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
          c(1) = boxes(p_nb_id)%cc(i_c1, ix_c, iv)
          c(2) = boxes(p_nb_id)%cc(i_c2, ix_c, iv)
          boxes(id)%cc(i, ix, iv) = 0.5_dp * c(1) + sixth * c(2) + &
               third * boxes(id)%cc(i, ix_f, iv)
          if (boxes(id)%cc(i, ix, iv) > 2 * c(1)) boxes(id)%cc(i, ix, iv) = 2 * c(1)
       end do
#elif NDIM==3
    case (1)
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)          ! even: +1, odd: -1
          do j = 1, nc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
             c(1) = boxes(p_nb_id)%cc(ix_c, j_c1, k_c1, iv)
             c(2) = boxes(p_nb_id)%cc(ix_c, j_c2, k_c1, iv)
             c(3) = boxes(p_nb_id)%cc(ix_c, j_c1, k_c2, iv)
             boxes(id)%cc(ix, j, k, iv) = third * c(1) + sixth * c(2) + &
                  sixth * c(3) + third * boxes(id)%cc(ix_f, j, k, iv)
             if (boxes(id)%cc(ix, j, k, iv) > 2 * c(1)) &
                  boxes(id)%cc(ix, j, k, iv) = 2 * c(1)
          end do
       end do
    case (2)
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             c(1) = boxes(p_nb_id)%cc(i_c1, ix_c, k_c1, iv)
             c(2) = boxes(p_nb_id)%cc(i_c2, ix_c, k_c1, iv)
             c(3) = boxes(p_nb_id)%cc(i_c1, ix_c, k_c2, iv)
             boxes(id)%cc(i, ix, k, iv) = third * c(1) + sixth * c(2) + &
                  sixth * c(3) + third * boxes(id)%cc(i, ix_f, k, iv)
             if (boxes(id)%cc(i, ix, k, iv) > 2 * c(1)) &
                  boxes(id)%cc(i, ix, k, iv) = 2 * c(1)
          end do
       end do
    case (3)
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             c(1) = boxes(p_nb_id)%cc(i_c1, j_c1, ix_c, iv)
             c(2) = boxes(p_nb_id)%cc(i_c1, j_c2, ix_c, iv)
             c(3) = boxes(p_nb_id)%cc(i_c2, j_c1, ix_c, iv)
             boxes(id)%cc(i, j, ix, iv) = third * c(1) + sixth * c(2) + &
                  sixth * c(3) + third * boxes(id)%cc(i, j, ix_f, iv)
             if (boxes(id)%cc(i, j, ix, iv) > 2 * c(1)) &
                  boxes(id)%cc(i, j, ix, iv) = 2 * c(1)
          end do
       end do
#endif
    end select

  end subroutine af_gc_interp_lim

  !> This fills ghost cells near physical boundaries using Neumann zero
  subroutine af_bc_neumann_zero(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    bc_type = af_bc_neumann
    bc_val  = 0.0_dp
  end subroutine af_bc_neumann_zero

  !> This fills ghost cells near physical boundaries using Dirichlet zero
  subroutine af_bc_dirichlet_zero(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    bc_type = af_bc_dirichlet
    bc_val  = 0.0_dp
  end subroutine af_bc_dirichlet_zero

  !> This fills ghost cells near physical boundaries using the same slope
  subroutine af_bc_set_continuous(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    bc_type = af_bc_continuous
    ! Set values to zero (to prevent problems with NaN)
    bc_val  = 0.0_dp
  end subroutine af_bc_set_continuous

  subroutine copy_from_nb(box, box_nb, dnb, lo, hi, ivs)
    type(box_t), intent(inout) :: box     !< Box on which to fill ghost cells
    type(box_t), intent(in)    :: box_nb  !< Neighbouring box
    integer, intent(in)          :: dnb(NDIM) !< Neighbor spatial index offset
    integer, intent(in)          :: lo(NDIM)  !< Ghost cell low index
    integer, intent(in)          :: hi(NDIM)  !< Ghost cell high index
    integer, intent(in)          :: ivs(:)      !< Ghost cell variable
    integer                      :: nlo(NDIM), nhi(NDIM)

    ! Get indices on neighbor
    nlo = lo - dnb * box%n_cell
    nhi = hi - dnb * box%n_cell

    box%cc(DSLICE(lo, hi), ivs) = &
         box_nb%cc(DSLICE(nlo, nhi), ivs)
  end subroutine copy_from_nb

  !> Get array of cell-centered variables with multiple ghost cells, excluding corners
  subroutine af_gc2_box(tree, id, ivs, cc)
    type(af_t), intent(inout) :: tree   !< Tree to fill ghost cells on
    integer, intent(in)       :: id     !< Id of box for which we set ghost cells
    integer, intent(in)       :: ivs(:) !< Variables for which ghost cells are set
    real(dp)                  :: cc(DTIMES(-1:tree%n_cell+2), size(ivs))
    integer                   :: i, iv, nb, nc, nb_id, bc_type
    integer                   :: lo(NDIM), hi(NDIM), dnb(NDIM)
    integer                   :: nlo(NDIM), nhi(NDIM)
    real(dp)                  :: coords(NDIM, tree%n_cell**(NDIM-1))
    real(dp)                  :: bc_val(tree%n_cell**(NDIM-1))

    nc = tree%n_cell

    do i = 1, size(ivs)
       iv = ivs(i)

       ! Copy interior region
       cc(DTIMES(0:nc+1), i) = tree%boxes(id)%cc(DTIMES(:), iv)

       if (.not. tree%has_cc_method(iv)) then
          print *, "For variable ", trim(tree%cc_names(iv))
          error stop "No methods stored"
       end if
    end do

    do nb = 1, af_num_neighbors
       nb_id = tree%boxes(id)%neighbors(nb)

       if (nb_id > af_no_box) then
          ! There is a neighbor
          call af_get_index_bc_outside(nb, tree%n_cell, 2, lo, hi)

          dnb = af_neighb_offset([nb])
          nlo = lo - dnb * tree%n_cell
          nhi = hi - dnb * tree%n_cell

          cc(DSLICE(lo, hi), :) = &
               tree%boxes(nb_id)%cc(DSLICE(nlo, nhi), ivs)

       else if (nb_id == af_no_box) then
          ! Refinement boundary
          do i = 1, size(ivs)
             iv = ivs(i)
             call gc2_prolong_rb(tree%boxes, id, nb, iv, tree%n_cell, &
                  cc(DTIMES(:), i))
          end do
       else
          ! Physical boundary
          call af_gc_get_boundary_coords(tree%boxes(id), nb, coords)
          do i = 1, size(ivs)
             iv = ivs(i)
             if (associated(tree%cc_methods(iv)%bc_custom)) then
                call tree%cc_methods(iv)%bc_custom(tree%boxes(id), &
                     nb, iv, 2, cc(DTIMES(:), i))
             else
                call tree%cc_methods(iv)%bc(tree%boxes(id), &
                     nb, iv, coords, bc_val, bc_type)
                call bc_to_gc2(nc, cc(DTIMES(:), i), nb, bc_val, &
                     bc_type, tree%boxes(id)%dr)
             end if
          end do
       end if
    end do

    do i = 1, size(ivs)
       iv = ivs(i)

       ! Copy back updated ghost cells
       tree%boxes(id)%cc(DTIMES(:), iv) = cc(DTIMES(0:nc+1), i)
    end do

  end subroutine af_gc2_box

  !> Conservative prolongation with a limited slope for ghost cells near
  !> refinement boundaries
  !>
  !> @todo make compatible with arbitrary number of ghost cells
  subroutine gc2_prolong_rb(boxes, id, nb, iv, nc, cc)
    type(box_t), intent(in) :: boxes(:)            !< List of all boxes
    integer, intent(in)     :: id                  !< Index of box to fill ghost cells for
    integer, intent(in)     :: nb                  !< Neighbor direction
    integer, intent(in)     :: iv                  !< Index of variable
    integer, intent(in)     :: nc                  !< Number of cells
    real(dp), intent(inout) :: cc(DTIMES(-1:nc+2)) !< Enlarged array
    integer                 :: IJK, IJK_(f), p_nb_id
    integer                 :: lo_c(NDIM), hi_c(NDIM), ix_offset(NDIM)
    integer                 :: lo(NDIM), hi(NDIM)
    real(dp)                :: f(0:NDIM)

    p_nb_id = boxes(boxes(id)%parent)%neighbors(nb)

    ! Get index in current box
    call af_get_index_bc_outside(nb, nc, 2, lo, hi)

    ! Convert to index on parent box
    ix_offset = af_get_child_offset(boxes(id))
    lo_c = ix_offset + (lo+1)/2
    hi_c = ix_offset + (hi+1)/2

    ! Convert to index on neighbor of parent box
    lo_c = lo_c - af_neighb_dix(:, nb) * nc
    hi_c = hi_c - af_neighb_dix(:, nb) * nc

    associate(cc_p => boxes(p_nb_id)%cc)
#if NDIM == 1
      do i = lo_c(1), hi_c(1)
         i_f = lo(1) + 2 * (i - lo_c(1))

         ! Compute slopes on parent
         f(0) = cc_p(i, iv)
         f(1) = 0.25_dp * limit_slope( &
              cc_p(i, iv) - cc_p(i-1, iv), &
              cc_p(i+1, iv) - cc_p(i, iv))

         ! Prolong to fine cells
         cc(i_f) = f(0) - f(1)
         cc(i_f+1) = f(0) + f(1)
      end do
#elif NDIM == 2
      do j = lo_c(2), hi_c(2)
         j_f = lo(2) + 2 * (j - lo_c(2))
         do i = lo_c(1), hi_c(1)
            i_f = lo(1) + 2 * (i - lo_c(1))

            ! Compute slopes on parent
            f(0) = cc_p(i, j, iv)
            f(1) = 0.25_dp * limit_slope( &
                 cc_p(i, j, iv) - cc_p(i-1, j, iv), &
                 cc_p(i+1, j, iv) - cc_p(i, j, iv))
            f(2) = 0.25_dp * limit_slope( &
                 cc_p(i, j, iv) - cc_p(i, j-1, iv), &
                 cc_p(i, j+1, iv) - cc_p(i, j, iv))

            ! Prolong to fine cells
            cc(i_f,   j_f) = f(0) - f(1) - f(2)
            cc(i_f,   j_f+1) = f(0) - f(1) + f(2)
            cc(i_f+1, j_f) = f(0) + f(1) - f(2)
            cc(i_f+1, j_f+1) = f(0) + f(1) + f(2)
         end do
      end do
#elif NDIM == 3
      do k = lo_c(3), hi_c(3)
         k_f = lo(3) + 2 * (k - lo_c(3))
         do j = lo_c(2), hi_c(2)
            j_f = lo(2) + 2 * (j - lo_c(2))
            do i = lo_c(1), hi_c(1)
               i_f = lo(1) + 2 * (i - lo_c(1))

               ! Compute slopes on parent
               f(0) = cc_p(i, j, k, iv)
               f(1) = 0.25_dp * limit_slope( &
                    cc_p(i, j, k, iv) - cc_p(i-1, j, k, iv), &
                    cc_p(i+1, j, k, iv) - cc_p(i, j, k, iv))
               f(2) = 0.25_dp * limit_slope( &
                    cc_p(i, j, k, iv) - cc_p(i, j-1, k, iv), &
                    cc_p(i, j+1, k, iv) - cc_p(i, j, k, iv))
               f(3) = 0.25_dp * limit_slope( &
                    cc_p(i, j, k, iv) - cc_p(i, j, k-1, iv), &
                    cc_p(i, j, k+1, iv) - cc_p(i, j, k, iv))

               ! Prolong to fine cells
               cc(i_f,   j_f,   k_f)   = f(0) - f(1) - f(2) - f(3)
               cc(i_f,   j_f,   k_f+1) = f(0) - f(1) - f(2) + f(3)
               cc(i_f,   j_f+1, k_f)   = f(0) - f(1) + f(2) - f(3)
               cc(i_f,   j_f+1, k_f+1) = f(0) - f(1) + f(2) + f(3)
               cc(i_f+1, j_f,   k_f)   = f(0) + f(1) - f(2) - f(3)
               cc(i_f+1, j_f,   k_f+1) = f(0) + f(1) - f(2) + f(3)
               cc(i_f+1, j_f+1, k_f)   = f(0) + f(1) + f(2) - f(3)
               cc(i_f+1, j_f+1, k_f+1) = f(0) + f(1) + f(2) + f(3)
            end do
         end do
      end do
#endif
    end associate

  contains

    ! Take minimum of two slopes if they have the same sign, else take zero
    elemental function limit_slope(ll, rr) result(slope)
      real(dp), intent(in) :: ll, rr
      real(dp)             :: slope

      if (ll * rr < 0) then
         slope = 0.0_dp
      else
         ! MC limiter
         slope = sign(minval(abs([2 * ll, 2 * rr, 0.5_dp * (ll + rr)])), ll)
      end if
    end function limit_slope

  end subroutine gc2_prolong_rb

  !> This fills corner ghost cells using linear extrapolation. The ghost cells
  !> on the sides already need to be filled.
  subroutine af_corner_gc_extrap(box, ix, ivs)
    type(box_t), intent(inout) :: box !< Box to fill ghost cells for
    integer, intent(in)          :: ix(NDIM) !< Cell-index of corner
    integer, intent(in)          :: ivs(:)     !< Variable to fill
    integer                      :: di(NDIM)

    di = 1 - 2 * iand(ix, 1)    ! 0 -> di = 1, nc+1 -> di = -1

#if NDIM == 2
    box%cc(ix(1), ix(2), ivs) = box%cc(ix(1)+di(1), ix(2), ivs) &
         + box%cc(ix(1), ix(2)+di(2), ivs) &
         - box%cc(ix(1)+di(1), ix(2)+di(2), ivs)
#elif NDIM == 3
    box%cc(ix(1), ix(2), ix(3), ivs) = &
         box%cc(ix(1), ix(2)+di(2), ix(3)+di(3), ivs) + &
         box%cc(ix(1)+di(1), ix(2), ix(3)+di(3), ivs) + &
         box%cc(ix(1)+di(1), ix(2)+di(2), ix(3), ivs) - 2 * &
         box%cc(ix(1)+di(1), ix(2)+di(2), ix(3)+di(3), ivs)
#endif
  end subroutine af_corner_gc_extrap

#if NDIM == 3
  !> This fills edge ghost cells using linear extrapolation. The ghost cells on
  !> the sides already need to be filled. This routine basically performs the
  !> same operation as af_corner_gc_extrap does in 2D.
  subroutine af_edge_gc_extrap(box, lo, dim, ivs)
    type(box_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)          :: lo(NDIM) !< Lowest index of edge ghost cells
    integer, intent(in)          :: dim !< Dimension parallel to edge
    integer, intent(in)          :: ivs(:) !< Variable to fill
    integer                      :: di(NDIM), ix(NDIM), ia(NDIM), ib(NDIM), ic(NDIM)
    integer                      :: n, o_dims(NDIM-1)

    ! Dimensions other than/perpendicular to dim
    o_dims = [1 + mod(dim, NDIM), 1 + mod(dim + 1, NDIM)]

    ! Index offsets
    di = 1 - 2 * iand(lo, 1)    ! 0 -> di = 1, nc+1 -> di = -1
    di(dim) = 0

    ! Neighbor index in direction one
    ia = lo
    ia(o_dims(1)) = ia(o_dims(1)) + di(o_dims(1))

    ! Neighbor index in direction two
    ib = lo
    ib(o_dims(2)) = ib(o_dims(2)) + di(o_dims(2))

    ! Diagional neighbor index
    ic = lo + di
    ix = lo

    do n = 1, box%n_cell
       ia(dim) = n
       ib(dim) = n
       ic(dim) = n
       ix(dim) = n

       box%cc(ix(1), ix(2), ix(3), ivs) = &
            box%cc(ia(1), ia(2), ia(3), ivs) + &
            box%cc(ib(1), ib(2), ib(3), ivs) - &
            box%cc(ic(1), ic(2), ic(3), ivs)
    end do

  end subroutine af_edge_gc_extrap
#endif

end module m_af_ghostcell
