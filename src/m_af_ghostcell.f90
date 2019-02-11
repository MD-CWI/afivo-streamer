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
  public :: af_gc2_prolong_linear
  public :: af_bc2_neumann_zero
  public :: af_bc2_dirichlet_zero

  ! Define interfaces so ghost cell routines can be called for a single variable
  ! or for multiple variables
  interface af_gc_tree
     module procedure af_gc_tree_iv, af_gc_tree_ivs
  end interface af_gc_tree

  interface af_gc_ids
     module procedure af_gc_ids_iv, af_gc_ids_ivs, af_gc_ids_v1
  end interface af_gc_ids

  interface af_gc_box
     module procedure af_gc_box_iv, af_gc_box_ivs, af_gc_box_v1
  end interface

contains

  !> Fill ghost cells for variables ivs on the sides of all boxes, using
  !> subr_rb on refinement boundaries and subr_bc on physical boundaries
  subroutine af_gc_tree_ivs(tree, ivs, subr_rb, subr_bc, corners, leaves_only)
    type(af_t), intent(inout)       :: tree        !< Tree to fill ghost cells on
    integer, intent(in)              :: ivs(:)      !< Variables for which ghost cells are set
    procedure(af_subr_rb), optional :: subr_rb     !< Procedure called at refinement boundaries
    procedure(af_subr_bc), optional :: subr_bc     !< Procedure called at physical boundaries
    logical, intent(in), optional    :: corners     !< Fill corner ghost cells (default: yes)
    logical, intent(in), optional    :: leaves_only !< Fill only leaves' ghost cells (default: false)
    integer                          :: lvl
    logical                          :: all_ids

    if (.not. tree%ready) error stop "af_gc_tree: tree not ready"
    all_ids = .true.
    if (present(leaves_only)) all_ids = .not. leaves_only

    if (all_ids) then
       do lvl = 1, tree%highest_lvl
          call af_gc_ids(tree, tree%lvls(lvl)%ids, ivs, &
               subr_rb, subr_bc, corners)
       end do
    else
       do lvl = 1, tree%highest_lvl
          call af_gc_ids(tree, tree%lvls(lvl)%leaves, ivs, &
               subr_rb, subr_bc, corners)
       end do
    end if
  end subroutine af_gc_tree_ivs

  !> Fill ghost cells for variables ivs on the sides of all boxes, using
  !> subr_rb on refinement boundaries and subr_bc on physical boundaries
  subroutine af_gc_tree_iv(tree, iv, subr_rb, subr_bc, corners, leaves_only)
    type(af_t), intent(inout)       :: tree        !< Tree to fill ghost cells on
    integer, intent(in)              :: iv          !< Variable for which ghost cells are set
    procedure(af_subr_rb), optional :: subr_rb     !< Procedure called at refinement boundaries
    procedure(af_subr_bc), optional :: subr_bc     !< Procedure called at physical boundaries
    logical, intent(in), optional    :: corners     !< Fill corner ghost cells (default: yes)
    logical, intent(in), optional    :: leaves_only !< Fill only leaves' ghost cells (default: false)
    integer                          :: lvl
    logical                          :: all_ids

    if (.not. tree%ready) error stop "af_gc_tree: tree not ready"
    all_ids = .true.
    if (present(leaves_only)) all_ids = .not. leaves_only

    if (all_ids) then
       do lvl = 1, tree%highest_lvl
          call af_gc_ids(tree, tree%lvls(lvl)%ids, iv, &
               subr_rb, subr_bc, corners)
       end do
    else
       do lvl = 1, tree%highest_lvl
          call af_gc_ids(tree, tree%lvls(lvl)%leaves, iv, &
               subr_rb, subr_bc, corners)
       end do
    end if
  end subroutine af_gc_tree_iv

  !> Fill ghost cells for variables ivs on the sides of all boxes, using subr_rb
  !> on refinement boundaries and subr_bc on physical boundaries. This routine
  !> assumes that ghost cells on other ids have been set already.
  subroutine af_gc_ids_ivs(tree, ids, ivs, subr_rb, subr_bc, corners)
    type(af_t), intent(inout)       :: tree    !< Tree to fill ghost cells on
    integer, intent(in)              :: ids(:)  !< Ids of boxes for which we set ghost cells
    integer, intent(in)              :: ivs(:)  !< Variables for which ghost cells are set
    procedure(af_subr_rb), optional :: subr_rb !< Procedure called at refinement boundaries
    procedure(af_subr_bc), optional :: subr_bc !< Procedure called at physical boundaries
    logical, intent(in), optional    :: corners !< Fill corner ghost cells (default: yes)
    integer                          :: i

    !$omp parallel do
    do i = 1, size(ids)
       call af_gc_box(tree, ids(i), ivs, subr_rb, subr_bc, corners)
    end do
    !$omp end parallel do
  end subroutine af_gc_ids_ivs

  !> Fill ghost cells for variables iv on the sides of all boxes, using subr_rb
  !> on refinement boundaries and subr_bc on physical boundaries. This routine
  !> assumes that ghost cells on other ids have been set already.
  subroutine af_gc_ids_iv(tree, ids, iv, subr_rb, subr_bc, corners)
    type(af_t), intent(inout)       :: tree    !< Tree to fill ghost cells on
    integer, intent(in)              :: ids(:)  !< Ids of boxes for which we set ghost cells
    integer, intent(in)              :: iv      !< Variable for which ghost cells are set
    procedure(af_subr_rb), optional :: subr_rb !< Procedure called at refinement boundaries
    procedure(af_subr_bc), optional :: subr_bc !< Procedure called at physical boundaries
    logical, intent(in), optional    :: corners !< Fill corner ghost cells (default: yes)
    integer                          :: i

    !$omp parallel do
    do i = 1, size(ids)
       call af_gc_box(tree, ids(i), iv, subr_rb, subr_bc, corners)
    end do
    !$omp end parallel do
  end subroutine af_gc_ids_iv

  !> Fill ghost cells for variables iv on the sides of all boxes, using subr_rb
  !> on refinement boundaries and subr_bc on physical boundaries. This routine
  !> assumes that ghost cells on other ids have been set already.
  subroutine af_gc_ids_v1(boxes, ids, iv, subr_rb, subr_bc, corners)
    type(box_t), intent(inout)     :: boxes(:) !< List of all the boxes
    integer, intent(in)              :: ids(:)   !< Ids of boxes for which we set ghost cells
    integer, intent(in)              :: iv       !< Variable for which ghost cells are set
    procedure(af_subr_rb), optional :: subr_rb  !< Procedure called at refinement boundaries
    procedure(af_subr_bc), optional :: subr_bc  !< Procedure called at physical boundaries
    logical, intent(in), optional    :: corners  !< Fill corner ghost cells (default: yes)
    integer                          :: i

    !$omp parallel do
    do i = 1, size(ids)
       call af_gc_box(boxes, ids(i), iv, subr_rb, subr_bc, corners)
    end do
    !$omp end parallel do
  end subroutine af_gc_ids_v1

  !> Fill ghost cells for variables ivs
  subroutine af_gc_box_ivs(tree, id, ivs, subr_rb, subr_bc, corners)
    type(af_t), intent(inout)       :: tree    !< Tree to fill ghost cells on
    integer, intent(in)              :: id      !< Id of box for which we set ghost cells
    integer, intent(in)              :: ivs(:)  !< Variables for which ghost cells are set
    procedure(af_subr_rb), optional :: subr_rb !< Procedure called at refinement boundaries
    procedure(af_subr_bc), optional :: subr_bc !< Procedure called at physical boundaries
    logical, intent(in), optional    :: corners !< Fill corner ghost cells (default: yes)
    logical                          :: do_corners
    integer                          :: i, iv
    procedure(af_subr_rb), pointer  :: use_rb
    procedure(af_subr_bc), pointer  :: use_bc

    do_corners = .true.
    if (present(corners)) do_corners = corners

    do i = 1, size(ivs)
       iv = ivs(i)

       if (present(subr_rb)) then
          use_rb => subr_rb
       else if (tree%has_cc_method(iv)) then
          use_rb => tree%cc_methods(iv)%rb
       else
          print *, "For variable ", trim(tree%cc_names(iv))
          error stop "af_gc_box: no refinement boundary method stored"
       end if

       if (present(subr_bc)) then
          use_bc => subr_bc
       else if (tree%has_cc_method(iv)) then
          use_bc => tree%cc_methods(iv)%bc
       else
          print *, "For variable ", trim(tree%cc_names(iv))
          error stop "af_gc_box: no boundary condition stored"
       end if

       call af_gc_box_sides(tree%boxes, id, iv, use_rb, use_bc)
       if (do_corners) call af_gc_box_corner(tree%boxes, id, iv)
    end do
  end subroutine af_gc_box_ivs

  !> Fill ghost cells for variable iv
  subroutine af_gc_box_iv(tree, id, iv, subr_rb, subr_bc, corners)
    type(af_t), intent(inout)       :: tree    !< Tree to fill ghost cells on
    integer, intent(in)              :: id      !< Id of box for which we set ghost cells
    integer, intent(in)              :: iv      !< Variable for which ghost cells are set
    procedure(af_subr_rb), optional :: subr_rb !< Procedure called at refinement boundaries
    procedure(af_subr_bc), optional :: subr_bc !< Procedure called at physical boundaries
    logical, intent(in), optional    :: corners !< Fill corner ghost cells (default: yes)
    logical                          :: do_corners
    procedure(af_subr_rb), pointer  :: use_rb
    procedure(af_subr_bc), pointer  :: use_bc

    if (present(subr_rb)) then
       use_rb => subr_rb
    else if (tree%has_cc_method(iv)) then
       use_rb => tree%cc_methods(iv)%rb
    else
       print *, "For variable ", trim(tree%cc_names(iv))
       error stop "af_gc_box: no refinement boundary method stored"
    end if

    if (present(subr_bc)) then
       use_bc => subr_bc
    else if (tree%has_cc_method(iv)) then
       use_bc => tree%cc_methods(iv)%bc
    else
       print *, "For variable ", trim(tree%cc_names(iv))
       error stop "af_gc_box: no boundary condition stored"
    end if

    do_corners = .true.
    if (present(corners)) do_corners = corners

    call af_gc_box_sides(tree%boxes, id, iv, use_rb, use_bc)
    if (do_corners) call af_gc_box_corner(tree%boxes, id, iv)
  end subroutine af_gc_box_iv

  !> Fill ghost cells for variable iv
  subroutine af_gc_box_v1(boxes, id, iv, subr_rb, subr_bc, corners)
    type(box_t), intent(inout)    :: boxes(:) !< List of all the boxes
    integer, intent(in)             :: id       !< Id of box for which we set ghost cells
    integer, intent(in)             :: iv       !< Variable for which ghost cells are set
    procedure(af_subr_rb)          :: subr_rb  !< Procedure called at refinement boundaries
    procedure(af_subr_bc)          :: subr_bc  !< Procedure called at physical boundaries
    logical, intent(in), optional   :: corners  !< Fill corner ghost cells (default: yes)
    logical                         :: do_corners

    do_corners = .true.
    if (present(corners)) do_corners = corners

    call af_gc_box_sides(boxes, id, iv, subr_rb, subr_bc)
    if (do_corners) call af_gc_box_corner(boxes, id, iv)
  end subroutine af_gc_box_v1

  !> Fill ghost cells for variable iv on the sides of a box, using subr_rb on
  !> refinement boundaries and subr_bc on physical boundaries.
  subroutine af_gc_box_sides(boxes, id, iv, subr_rb, subr_bc)
    type(box_t), intent(inout) :: boxes(:) !< List of all the boxes
    integer, intent(in)        :: id       !< Id of box for which we set ghost cells
    integer, intent(in)        :: iv       !< Variable for which ghost cells are set
    procedure(af_subr_rb)      :: subr_rb  !< Procedure called at refinement boundaries
    procedure(af_subr_bc)      :: subr_bc  !< Procedure called at physical boundaries
    integer                    :: nb, nb_id, bc_type
    integer                    :: lo(NDIM), hi(NDIM), dnb(NDIM)
    real(dp)                   :: coords(NDIM, boxes(id)%n_cell**(NDIM-1))
    real(dp)                   :: bc_val(boxes(id)%n_cell**(NDIM-1))

    do nb = 1, af_num_neighbors
       nb_id = boxes(id)%neighbors(nb)

       if (nb_id > af_no_box) then
          ! There is a neighbor
          call af_get_index_bc_outside(nb, boxes(id)%n_cell, lo, hi)
          dnb = af_neighb_offset([nb])
          call copy_from_nb(boxes(id), boxes(nb_id), dnb, lo, hi, iv)
       else if (nb_id == af_no_box) then
          ! Refinement boundary
          call subr_rb(boxes, id, nb, iv)
       else
          ! Physical boundary
          call af_gc_get_boundary_coords(boxes(id), nb, coords)
          call subr_bc(boxes(id), nb, iv, coords, bc_val, bc_type)
          call bc_to_gc(boxes(id), nb, iv, bc_val, bc_type)
       end if
    end do

  end subroutine af_gc_box_sides

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

#if NDIM == 2
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
  subroutine af_gc_box_corner(boxes, id, iv)
    type(box_t), intent(inout)  :: boxes(:)              !< List of all the boxes
    integer, intent(in)          :: id                    !< Id of box for which we set ghost cells
    integer, intent(in)          :: iv                    !< Variable for which ghost cells are set
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
          call copy_from_nb(boxes(id), boxes(nb_id), dnb, lo, hi, iv)
       else
          call af_edge_gc_extrap(boxes(id), lo, dim, iv)
       end if
    end do
#endif

    do n = 1, af_num_children
       ! Check whether there is a neighbor, and find its index
       dnb   = 2 * af_child_dix(:, n) - 1
#if NDIM == 2
       nb_id = boxes(id)%neighbor_mat(dnb(1), dnb(2))
#elif NDIM == 3
       nb_id = boxes(id)%neighbor_mat(dnb(1), dnb(2), dnb(3))
#endif
       lo    = af_child_dix(:, n) * (boxes(id)%n_cell + 1)

       if (nb_id > af_no_box) then
          call copy_from_nb(boxes(id), boxes(nb_id), dnb, lo, lo, iv)
       else
          call af_corner_gc_extrap(boxes(id), lo, iv)
       end if
    end do
  end subroutine af_gc_box_corner

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
    case default
       stop "fill_bc: unknown boundary condition"
    end select

    select case (nb)
#if NDIM == 2
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

  !> Partial prolongation to the ghost cells of box id from parent
  subroutine af_gc_prolong_copy(boxes, id, nb, iv)
    use m_af_prolong, only: af_prolong_copy
    type(box_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Id of child
    integer, intent(in)           :: iv       !< Variable to fill
    integer, intent(in)           :: nb       !< Neighbor to get data from
    integer                       :: p_id, lo(NDIM), hi(NDIM)

    p_id = boxes(id)%parent
    call af_get_index_bc_outside(nb, boxes(id)%n_cell, lo, hi)
    call af_prolong_copy(boxes(p_id), boxes(id), iv, low=lo, high=hi)
  end subroutine af_gc_prolong_copy

  !> Interpolate between fine points and coarse neighbors to fill ghost cells
  !> near refinement boundaries
  subroutine af_gc_interp(boxes, id, nb, iv)
    type(box_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc, ix, ix_c, ix_f, i, j
    integer                     :: i_c1, i_c2, j_c1, j_c2, p_nb_id
    integer                     :: p_id, ix_offset(NDIM)
    real(dp), parameter         :: sixth=1/6.0_dp, third=1/3.0_dp
#if NDIM == 3
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
#if NDIM == 2
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
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc, ix, ix_c, ix_f, i, j
    integer                     :: i_c1, i_c2, j_c1, j_c2, p_nb_id
    integer                     :: p_id, ix_offset(NDIM)
    real(dp)                    :: c1, c2
    real(dp), parameter         :: sixth=1/6.0_dp, third=1/3.0_dp
#if NDIM == 3
    integer                     :: k_c1, k_c2, k
    real(dp)                    :: c3
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
#if NDIM == 2
    case (1)
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          c1 = boxes(p_nb_id)%cc(ix_c, j_c1, iv)
          c2 = boxes(p_nb_id)%cc(ix_c, j_c2, iv)
          boxes(id)%cc(ix, j, iv) = 0.5_dp * c1 + sixth * c2 + &
               third * boxes(id)%cc(ix_f, j, iv)
          if (boxes(id)%cc(ix, j, iv) > 2 * c1) boxes(id)%cc(ix, j, iv) = 2 * c1
       end do
    case (2)
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
          c1 = boxes(p_nb_id)%cc(i_c1, ix_c, iv)
          c2 = boxes(p_nb_id)%cc(i_c2, ix_c, iv)
          boxes(id)%cc(i, ix, iv) = 0.5_dp * c1 + sixth * c2 + &
               third * boxes(id)%cc(i, ix_f, iv)
          if (boxes(id)%cc(i, ix, iv) > 2 * c1) boxes(id)%cc(i, ix, iv) = 2 * c1
       end do
#elif NDIM==3
    case (1)
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)          ! even: +1, odd: -1
          do j = 1, nc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
             c1 = boxes(p_nb_id)%cc(ix_c, j_c1, k_c1, iv)
             c2 = boxes(p_nb_id)%cc(ix_c, j_c2, k_c1, iv)
             c3 = boxes(p_nb_id)%cc(ix_c, j_c1, k_c2, iv)
             boxes(id)%cc(ix, j, k, iv) = third * c1 + sixth * c2 + &
                  sixth * c3 + third * boxes(id)%cc(ix_f, j, k, iv)
             if (boxes(id)%cc(ix, j, k, iv) > 2 * c1) &
                  boxes(id)%cc(ix, j, k, iv) = 2 * c1
          end do
       end do
    case (2)
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             c1 = boxes(p_nb_id)%cc(i_c1, ix_c, k_c1, iv)
             c2 = boxes(p_nb_id)%cc(i_c2, ix_c, k_c1, iv)
             c3 = boxes(p_nb_id)%cc(i_c1, ix_c, k_c2, iv)
             boxes(id)%cc(i, ix, k, iv) = third * c1 + sixth * c2 + &
                  sixth * c3 + third * boxes(id)%cc(i, ix_f, k, iv)
             if (boxes(id)%cc(i, ix, k, iv) > 2 * c1) &
                  boxes(id)%cc(i, ix, k, iv) = 2 * c1
          end do
       end do
    case (3)
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1
             c1 = boxes(p_nb_id)%cc(i_c1, j_c1, ix_c, iv)
             c2 = boxes(p_nb_id)%cc(i_c1, j_c2, ix_c, iv)
             c3 = boxes(p_nb_id)%cc(i_c2, j_c1, ix_c, iv)
             boxes(id)%cc(i, j, ix, iv) = third * c1 + sixth * c2 + &
                  sixth * c3 + third * boxes(id)%cc(i, j, ix_f, iv)
             if (boxes(id)%cc(i, j, ix, iv) > 2 * c1) &
                  boxes(id)%cc(i, j, ix, iv) = 2 * c1
          end do
       end do
#endif
    end select

  end subroutine af_gc_interp_lim

  ! This fills ghost cells near physical boundaries using Neumann zero
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

  ! This fills ghost cells near physical boundaries using Neumann zero
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

  ! This fills ghost cells near physical boundaries using the same slope
  subroutine af_bc_set_continuous(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    bc_type = af_bc_dirichlet
    ! Set values to zero (to prevent problems with NaN)
    bc_val  = 0.0_dp
  end subroutine af_bc_set_continuous

  subroutine copy_from_nb(box, box_nb, dnb, lo, hi, iv)
    type(box_t), intent(inout) :: box     !< Box on which to fill ghost cells
    type(box_t), intent(in)    :: box_nb  !< Neighbouring box
    integer, intent(in)          :: dnb(NDIM) !< Neighbor spatial index offset
    integer, intent(in)          :: lo(NDIM)  !< Ghost cell low index
    integer, intent(in)          :: hi(NDIM)  !< Ghost cell high index
    integer, intent(in)          :: iv      !< Ghost cell variable
    integer                      :: nlo(NDIM), nhi(NDIM)

    ! Get indices on neighbor
    nlo = lo - dnb * box%n_cell
    nhi = hi - dnb * box%n_cell

#if NDIM == 2
    box%cc(lo(1):hi(1), lo(2):hi(2), iv) = &
         box_nb%cc(nlo(1):nhi(1), nlo(2):nhi(2), iv)
#elif NDIM == 3
    box%cc(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), iv) = &
         box_nb%cc(nlo(1):nhi(1), nlo(2):nhi(2), nlo(3):nhi(3), iv)
#endif
  end subroutine copy_from_nb

  !> Get a second layer of ghost cell data (the 'normal' routines give just one
  !> layer of ghost cells). Use subr_rb on refinement boundaries and subr_bc
  !> on physical boundaries.
  subroutine af_gc2_box(boxes, id, iv, subr_rb, subr_bc, cc, nc)
    type(box_t), intent(inout) :: boxes(:) !< List of all the boxes
    integer, intent(in)          :: id       !< Id of box for which we set ghost cells
    integer, intent(in)          :: iv       !< Variable for which ghost cells are set
    procedure(af_subr_egc)      :: subr_rb  !< Procedure called at refinement boundaries
    procedure(af_subr_egc)      :: subr_bc  !< Procedure called at physical boundaries
    integer, intent(in)          :: nc       !< box%n_cell
    !> The data with extra ghost cells
#if NDIM   == 2
    real(dp), intent(out)        :: cc(-1:nc+2, -1:nc+2)
    real(dp)                     :: gc(1:nc)
#elif NDIM == 3
    real(dp), intent(out)        :: cc(-1:nc+2, -1:nc+2, -1:nc+2)
    real(dp)                     :: gc(1:nc, 1:nc)
#endif
    integer                      :: nb, nb_id, nb_dim, lo(NDIM), hi(NDIM)

#if NDIM == 2
    cc(0:nc+1, 0:nc+1) = boxes(id)%cc(0:nc+1, 0:nc+1, iv)
#elif NDIM == 3
    cc(0:nc+1, 0:nc+1, 0:nc+1) = boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iv)
#endif

    do nb = 1, af_num_neighbors
       nb_id = boxes(id)%neighbors(nb)

       if (nb_id > af_no_box) then
          call sides2_from_nb(boxes(nb_id), nb, iv, gc, nc)
       else if (nb_id == af_no_box) then
          call subr_rb(boxes, id, nb, iv, gc, nc)
       else
          call subr_bc(boxes, id, nb, iv, gc, nc)
       end if

       ! Determine ghost cell indices
       nb_dim = af_neighb_dim(nb)
       lo(:) = 1
       hi(:) = boxes(id)%n_cell
       lo(nb_dim) = -1 + af_neighb_high_01(nb) * (boxes(id)%n_cell + 3)
       hi(nb_dim) = lo(nb_dim)

#if NDIM == 2
       cc(lo(1):hi(1), lo(2):hi(2)) = reshape(gc, 1 + hi - lo)
#elif NDIM == 3
       cc(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)) = reshape(gc, 1 + hi - lo)
#endif
    end do
  end subroutine af_gc2_box

  !> Fill values on the side of a box from a neighbor nb
  subroutine sides2_from_nb(box_nb, nb, iv, gc_side, nc)
    type(box_t), intent(in) :: box_nb !< Neighbouring box
    integer, intent(in)       :: nb     !< Ghost cell / neighbor direction
    integer, intent(in)       :: iv     !< Ghost cell variable
    integer, intent(in)       :: nc
#if NDIM == 2
    real(dp), intent(out)     :: gc_side(nc)
#elif NDIM == 3
    real(dp), intent(out)     :: gc_side(nc, nc)
#endif

    select case (nb)
#if NDIM == 2
    case (af_neighb_lowx)
       gc_side = box_nb%cc(nc-1, 1:nc, iv)
    case (af_neighb_highx)
       gc_side = box_nb%cc(2, 1:nc, iv)
    case (af_neighb_lowy)
       gc_side = box_nb%cc(1:nc, nc-1, iv)
    case (af_neighb_highy)
       gc_side = box_nb%cc(1:nc, 2, iv)
#elif NDIM == 3
    case (af_neighb_lowx)
       gc_side = box_nb%cc(nc-1, 1:nc, 1:nc, iv)
    case (af_neighb_highx)
       gc_side = box_nb%cc(2, 1:nc, 1:nc, iv)
    case (af_neighb_lowy)
       gc_side = box_nb%cc(1:nc, nc-1, 1:nc, iv)
    case (af_neighb_highy)
       gc_side = box_nb%cc(1:nc, 2, 1:nc, iv)
    case (af_neighb_lowz)
       gc_side = box_nb%cc(1:nc, 1:nc, nc-1, iv)
    case (af_neighb_highz)
       gc_side = box_nb%cc(1:nc, 1:nc, 2, iv)
#endif
    end select
  end subroutine sides2_from_nb

  !> Linear interpolation (using data from neighbor) to fill ghost cells
  subroutine af_gc2_prolong_linear(boxes, id, nb, iv, gc_side, nc)
    type(box_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer, intent(in)         :: nc        !< Box n_cell
#if NDIM == 2
    real(dp), intent(out)       :: gc_side(nc) !< Ghost cells on side
#elif NDIM == 3
    real(dp), intent(out)       :: gc_side(nc, nc) !< Ghost cells on side
#endif
    integer                     :: ix, i, j
    integer                     :: i_c1, i_c2, j_c1, j_c2, p_nb_id
    integer                     :: p_id, ix_offset(NDIM)
#if NDIM == 3
    integer                     :: k, k_c1, k_c2
#endif

    p_id      = boxes(id)%parent
    p_nb_id   = boxes(p_id)%neighbors(nb)
    ix_offset = af_get_child_offset(boxes(id), nb)
    ix        = af_neighb_high_01(nb) * (nc+3) - 1 ! -1 or nc+2

    select case (af_neighb_dim(nb))
#if NDIM == 2
    case (1)
       i_c1 = ix_offset(1) + ishft(ix+1, -1) ! (ix+1)/2
       i_c2 = i_c1 + 1 - 2 * iand(ix, 1)     ! even: +1, odd: -1
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          gc_side(j) = &
               0.5_dp * boxes(p_nb_id)%cc(i_c1, j_c1, iv) + &
               0.25_dp * boxes(p_nb_id)%cc(i_c2, j_c1, iv) + &
               0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c2, iv)
       end do
    case (2)
       j_c1 = ix_offset(2) + ishft(ix+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(ix, 1)     ! even: +1, odd: -1
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
          gc_side(i) = &
               0.5_dp * boxes(p_nb_id)%cc(i_c1, j_c1, iv) + &
               0.25_dp * boxes(p_nb_id)%cc(i_c2, j_c1, iv) + &
               0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c2, iv)
       end do
#elif NDIM==3
    case (1)
       i_c1 = ix_offset(1) + ishft(ix+1, -1) ! (ix+1)/2
       i_c2 = i_c1 + 1 - 2 * iand(ix, 1)     ! even: +1, odd: -1
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
          do j = 1, nc
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
             gc_side(j, k) = &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c2, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c2, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c2, iv)
          end do
       end do
    case (2)
       j_c1 = ix_offset(2) + ishft(ix+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(ix, 1)     ! even: +1, odd: -1
       do k = 1, nc
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
             gc_side(i, k) = &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c2, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c2, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c2, iv)
          end do
       end do
    case (3)
       k_c1 = ix_offset(3) + ishft(ix+1, -1) ! (k+1)/2
       k_c2 = k_c1 + 1 - 2 * iand(ix, 1)     ! even: +1, odd: -1
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1
             gc_side(i, j) = &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c2, j_c1, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c2, k_c1, iv) + &
                  0.25_dp * boxes(p_nb_id)%cc(i_c1, j_c1, k_c2, iv)
          end do
       end do
#endif
    end select

  end subroutine af_gc2_prolong_linear

  ! This fills second ghost cells near physical boundaries using Neumann zero
  subroutine af_bc2_neumann_zero(boxes, id, nb, iv, gc_side, nc)
    type(box_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer, intent(in)         :: nc        !< Box n_cell
#if NDIM == 2
    real(dp), intent(out)       :: gc_side(nc) !< Ghost cells on side
#elif NDIM == 3
    real(dp), intent(out)       :: gc_side(nc, nc) !< Ghost cells on side
#endif

    select case (nb)
#if NDIM == 2
    case (af_neighb_lowx)
       gc_side = boxes(id)%cc(2, 1:nc, iv)
    case (af_neighb_highx)
       gc_side = boxes(id)%cc(nc-1, 1:nc, iv)
    case (af_neighb_lowy)
       gc_side = boxes(id)%cc(1:nc, 2, iv)
    case (af_neighb_highy)
       gc_side = boxes(id)%cc(1:nc, nc-1, iv)
#elif NDIM == 3
    case (af_neighb_lowx)
       gc_side = boxes(id)%cc(2, 1:nc, 1:nc, iv)
    case (af_neighb_highx)
       gc_side = boxes(id)%cc(nc-1, 1:nc, 1:nc, iv)
    case (af_neighb_lowy)
       gc_side = boxes(id)%cc(1:nc, 2, 1:nc, iv)
    case (af_neighb_highy)
       gc_side = boxes(id)%cc(1:nc, nc-1, 1:nc, iv)
    case (af_neighb_lowz)
       gc_side = boxes(id)%cc(1:nc, 1:nc, 2, iv)
    case (af_neighb_highz)
       gc_side = boxes(id)%cc(1:nc, 1:nc, nc-1, iv)
#endif
    end select
  end subroutine af_bc2_neumann_zero

  ! This fills second ghost cells near physical boundaries using Neumann zero
  subroutine af_bc2_dirichlet_zero(boxes, id, nb, iv, gc_side, nc)
    type(box_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv, nc
#if NDIM == 2
    real(dp), intent(out)       :: gc_side(nc)
#elif NDIM == 3
    real(dp), intent(out)       :: gc_side(nc, nc)
#endif

    select case (nb)
#if NDIM == 2
    case (af_neighb_lowx)
       gc_side = -boxes(id)%cc(2, 1:nc, iv)
    case (af_neighb_highx)
       gc_side = -boxes(id)%cc(nc-1, 1:nc, iv)
    case (af_neighb_lowy)
       gc_side = -boxes(id)%cc(1:nc, 2, iv)
    case (af_neighb_highy)
       gc_side = -boxes(id)%cc(1:nc, nc-1, iv)
#elif NDIM == 3
    case (af_neighb_lowx)
       gc_side = -boxes(id)%cc(2, 1:nc, 1:nc, iv)
    case (af_neighb_highx)
       gc_side = -boxes(id)%cc(nc-1, 1:nc, 1:nc, iv)
    case (af_neighb_lowy)
       gc_side = -boxes(id)%cc(1:nc, 2, 1:nc, iv)
    case (af_neighb_highy)
       gc_side = -boxes(id)%cc(1:nc, nc-1, 1:nc, iv)
    case (af_neighb_lowz)
       gc_side = -boxes(id)%cc(1:nc, 1:nc, 2, iv)
    case (af_neighb_highz)
       gc_side = -boxes(id)%cc(1:nc, 1:nc, nc-1, iv)
#endif
    end select
  end subroutine af_bc2_dirichlet_zero

  !> This fills corner ghost cells using linear extrapolation. The ghost cells
  !> on the sides already need to be filled.
  subroutine af_corner_gc_extrap(box, ix, iv)
    type(box_t), intent(inout) :: box !< Box to fill ghost cells for
    integer, intent(in)          :: ix(NDIM) !< Cell-index of corner
    integer, intent(in)          :: iv     !< Variable to fill
    integer                      :: di(NDIM)

    di = 1 - 2 * iand(ix, 1)    ! 0 -> di = 1, nc+1 -> di = -1

#if NDIM == 2
    box%cc(ix(1), ix(2), iv) = box%cc(ix(1)+di(1), ix(2), iv) &
         + box%cc(ix(1), ix(2)+di(2), iv) &
         - box%cc(ix(1)+di(1), ix(2)+di(2), iv)
#elif NDIM == 3
    box%cc(ix(1), ix(2), ix(3), iv) = &
         box%cc(ix(1), ix(2)+di(2), ix(3)+di(3), iv) + &
         box%cc(ix(1)+di(1), ix(2), ix(3)+di(3), iv) + &
         box%cc(ix(1)+di(1), ix(2)+di(2), ix(3), iv) - 2 * &
         box%cc(ix(1)+di(1), ix(2)+di(2), ix(3)+di(3), iv)
#endif
  end subroutine af_corner_gc_extrap

#if NDIM == 3
  !> This fills edge ghost cells using linear extrapolation. The ghost cells on
  !> the sides already need to be filled. This routine basically performs the
  !> same operation as af_corner_gc_extrap does in 2D.
  subroutine af_edge_gc_extrap(box, lo, dim, iv)
    type(box_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)          :: lo(NDIM) !< Lowest index of edge ghost cells
    integer, intent(in)          :: dim !< Dimension parallel to edge
    integer, intent(in)          :: iv !< Variable to fill
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

       box%cc(ix(1), ix(2), ix(3), iv) = &
            box%cc(ia(1), ia(2), ia(3), iv) + &
            box%cc(ib(1), ib(2), ib(3), iv) - &
            box%cc(ic(1), ic(2), ic(3), iv)
    end do

  end subroutine af_edge_gc_extrap
#endif

end module m_af_ghostcell
