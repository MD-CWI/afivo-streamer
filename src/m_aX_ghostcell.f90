!> This module contains routines related to the filling of ghost cells. Note that
!> corner ghost cells are not used in Afivo.
module m_a$D_ghostcell
  use m_a$D_types

  implicit none
  private

  public :: a$D_gc_tree
  public :: a$D_gc_ids
  public :: a$D_gc_box
  public :: a$D_bc_dirichlet_zero
  public :: a$D_bc_neumann_zero
  public :: a$D_gc_interp
  public :: a$D_gc_prolong_copy
  public :: a$D_gc_interp_lim
  public :: a$D_gc2_box
  public :: a$D_gc2_prolong_linear
  public :: a$D_bc2_neumann_zero
  public :: a$D_bc2_dirichlet_zero

contains

  !> Fill ghost cells for variables iv on the sides of all boxes, using
  !> subr_rb on refinement boundaries and subr_bc on physical boundaries
  subroutine a$D_gc_tree(tree, iv, subr_rb, subr_bc)
    type(a$D_t), intent(inout) :: tree    !< Tree to fill ghost cells on
    integer, intent(in)        :: iv      !< Variable for which ghost cells are set
    procedure(a$D_subr_rb)     :: subr_rb !< Procedure called at refinement boundaries
    procedure(a$D_subr_bc)     :: subr_bc !< Procedure called at physical boundaries
    integer                    :: lvl

    if (.not. tree%ready) stop "Tree not ready"

    do lvl = lbound(tree%lvls, 1), tree%highest_lvl
       call a$D_gc_ids(tree%boxes, tree%lvls(lvl)%ids, iv, subr_rb, subr_bc)
    end do
  end subroutine a$D_gc_tree

  !> Fill ghost cells for variables iv on the sides of all boxes, using subr_rb
  !> on refinement boundaries and subr_bc on physical boundaries. This routine
  !> assumes that ghost cells on other ids have been set already.
  subroutine a$D_gc_ids(boxes, ids, iv, subr_rb, subr_bc)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all the boxes
    integer, intent(in)          :: ids(:)   !< Ids of boxes for which we set ghost cells
    integer, intent(in)          :: iv       !< Variable for which ghost cells are set
    procedure(a$D_subr_rb)       :: subr_rb  !< Procedure called at refinement boundaries
    procedure(a$D_subr_bc)       :: subr_bc  !< Procedure called at physical boundaries
    integer                      :: i

    !$omp do
    do i = 1, size(ids)
       call a$D_gc_box(boxes, ids(i), iv, subr_rb, subr_bc)
    end do
    !$omp end do
  end subroutine a$D_gc_ids

  !> Fill ghost cells for variable iv
  subroutine a$D_gc_box(boxes, id, iv, subr_rb, subr_bc)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all the boxes
    integer, intent(in)          :: id       !< Id of box for which we set ghost cells
    integer, intent(in)          :: iv       !< Variable for which ghost cells are set
    procedure(a$D_subr_rb)       :: subr_rb  !< Procedure called at refinement boundaries
    procedure(a$D_subr_bc)       :: subr_bc  !< Procedure called at physical boundaries

    call a$D_gc_box_sides(boxes, id, iv, subr_rb, subr_bc)
    call a$D_gc_box_corner(boxes, id, iv)
  end subroutine a$D_gc_box

  !> Fill ghost cells for variable iv on the sides of a box, using subr_rb on
  !> refinement boundaries and subr_bc on physical boundaries.
  subroutine a$D_gc_box_sides(boxes, id, iv, subr_rb, subr_bc)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all the boxes
    integer, intent(in)          :: id       !< Id of box for which we set ghost cells
    integer, intent(in)          :: iv       !< Variable for which ghost cells are set
    procedure(a$D_subr_rb)       :: subr_rb  !< Procedure called at refinement boundaries
    procedure(a$D_subr_bc)       :: subr_bc  !< Procedure called at physical boundaries
    integer                      :: nb, nb_id, bc_type
    integer                      :: nb_dim, lo($D), hi($D), dnb($D)

    do nb = 1, a$D_num_neighbors
       nb_id = boxes(id)%neighbors(nb)

       if (nb_id > af_no_box) then
          ! Compute index range
          nb_dim = a$D_neighb_dim(nb)
          lo(:) = 1
          hi(:) = boxes(id)%n_cell
          lo(nb_dim) = a$D_neighb_high_01(nb) * (boxes(id)%n_cell + 1)
          hi(nb_dim) = lo(nb_dim)

          dnb = a$D_neighb_offset([nb])
          call copy_from_nb(boxes(id), boxes(nb_id), dnb, lo, hi, iv)
       else if (nb_id == af_no_box) then
          call subr_rb(boxes, id, nb, iv)
       else
          call subr_bc(boxes(id), nb, iv, bc_type)
          call bc_to_gc(boxes(id), nb, iv, bc_type)
       end if
    end do

  end subroutine a$D_gc_box_sides

  !> Fill corner ghost cells for variable iv on corners/edges of a box. If there
  !> is no box to copy the data from, use linear extrapolation. This routine
  !> assumes ghost cells on the sides of the box are available.
  subroutine a$D_gc_box_corner(boxes, id, iv)
    type(box$D_t), intent(inout)  :: boxes(:)              !< List of all the boxes
    integer, intent(in)          :: id                    !< Id of box for which we set ghost cells
    integer, intent(in)          :: iv                    !< Variable for which ghost cells are set
    integer                      :: n, nb_id, dnb($D), lo($D)
#if $D == 3
    integer                      :: hi($D), dim
#endif

#if $D == 3
    ! Have to do edges before corners (since extrapolation can depend on edge values)
    do n = 1, a3_num_edges
       dim = a3_edge_dim(n)

       ! Check whether there is a neighbor, and find its index
       nb_id = a$D_diag_neighb_id(boxes, id, a3_nb_adj_edge(:, n))

       lo = a3_edge_min_ix(:, n) * (boxes(id)%n_cell + 1)
       lo(dim) = 1

       if (nb_id > af_no_box) then
          hi      = lo
          hi(dim) = boxes(id)%n_cell
          dnb   = a$D_neighb_offset(a$D_nb_adj_edge(:, n))
          call copy_from_nb(boxes(id), boxes(nb_id), dnb, lo, hi, iv)
       else
          call a3_edge_gc_extrap(boxes(id), lo, dim, iv)
       end if
    end do
#endif

    do n = 1, a$D_num_children
       ! Check whether there is a neighbor, and find its index
       nb_id = a$D_diag_neighb_id(boxes, id, a$D_nb_adj_child(:, n))
       lo    = a$D_child_dix(:, n) * (boxes(id)%n_cell + 1)

       if (nb_id > af_no_box) then
          dnb = 2 * a$D_child_dix(:, n) - 1
          call copy_from_nb(boxes(id), boxes(nb_id), dnb, lo, lo, iv)
       else
          call a$D_corner_gc_extrap(boxes(id), lo, iv)
       end if
    end do
  end subroutine a$D_gc_box_corner

  subroutine bc_to_gc(box, nb, iv, bc_type)
    type(box$D_t), intent(inout)  :: box
    integer, intent(in)          :: iv                    !< Variable to fill
    integer, intent(in)          :: nb                    !< Neighbor direction
    integer, intent(in)          :: bc_type               !< Type of b.c.
    real(dp)                     :: c1, c2
    integer                      :: nc

    nc = box%n_cell

    ! If we call the interior point x, and the ghost point y, then a Dirichlet
    ! boundary value b can be imposed as:
    ! y = -x + 2*b
    ! A Neumann b.c. can be imposed as y = x +/- dx * b
    ! Below, we set coefficients c1 and c2 to handle both cases
    select case (bc_type)
    case (af_bc_dirichlet)
       c1 = -1
       c2 = 2
    case (af_bc_neumann)
       c1 = 1
       c2 = box%dr * a$D_neighb_high_pm(nb) ! This gives a + or - sign
    case default
       stop "fill_bc: unknown boundary condition"
    end select

    select case (nb)
#if $D == 2
    case (a2_neighb_lowx)
       box%cc(0, 1:nc, iv) = &
            c2 * box%cc(0, 1:nc, iv) + c1 * box%cc(1, 1:nc, iv)
    case (a2_neighb_highx)
       box%cc(nc+1, 1:nc, iv) = &
            c2 * box%cc(nc+1, 1:nc, iv) + c1 * box%cc(nc, 1:nc, iv)
    case (a2_neighb_lowy)
       box%cc(1:nc, 0, iv) = &
            c2 * box%cc(1:nc, 0, iv) + c1 * box%cc(1:nc, 1, iv)
    case (a2_neighb_highy)
       box%cc(1:nc, nc+1, iv) = &
            c2 * box%cc(1:nc, nc+1, iv) + c1 * box%cc(1:nc, nc, iv)
#elif $D == 3
    case (a3_neighb_lowx)
       box%cc(0, 1:nc, 1:nc, iv) = &
            c2 * box%cc(0, 1:nc, 1:nc, iv) + c1 * box%cc(1, 1:nc, 1:nc, iv)
    case (a3_neighb_highx)
       box%cc(nc+1, 1:nc, 1:nc, iv) = &
            c2 * box%cc(nc+1, 1:nc, 1:nc, iv) + c1 * box%cc(nc, 1:nc, 1:nc, iv)
    case (a3_neighb_lowy)
       box%cc(1:nc, 0, 1:nc, iv) = &
            c2 * box%cc(1:nc, 0, 1:nc, iv) + c1 * box%cc(1:nc, 1, 1:nc, iv)
    case (a3_neighb_highy)
       box%cc(1:nc, nc+1, 1:nc, iv) = &
            c2 * box%cc(1:nc, nc+1, 1:nc, iv) + c1 * box%cc(1:nc, nc, 1:nc, iv)
    case (a3_neighb_lowz)
       box%cc(1:nc, 1:nc, 0, iv) = &
            c2 * box%cc(1:nc, 1:nc, 0, iv) + c1 * box%cc(1:nc, 1:nc, 1, iv)
    case (a3_neighb_highz)
       box%cc(1:nc, 1:nc, nc+1, iv) = &
            c2 * box%cc(1:nc, 1:nc, nc+1, iv) + c1 * box%cc(1:nc, 1:nc, nc, iv)
#endif
    end select
  end subroutine bc_to_gc

  !> Partial prolongation to the ghost cells of box id from parent
  subroutine a$D_gc_prolong_copy(boxes, id, nb, iv)
    use m_a$D_prolong, only: a$D_prolong_copy
    type(box$D_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Id of child
    integer, intent(in)           :: iv       !< Variable to fill
    integer, intent(in)           :: nb       !< Neighbor to get data from
    integer                       :: p_id, nb_dim, lo($D), hi($D)

    p_id       = boxes(id)%parent
    nb_dim     = a$D_neighb_dim(nb)
    lo(:)      = 1
    hi(:)      = boxes(id)%n_cell
    lo(nb_dim) = a$D_neighb_high_01(nb) * (boxes(id)%n_cell+1)
    hi(nb_dim) = a$D_neighb_high_01(nb) * (boxes(id)%n_cell+1)

    call a$D_prolong_copy(boxes(p_id), boxes(id), iv, low=lo, high=hi)
  end subroutine a$D_gc_prolong_copy

  !> Interpolate between fine points and coarse neighbors to fill ghost cells
  !> near refinement boundaries
  subroutine a$D_gc_interp(boxes, id, nb, iv)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc, ix, ix_c, ix_f, i, j
    integer                     :: i_c1, i_c2, j_c1, j_c2, p_nb_id
    integer                     :: p_id, ix_offset($D)
    real(dp), parameter         :: sixth=1/6.0_dp, third=1/3.0_dp
#if $D == 3
    integer                     :: k_c1, k_c2, k
#endif

    nc        = boxes(id)%n_cell
    p_id      = boxes(id)%parent
    p_nb_id   = boxes(p_id)%neighbors(nb)
    ix_offset = a$D_get_child_offset(boxes(id), nb)

    if (a$D_neighb_low(nb)) then
       ix = 0
       ix_f = 1
       ix_c = nc
    else
       ix = nc+1
       ix_f = nc
       ix_c = 1
    end if

    select case (a$D_neighb_dim(nb))
#if $D == 2
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
#elif $D==3
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

  end subroutine a$D_gc_interp

  !> Interpolate between fine points and coarse neighbors to fill ghost cells
  !> near refinement boundaries. The ghost values are less than twice the coarse
  !> values.
  subroutine a$D_gc_interp_lim(boxes, id, nb, iv)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc, ix, ix_c, ix_f, i, j
    integer                     :: i_c1, i_c2, j_c1, j_c2, p_nb_id
    integer                     :: p_id, ix_offset($D)
    real(dp)                    :: c1, c2
    real(dp), parameter         :: sixth=1/6.0_dp, third=1/3.0_dp
#if $D == 3
    integer                     :: k_c1, k_c2, k
    real(dp)                    :: c3
#endif

    nc        = boxes(id)%n_cell
    p_id      = boxes(id)%parent
    p_nb_id   = boxes(p_id)%neighbors(nb)
    ix_offset = a$D_get_child_offset(boxes(id), nb)

    if (a$D_neighb_low(nb)) then
       ix = 0
       ix_f = 1
       ix_c = nc
    else
       ix = nc+1
       ix_f = nc
       ix_c = 1
    end if

    select case (a$D_neighb_dim(nb))
#if $D == 2
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
#elif $D==3
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

  end subroutine a$D_gc_interp_lim

  ! This fills ghost cells near physical boundaries using Neumann zero
  subroutine a$D_bc_neumann_zero(box, nb, iv, bc_type)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: nb, iv
    integer, intent(out)        :: bc_type

    bc_type = af_bc_neumann
    call a$D_set_box_gc(box, nb, iv, 0.0_dp)
  end subroutine a$D_bc_neumann_zero

  ! This fills ghost cells near physical boundaries using Neumann zero
  subroutine a$D_bc_dirichlet_zero(box, nb, iv, bc_type)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: nb, iv
    integer, intent(out)        :: bc_type

    bc_type = af_bc_dirichlet
    call a$D_set_box_gc(box, nb, iv, 0.0_dp)
  end subroutine a$D_bc_dirichlet_zero

  subroutine copy_from_nb(box, box_nb, dnb, lo, hi, iv)
    type(box$D_t), intent(inout) :: box     !< Box on which to fill ghost cells
    type(box$D_t), intent(in)    :: box_nb  !< Neighbouring box
    integer, intent(in)          :: dnb($D) !< Neighbor spatial index offset
    integer, intent(in)          :: lo($D)  !< Ghost cell low index
    integer, intent(in)          :: hi($D)  !< Ghost cell high index
    integer, intent(in)          :: iv      !< Ghost cell variable
    integer                      :: nlo($D), nhi($D)

    ! Get indices on neighbor
    nlo = lo - dnb * box%n_cell
    nhi = hi - dnb * box%n_cell

#if $D == 2
    box%cc(lo(1):hi(1), lo(2):hi(2), iv) = &
         box_nb%cc(nlo(1):nhi(1), nlo(2):nhi(2), iv)
#elif $D == 3
    box%cc(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), iv) = &
         box_nb%cc(nlo(1):nhi(1), nlo(2):nhi(2), nlo(3):nhi(3), iv)
#endif
  end subroutine copy_from_nb

  !> Get a second layer of ghost cell data (the 'normal' routines give just one
  !> layer of ghost cells). Use subr_rb on refinement boundaries and subr_bc
  !> on physical boundaries.
  subroutine a$D_gc2_box(boxes, id, iv, subr_rb, subr_bc, cc, nc)
    type(box$D_t), intent(inout) :: boxes(:)        !< List of all the boxes
    integer, intent(in)          :: id              !< Id of box for which we set ghost cells
    integer, intent(in)          :: iv              !< Variable for which ghost cells are set
    procedure(a$D_subr_egc)      :: subr_rb         !< Procedure called at refinement boundaries
    procedure(a$D_subr_egc)      :: subr_bc         !< Procedure called at physical boundaries
    integer, intent(in)          :: nc              !< box%n_cell
    !> The data with extra ghost cells
#if $D   == 2
    real(dp), intent(out)        :: cc(-1:nc+2, -1:nc+2)
#elif $D == 3
    real(dp), intent(out)        :: cc(-1:nc+2, -1:nc+2, -1:nc+2)
#endif
    integer                      :: nb, nb_id, nb_dim, lo($D), hi($D)

#if $D == 2
    cc(0:nc+1, 0:nc+1) = boxes(id)%cc(0:nc+1, 0:nc+1, iv)
#elif $D == 3
    cc(0:nc+1, 0:nc+1, 0:nc+1) = boxes(id)%cc(0:nc+1, 0:nc+1, 0:nc+1, iv)
#endif

    do nb = 1, a$D_num_neighbors
       nb_dim = a$D_neighb_dim(nb)
       lo(:) = 1
       hi(:) = boxes(id)%n_cell
       lo(nb_dim) = -1 + a$D_neighb_high_01(nb) * (boxes(id)%n_cell + 3)
       hi(nb_dim) = lo(nb_dim)

       nb_id = boxes(id)%neighbors(nb)

       if (nb_id > af_no_box) then
#if $D == 2
          call sides2_from_nb(boxes(nb_id), nb, iv, &
               cc(lo(1):hi(1), lo(2):hi(2)), nc)
#elif $D == 3
          call sides2_from_nb(boxes(nb_id), nb, iv, &
               cc(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)), nc)
#endif
       else if (nb_id == af_no_box) then
#if $D == 2
          call subr_rb(boxes, id, nb, iv, &
               cc(lo(1):hi(1), lo(2):hi(2)), nc)
#elif $D == 3
          call subr_rb(boxes, id, nb, iv, &
               cc(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)), nc)
#endif
       else
#if $D == 2
          call subr_bc(boxes, id, nb, iv, &
               cc(lo(1):hi(1), lo(2):hi(2)), nc)
#elif $D == 3
          call subr_bc(boxes, id, nb, iv, &
               cc(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3)), nc)
#endif
       end if
    end do
  end subroutine a$D_gc2_box

  !> Fill values on the side of a box from a neighbor nb
  subroutine sides2_from_nb(box_nb, nb, iv, gc_side, nc)
    type(box$D_t), intent(in) :: box_nb !< Neighbouring box
    integer, intent(in)       :: nb     !< Ghost cell / neighbor direction
    integer, intent(in)       :: iv     !< Ghost cell variable
    integer, intent(in)       :: nc
#if $D == 2
    real(dp), intent(out)     :: gc_side(nc)
#elif $D == 3
    real(dp), intent(out)     :: gc_side(nc, nc)
#endif

    select case (nb)
#if $D == 2
    case (a2_neighb_lowx)
       gc_side = box_nb%cc(nc-1, 1:nc, iv)
    case (a2_neighb_highx)
       gc_side = box_nb%cc(2, 1:nc, iv)
    case (a2_neighb_lowy)
       gc_side = box_nb%cc(1:nc, nc-1, iv)
    case (a2_neighb_highy)
       gc_side = box_nb%cc(1:nc, 2, iv)
#elif $D == 3
    case (a3_neighb_lowx)
       gc_side = box_nb%cc(nc-1, 1:nc, 1:nc, iv)
    case (a3_neighb_highx)
       gc_side = box_nb%cc(2, 1:nc, 1:nc, iv)
    case (a3_neighb_lowy)
       gc_side = box_nb%cc(1:nc, nc-1, 1:nc, iv)
    case (a3_neighb_highy)
       gc_side = box_nb%cc(1:nc, 2, 1:nc, iv)
    case (a3_neighb_lowz)
       gc_side = box_nb%cc(1:nc, 1:nc, nc-1, iv)
    case (a3_neighb_highz)
       gc_side = box_nb%cc(1:nc, 1:nc, 2, iv)
#endif
    end select
  end subroutine sides2_from_nb

  !> Linear interpolation (using data from neighbor) to fill ghost cells
  subroutine a$D_gc2_prolong_linear(boxes, id, nb, iv, gc_side, nc)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer, intent(in)         :: nc
#if $D == 2
    real(dp), intent(out)       :: gc_side(nc)
#elif $D == 3
    real(dp), intent(out)       :: gc_side(nc, nc)
#endif
    integer                     :: ix, i, j
    integer                     :: i_c1, i_c2, j_c1, j_c2, p_nb_id
    integer                     :: p_id, ix_offset($D)
#if $D == 3
    integer                     :: k, k_c1, k_c2
#endif

    p_id      = boxes(id)%parent
    p_nb_id   = boxes(p_id)%neighbors(nb)
    ix_offset = a$D_get_child_offset(boxes(id), nb)
    ix        = a$D_neighb_high_01(nb) * (nc+3) - 1 ! -1 or nc+2

    select case (a$D_neighb_dim(nb))
#if $D == 2
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
#elif $D==3
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

  end subroutine a$D_gc2_prolong_linear

  ! This fills second ghost cells near physical boundaries using Neumann zero
  subroutine a$D_bc2_neumann_zero(boxes, id, nb, iv, gc_side, nc)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv, nc
#if $D == 2
    real(dp), intent(out)       :: gc_side(nc)
#elif $D == 3
    real(dp), intent(out)       :: gc_side(nc, nc)
#endif

    select case (nb)
#if $D == 2
    case (a2_neighb_lowx)
       gc_side = boxes(id)%cc(2, 1:nc, iv)
    case (a2_neighb_highx)
       gc_side = boxes(id)%cc(nc-1, 1:nc, iv)
    case (a2_neighb_lowy)
       gc_side = boxes(id)%cc(1:nc, 2, iv)
    case (a2_neighb_highy)
       gc_side = boxes(id)%cc(1:nc, nc-1, iv)
#elif $D == 3
    case (a3_neighb_lowx)
       gc_side = boxes(id)%cc(2, 1:nc, 1:nc, iv)
    case (a3_neighb_highx)
       gc_side = boxes(id)%cc(nc-1, 1:nc, 1:nc, iv)
    case (a3_neighb_lowy)
       gc_side = boxes(id)%cc(1:nc, 2, 1:nc, iv)
    case (a3_neighb_highy)
       gc_side = boxes(id)%cc(1:nc, nc-1, 1:nc, iv)
    case (a3_neighb_lowz)
       gc_side = boxes(id)%cc(1:nc, 1:nc, 2, iv)
    case (a3_neighb_highz)
       gc_side = boxes(id)%cc(1:nc, 1:nc, nc-1, iv)
#endif
    end select
  end subroutine a$D_bc2_neumann_zero

  ! This fills second ghost cells near physical boundaries using Neumann zero
  subroutine a$D_bc2_dirichlet_zero(boxes, id, nb, iv, gc_side, nc)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv, nc
#if $D == 2
    real(dp), intent(out)       :: gc_side(nc)
#elif $D == 3
    real(dp), intent(out)       :: gc_side(nc, nc)
#endif

    select case (nb)
#if $D == 2
    case (a2_neighb_lowx)
       gc_side = -boxes(id)%cc(2, 1:nc, iv)
    case (a2_neighb_highx)
       gc_side = -boxes(id)%cc(nc-1, 1:nc, iv)
    case (a2_neighb_lowy)
       gc_side = -boxes(id)%cc(1:nc, 2, iv)
    case (a2_neighb_highy)
       gc_side = -boxes(id)%cc(1:nc, nc-1, iv)
#elif $D == 3
    case (a3_neighb_lowx)
       gc_side = -boxes(id)%cc(2, 1:nc, 1:nc, iv)
    case (a3_neighb_highx)
       gc_side = -boxes(id)%cc(nc-1, 1:nc, 1:nc, iv)
    case (a3_neighb_lowy)
       gc_side = -boxes(id)%cc(1:nc, 2, 1:nc, iv)
    case (a3_neighb_highy)
       gc_side = -boxes(id)%cc(1:nc, nc-1, 1:nc, iv)
    case (a3_neighb_lowz)
       gc_side = -boxes(id)%cc(1:nc, 1:nc, 2, iv)
    case (a3_neighb_highz)
       gc_side = -boxes(id)%cc(1:nc, 1:nc, nc-1, iv)
#endif
    end select
  end subroutine a$D_bc2_dirichlet_zero

  !> This fills corner ghost cells using linear extrapolation. The ghost cells
  !> on the sides already need to be filled.
  subroutine a$D_corner_gc_extrap(box, ix, iv)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)          :: ix($D)
    integer, intent(in)          :: iv
    integer                      :: di($D)

    di = 1 - 2 * iand(ix, 1)    ! 0 -> di = 1, nc+1 -> di = -1

#if $D == 2
    box%cc(ix(1), ix(2), iv) = box%cc(ix(1)+di(1), ix(2), iv) &
         + box%cc(ix(1), ix(2)+di(2), iv) &
         - box%cc(ix(1)+di(1), ix(2)+di(2), iv)
#elif $D == 3
    box%cc(ix(1), ix(2), ix(3), iv) = box%cc(ix(1)+di(1), ix(2), ix(3), iv) &
         + box%cc(ix(1), ix(2)+di(2), ix(3), iv) &
         + box%cc(ix(1), ix(2), ix(3)+di(3), iv) &
         - 2 * box%cc(ix(1)+di(1), ix(2)+di(2), ix(3)+di(3), iv)
#endif

  end subroutine a$D_corner_gc_extrap

#if $D == 3
  !> This fills edge ghost cells using linear extrapolation. The ghost cells on
  !> the sides already need to be filled. This routine basically performs the
  !> same operation as a$D_corner_gc_extrap does in 2D.
  subroutine a3_edge_gc_extrap(box, lo, dim, iv)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)          :: lo($D)
    integer, intent(in)          :: dim !< Dimension parallel to edge
    integer, intent(in)          :: iv
    integer                      :: di($D), ix($D), ia($D), ib($D), ic($D)
    integer                      :: n, o_dims($D-1)

    ! Dimensions other than/perpendicular to dim
    o_dims = [1 + mod(dim, $D), 1 + mod(dim + 1, $D)]

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

    ! print *, ix, di, dim, o_dims
    ! print *, lo

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

  end subroutine a3_edge_gc_extrap
#endif

end module m_a$D_ghostcell
