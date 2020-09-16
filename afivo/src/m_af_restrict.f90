#include "cpp_macros.h"
!> This module contains routines for restriction: going from fine to coarse
!> variables.
module m_af_restrict

  use m_af_types

  implicit none
  private

  public :: af_restrict_to_box
  public :: af_restrict_to_boxes
  public :: af_restrict_tree
  public :: af_restrict_box
  public :: af_restrict_ref_boundary
  ! public :: af_restrict_box_face

contains

  !> Restrict the children of a box to the box (e.g., in 2D, average the values
  !> at the four children to get the value for the parent)
  subroutine af_restrict_to_box(boxes, id, ivs)
    type(box_t), intent(inout) :: boxes(:) !< List of all the boxes
    integer, intent(in)        :: id       !< Box whose children will be restricted to it
    integer, intent(in)        :: ivs(:)   !< Variable to restrict
    integer                    :: nc, i_c, c_id

    nc = boxes(id)%n_cell
    do i_c = 1, af_num_children
       c_id = boxes(id)%children(i_c)
       if (c_id == af_no_box) cycle
       call af_restrict_box(boxes(c_id), boxes(id), ivs)
    end do
  end subroutine af_restrict_to_box

  !> Restrict the children of boxes ids(:) to them.
  subroutine af_restrict_to_boxes(boxes, ids, ivs)
    type(box_t), intent(inout) :: boxes(:) !< List of all the boxes
    integer, intent(in)        :: ids(:)   !< Boxes whose children will be restricted to it
    integer, intent(in)        :: ivs(:)   !< Variables to restrict
    integer                    :: i

    !$omp parallel do
    do i = 1, size(ids)
       call af_restrict_to_box(boxes, ids(i), ivs)
    end do
    !$omp end parallel do
  end subroutine af_restrict_to_boxes

  !> Restrict variables iv to all parent boxes, from the highest to the lowest level
  subroutine af_restrict_tree(tree, ivs)
    type(af_t), intent(inout) :: tree   !< Tree to restrict on
    integer, intent(in)       :: ivs(:) !< Variables to restrict
    integer                   :: lvl

    if (.not. tree%ready) stop "Tree not ready"
    do lvl = tree%highest_lvl-1, 1, -1
       call af_restrict_to_boxes(tree%boxes, tree%lvls(lvl)%parents, ivs)
    end do
  end subroutine af_restrict_tree

  !> Restriction of child box (box_c) to its parent (box_p)
  subroutine af_restrict_box(box_c, box_p, ivs, use_geometry)
    type(box_t), intent(in)       :: box_c        !< Child box to restrict
    type(box_t), intent(inout)    :: box_p        !< Parent box to restrict to
    integer, intent(in)           :: ivs(:)       !< Variable to restrict
    logical, intent(in), optional :: use_geometry !< If set to false, don't use geometry
    integer                       :: IJK, IJK_(c), IJK_(f), n, iv
    integer                       :: hnc, ix_offset(NDIM)
    logical                       :: use_geom
#if NDIM == 2
    real(dp)                      :: w1, w2
#endif

    hnc       = ishft(box_c%n_cell, -1) ! n_cell / 2
    ix_offset = af_get_child_offset(box_c)

    if (present(use_geometry)) then
       use_geom = use_geometry
    else
       use_geom = .true.
    end if

    do n = 1, size(ivs)
       iv = ivs(n)
#if NDIM == 1
       do i = 1, hnc
          i_c = ix_offset(1) + i
          i_f = 2 * i - 1
          box_p%cc(i_c, iv) = 0.5_dp * &
               sum(box_c%cc(i_f:i_f+1, iv))
       end do
#elif NDIM == 2
       if (box_p%coord_t == af_cyl .and. use_geom) then
          do j = 1, hnc
             j_c = ix_offset(2) + j
             j_f = 2 * j - 1
             do i = 1, hnc
                i_c = ix_offset(1) + i
                i_f = 2 * i - 1

                call af_cyl_child_weights(box_p, i_c, w1, w2)
                box_p%cc(i_c, j_c, iv) = 0.25_dp * (&
                     w1 * sum(box_c%cc(i_f, j_f:j_f+1, iv)) + &
                     w2 * sum(box_c%cc(i_f+1, j_f:j_f+1, iv)))
             end do
          end do
       else
          do j = 1, hnc
             j_c = ix_offset(2) + j
             j_f = 2 * j - 1
             do i = 1, hnc
                i_c = ix_offset(1) + i
                i_f = 2 * i - 1
                box_p%cc(i_c, j_c, iv) = 0.25_dp * &
                     sum(box_c%cc(i_f:i_f+1, j_f:j_f+1, iv))
             end do
          end do
       endif
#elif NDIM == 3
       do k = 1, hnc
          k_c = ix_offset(3) + k
          k_f = 2 * k - 1
          do j = 1, hnc
             j_c = ix_offset(2) + j
             j_f = 2 * j - 1
             do i = 1, hnc
                i_c = ix_offset(1) + i
                i_f = 2 * i - 1
                box_p%cc(i_c, j_c, k_c, iv) = 0.125_dp * &
                     sum(box_c%cc(i_f:i_f+1, j_f:j_f+1, k_f:k_f+1, iv))
             end do
          end do
       end do
#endif
    end do
  end subroutine af_restrict_box

  !> Restrict only next to refinement boundaries, which which can be required
  !> for filling coarse-grid ghost cells
  subroutine af_restrict_ref_boundary(tree, ivs)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: ivs(:)
    integer                   :: lvl, i, id, p_id

    !$omp parallel private(lvl, i, id, p_id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          p_id = tree%boxes(id)%parent

          ! Only restrict near refinement boundaries
          if (p_id > af_no_box .and. &
               any(tree%boxes(id)%neighbors == af_no_box)) then
             call af_restrict_box(tree%boxes(id), tree%boxes(p_id), ivs)
          end if
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine af_restrict_ref_boundary

!   !> Restriction of face-centered variables from child to parent
!   subroutine af_restrict_box_face(box_c, box_p, ivf, ivf_to)
!     type(box_t), intent(in)      :: box_c         !< Child box to restrict
!     type(box_t), intent(inout)   :: box_p         !< Parent box to restrict to
!     integer, intent(in)           :: ivf            !< Face-variable to restrict
!     integer, intent(in), optional :: ivf_to         !< Destination (if /= ivf)
!     integer                       :: i, j, i_f, j_f, i_c, j_c, i_dest
!     integer                       :: hnc, ix_offset(NDIM)
! #if NDIM == 3
!     integer                       :: k, k_f, k_c
! #endif

!     hnc       = ishft(box_c%n_cell, -1) ! n_cell / 2
!     ix_offset = af_get_child_offset(box_c)

!     if (present(ivf_to)) then
!        i_dest = ivf_to
!     else
!        i_dest = ivf
!     end if

!     if (box_p%coord_t == af_cyl) &
!          stop "restrict_box_face not implemented for cylindrical case"

! #if NDIM == 2
!     do j = 1, hnc
!        j_c = ix_offset(2) + j
!        j_f = 2 * j - 1
!        do i = 1, hnc
!           i_c = ix_offset(1) + i
!           i_f = 2 * i - 1

!           box_p%fc(i_c, j_c, 1, i_dest) = 0.5_dp * &
!                sum(box_c%fc(i_f, j_f:j_f+1, 1, ivf))
!           if (i == hnc) then
!              box_p%fc(i_c+1, j_c, 1, i_dest) = 0.5_dp * &
!                sum(box_c%fc(i_f+2, j_f:j_f+1, 1, ivf))
!           end if

!           box_p%fc(i_c, j_c, 2, i_dest) = 0.5_dp * &
!                sum(box_c%fc(i_f:i_f+1, j_f, 2, ivf))
!           if (j == hnc) then
!              box_p%fc(i_c, j_c+1, 2, i_dest) = 0.5_dp * &
!                sum(box_c%fc(i_f:i_f+1, j_f+2, 2, ivf))
!           end if

!        end do
!     end do
! #elif NDIM == 3
!     if (box_p%coord_t == af_cyl) &
!          stop "restrict_box_face not implemented for 3D case"
! #endif
!   end subroutine af_restrict_box_face

end module m_af_restrict
