! This module contains the routines related to prolongation or interpolation
! (going from a coarse to a fine variable).
!
! Author: Jannis Teunissen
! License: GPLv3

module m_a$D_prolong
  use m_a$D_t

  implicit none
  private

  public :: a$D_prolong0_from
  public :: a$D_prolong0_to
  public :: a$D_prolong1_from
  public :: a$D_prolong1_to
  public :: a$D_prolong2_from
  public :: a$D_prolong2_to

contains

    !> Zeroth-order prolongation to children.
  subroutine a$D_prolong0_from(boxes, id, iv)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Box whose children we will fill
    integer, intent(in)         :: iv        !< Variable that is filled
    integer                     :: i_c, c_id

    do i_c = 1, a$D_num_children
       c_id = boxes(id)%children(i_c)
       if (c_id == a5_no_box) cycle
       call a$D_prolong0_to(boxes, c_id, iv)
    end do
  end subroutine a$D_prolong0_from

  !> Partial prolongation to a child (from parent) using injection (simply copy value)
  subroutine a$D_prolong0_to(boxes, id, iv, lo_a, hi_a)
    use m_a$D_utils, only: a$D_get_child_offset
    type(box$D_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Id of child
    integer, intent(in)           :: iv       !< Variable to fill
    integer, intent(in), optional :: lo_a($D) !< Min cell index at child
    integer, intent(in), optional :: hi_a($D) !< Max cell index at child
    integer                       :: nc, p_id, ix_offset($D)
    integer                       :: i, j, i_c1, j_c1, lo($D), hi($D)
#if $D == 3
    integer                       :: k, k_c1
#endif

    nc   = boxes(id)%n_cell
    p_id = boxes(id)%parent
    lo   = 1; if (present(lo_a)) lo = lo_a
    hi   = nc; if (present(hi_a)) hi = hi_a

    ! Offset of child w.r.t. parent
    ix_offset = a$D_get_child_offset(boxes(id))

#if $D == 2
    do j = lo(2), hi(2)
       j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
       do i = lo(1), hi(1)
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          boxes(id)%cc(i, j, iv) = boxes(p_id)%cc(i_c1, j_c1, iv)
       end do
    end do
#elif $D == 3
    do k = lo(3), hi(3)
       k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
       do j = lo(2), hi(2)
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          do i = lo(1), hi(1)
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             boxes(id)%cc(i, j, k, iv) = boxes(p_id)%cc(i_c1, j_c1, k_c1, iv)
          end do
       end do
    end do
#endif
  end subroutine a$D_prolong0_to

  !> Linear prolongation to children. We use 2-1-1 interpolation (2d) and
  !> 1-1-1-1 interpolation (3D), which do not require corner ghost cells.
  subroutine a$D_prolong1_from(boxes, id, iv, add)
    type(box$D_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Box whose children we will fill
    integer, intent(in)           :: iv       !< Variable that is filled
    logical, intent(in), optional :: add      !< Add to old values
    integer                       :: i_c, c_id

    do i_c = 1, a$D_num_children
       c_id = boxes(id)%children(i_c)
       if (c_id == a5_no_box) cycle
       call a$D_prolong1_to(boxes, c_id, iv, add)
    end do
  end subroutine a$D_prolong1_from

  !> Prolongation to a child (from parent) using linear interpolation. We use
  !> 2-1-1 interpolation (2D) and 1-1-1-1 interpolation (3D) which do not need
  !> corner ghost cells.
  subroutine a$D_prolong1_to(boxes, id, iv, add)
    use m_a$D_utils, only: a$D_get_child_offset
    type(box$D_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Id of child
    integer, intent(in)           :: iv       !< Variable to fill
    logical, intent(in), optional :: add      !< Add to old values
    integer                       :: hnc, nc, p_id, ix_offset($D)
    integer                       :: i, j, i_c, i_f, j_c, j_f
    real(dp)                      :: f0, flx, fhx, fly, fhy
    logical                       :: add_to
#if $D == 3
    real(dp)                      :: flz, fhz
    integer                       :: k, k_c, k_f
#endif

    nc        = boxes(id)%n_cell
    hnc       = ishft(boxes(id)%n_cell, -1)
    p_id      = boxes(id)%parent
    ix_offset = a$D_get_child_offset(boxes(id))
    add_to    = .false.; if (present(add)) add_to = add

    if (.not. add_to) then
#if $D == 2
       boxes(id)%cc(1:nc, 1:nc, iv) = 0
#elif $D == 3
       boxes(id)%cc(1:nc, 1:nc, 1:nc, iv) = 0
#endif
    end if

#if $D == 2
    do j = 1, hnc
       j_c = j + ix_offset(2)
       j_f = 2 * j - 1
       do i = 1, hnc
          i_c = i + ix_offset(1)
          i_f = 2 * i - 1

          f0 = 0.5_dp * boxes(p_id)%cc(i_c, j_c, iv)
          flx = 0.25_dp * boxes(p_id)%cc(i_c-1, j_c, iv)
          fhx = 0.25_dp * boxes(p_id)%cc(i_c+1, j_c, iv)
          fly = 0.25_dp * boxes(p_id)%cc(i_c, j_c-1, iv)
          fhy = 0.25_dp * boxes(p_id)%cc(i_c, j_c+1, iv)

          boxes(id)%cc(i_f,   j_f,   iv) = f0 + flx + fly &
               + boxes(id)%cc(i_f,   j_f,   iv)
          boxes(id)%cc(i_f+1, j_f,   iv) = f0 + fhx + fly &
               + boxes(id)%cc(i_f+1, j_f,   iv)
          boxes(id)%cc(i_f,   j_f+1, iv) = f0 + flx + fhy &
               + boxes(id)%cc(i_f,   j_f+1, iv)
          boxes(id)%cc(i_f+1, j_f+1, iv) = f0 + fhx + fhy &
               + boxes(id)%cc(i_f+1, j_f+1, iv)
       end do
    end do
#elif $D == 3
    do k = 1, hnc
       k_c = k + ix_offset(3)
       k_f = 2 * k - 1
       do j = 1, hnc
          j_c = j + ix_offset(2)
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = i + ix_offset(1)
             i_f = 2 * i - 1

             f0  = 0.25_dp * boxes(p_id)%cc(i_c,   j_c,   k_c,   iv)
             flx = 0.25_dp * boxes(p_id)%cc(i_c-1, j_c,   k_c,   iv)
             fhx = 0.25_dp * boxes(p_id)%cc(i_c+1, j_c,   k_c,   iv)
             fly = 0.25_dp * boxes(p_id)%cc(i_c,   j_c-1, k_c,   iv)
             fhy = 0.25_dp * boxes(p_id)%cc(i_c,   j_c+1, k_c,   iv)
             flz = 0.25_dp * boxes(p_id)%cc(i_c,   j_c,   k_c-1, iv)
             fhz = 0.25_dp * boxes(p_id)%cc(i_c,   j_c,   k_c+1, iv)

             boxes(id)%cc(i_f,   j_f,   k_f,   iv) = f0 + flx + &
                  fly + flz + boxes(id)%cc(i_f,   j_f,   k_f,   iv)
             boxes(id)%cc(i_f+1, j_f,   k_f,   iv) = f0 + fhx + &
                  fly + flz + boxes(id)%cc(i_f+1, j_f,   k_f,   iv)
             boxes(id)%cc(i_f,   j_f+1, k_f,   iv) = f0 + flx + &
                  fhy + flz + boxes(id)%cc(i_f,   j_f+1, k_f,   iv)
             boxes(id)%cc(i_f+1, j_f+1, k_f,   iv) = f0 + fhx + &
                  fhy + flz + boxes(id)%cc(i_f+1, j_f+1, k_f,   iv)
             boxes(id)%cc(i_f,   j_f,   k_f+1, iv) = f0 + flx + &
                  fly + fhz + boxes(id)%cc(i_f,   j_f,   k_f+1, iv)
             boxes(id)%cc(i_f+1, j_f,   k_f+1, iv) = f0 + fhx + &
                  fly + fhz + boxes(id)%cc(i_f+1, j_f,   k_f+1, iv)
             boxes(id)%cc(i_f,   j_f+1, k_f+1, iv) = f0 + flx + &
                  fhy + fhz + boxes(id)%cc(i_f,   j_f+1, k_f+1, iv)
             boxes(id)%cc(i_f+1, j_f+1, k_f+1, iv) = f0 + fhx + &
                  fhy + fhz + boxes(id)%cc(i_f+1, j_f+1, k_f+1, iv)
          end do
       end do
    end do
#endif
  end subroutine a$D_prolong1_to

  !> Quadratic prolongation to children. We use stencils that do not require
  !> corner ghost cells.
  subroutine a$D_prolong2_from(boxes, id, iv)
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Box whose children we will fill
    integer, intent(in)         :: iv        !< Variable that is filled
    integer                     :: i_c, c_id

    do i_c = 1, a$D_num_children
       c_id = boxes(id)%children(i_c)
       if (c_id == a5_no_box) cycle
       call a$D_prolong2_to(boxes, c_id, iv)
    end do
  end subroutine a$D_prolong2_from

  !> Prolongation to a child (from parent) using quadratic interpolation. We use
  !> 5 / 7 point stencils which do not need corner ghost cells.
  !> @TODO 3D version
  subroutine a$D_prolong2_to(boxes, id, iv)
    use m_a$D_utils, only: a$D_get_child_offset
    type(box$D_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)          :: id       !< Id of child
    integer, intent(in)          :: iv       !< Variable to fill
    integer                      :: hnc, p_id, ix_offset($D)
    integer                      :: i, j
    integer                      :: i_c, i_f, j_c, j_f
    real(dp)                     :: f0, fx, fy, fxx, fyy, f2
#if $D == 3
    real(dp)                     :: fz, fzz
    integer                      :: k, k_c, k_f
#endif

    hnc       = ishft(boxes(id)%n_cell, -1)
    p_id      = boxes(id)%parent
    ix_offset = a$D_get_child_offset(boxes(id))

#if $D == 2
    do j = 1, hnc
       j_c = j + ix_offset(2)
       j_f = 2 * j - 1
       do i = 1, hnc
          i_c = i + ix_offset(1)
          i_f = 2 * i - 1

          f0 = boxes(p_id)%cc(i_c, j_c, iv)
          fx = 0.125_dp * (boxes(p_id)%cc(i_c+1, j_c, iv) - &
               boxes(p_id)%cc(i_c-1, j_c, iv))
          fy = 0.125_dp * (boxes(p_id)%cc(i_c, j_c+1, iv) - &
               boxes(p_id)%cc(i_c, j_c-1, iv))
          fxx = 0.03125_dp * (boxes(p_id)%cc(i_c-1, j_c, iv) - &
               2 * f0 + boxes(p_id)%cc(i_c+1, j_c, iv))
          fyy = 0.03125_dp * (boxes(p_id)%cc(i_c, j_c-1, iv) - &
               2 * f0 + boxes(p_id)%cc(i_c, j_c+1, iv))
          f2 = fxx + fyy

          boxes(id)%cc(i_f,   j_f,   iv) = f0 - fx - fy + f2
          boxes(id)%cc(i_f+1, j_f,   iv) = f0 + fx - fy + f2
          boxes(id)%cc(i_f,   j_f+1, iv) = f0 - fx + fy + f2
          boxes(id)%cc(i_f+1, j_f+1, iv) = f0 + fx + fy + f2
       end do
    end do
#elif $D == 3
    do k = 1, hnc
       k_c = k + ix_offset(3)
       k_f = 2 * k - 1
       do j = 1, hnc
          j_c = j + ix_offset(2)
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = i + ix_offset(1)
             i_f = 2 * i - 1

             f0 = boxes(p_id)%cc(i_c, j_c, k_c, iv)
             fx = 0.125_dp * (boxes(p_id)%cc(i_c+1, j_c, k_c, iv) - &
                  boxes(p_id)%cc(i_c-1, j_c, k_c, iv))
             fy = 0.125_dp * (boxes(p_id)%cc(i_c, j_c+1, k_c, iv) - &
                  boxes(p_id)%cc(i_c, j_c-1, k_c, iv))
             fz = 0.125_dp * (boxes(p_id)%cc(i_c, j_c, k_c+1, iv) - &
                  boxes(p_id)%cc(i_c, j_c, k_c-1, iv))
             fxx = 0.03125_dp * (boxes(p_id)%cc(i_c-1, j_c, k_c, iv) - &
                  2 * f0 + boxes(p_id)%cc(i_c+1, j_c, k_c, iv))
             fyy = 0.03125_dp * (boxes(p_id)%cc(i_c, j_c-1, k_c, iv) - &
                  2 * f0 + boxes(p_id)%cc(i_c, j_c+1, k_c, iv))
             fzz = 0.03125_dp * (boxes(p_id)%cc(i_c, j_c, k_c-1, iv) - &
                  2 * f0 + boxes(p_id)%cc(i_c, j_c, k_c+1, iv))
             f2 = fxx + fyy + fzz

             boxes(id)%cc(i_f,   j_f,   k_f,   iv) = f0 - fx - fy - fz + f2
             boxes(id)%cc(i_f+1, j_f,   k_f,   iv) = f0 + fx - fy - fz + f2
             boxes(id)%cc(i_f,   j_f+1, k_f,   iv) = f0 - fx + fy - fz + f2
             boxes(id)%cc(i_f+1, j_f+1, k_f,   iv) = f0 + fx + fy - fz + f2
             boxes(id)%cc(i_f,   j_f,   k_f+1, iv) = f0 - fx - fy + fz + f2
             boxes(id)%cc(i_f+1, j_f,   k_f+1, iv) = f0 + fx - fy + fz + f2
             boxes(id)%cc(i_f,   j_f+1, k_f+1, iv) = f0 - fx + fy + fz + f2
             boxes(id)%cc(i_f+1, j_f+1, k_f+1, iv) = f0 + fx + fy + fz + f2
          end do
       end do
    end do
#endif
  end subroutine a$D_prolong2_to

end module m_a$D_prolong
