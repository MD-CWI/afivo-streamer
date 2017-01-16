!> This module contains the routines related to prolongation: going from
!> coarse to fine variables.
module m_a$D_prolong
  use m_a$D_types

  implicit none
  private

  public :: a$D_prolong_copy_from
  public :: a$D_prolong_copy
  public :: a$D_prolong_linear_from
  public :: a$D_prolong_linear
  public :: a$D_prolong_quadratic_from
  public :: a$D_prolong_quadratic

contains

  !> Zeroth-order prolongation to children.
  subroutine a$D_prolong_copy_from(boxes, id, iv, iv_to, add)
    type(box$D_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Box whose children we will fill
    integer, intent(in)           :: iv       !< Variable that is prolonged
    integer, intent(in), optional :: iv_to    !< Destination variable
    logical, intent(in), optional :: add      !< Add to old values
    integer                       :: i_c, c_id

    do i_c = 1, a$D_num_children
       c_id = boxes(id)%children(i_c)
       if (c_id == af_no_box) cycle
       call a$D_prolong_copy(boxes(id), boxes(c_id), iv, iv_to=iv_to, add=add)
    end do
  end subroutine a$D_prolong_copy_from

  !> Partial prolongation to a child (from parent) using injection (simply copy value)
  subroutine a$D_prolong_copy(box_p, box_c, iv, iv_to, low, high, add)
    type(box$D_t), intent(in)      :: box_p !< Parent box
    type(box$D_t), intent(inout)   :: box_c !< Child box
    integer, intent(in)           :: iv       !< Variable to fill
    integer, intent(in), optional :: iv_to    !< Destination variable
    integer, intent(in), optional :: low($D) !< Min cell index at child
    integer, intent(in), optional :: high($D) !< Max cell index at child
    logical, intent(in), optional :: add      !< Add to old values
    logical                       :: add_to
    integer                       :: nc, ix_offset($D), ivc
    integer                       :: i, j, i_c1, j_c1, lo($D), hi($D)
#if $D == 3
    integer                       :: k, k_c1
#endif

    nc   = box_c%n_cell
    add_to = .false.; if (present(add)) add_to = add
    ivc = iv; if (present(iv_to)) ivc = iv_to
    lo   = 1; if (present(low)) lo = low
    hi   = nc; if (present(high)) hi = high

    ! Offset of child w.r.t. parent
    ix_offset = a$D_get_child_offset(box_c)

    if (add_to) then
#if $D == 2
       do j = lo(2), hi(2)
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          do i = lo(1), hi(1)
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             box_c%cc(i, j, ivc) = box_c%cc(i, j, ivc) + &
                  box_p%cc(i_c1, j_c1, iv)
          end do
       end do
#elif $D == 3
       do k = lo(3), hi(3)
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          do j = lo(2), hi(2)
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             do i = lo(1), hi(1)
                i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
                box_c%cc(i, j, k, ivc) = box_c%cc(i, j, k, ivc) + &
                     box_p%cc(i_c1, j_c1, k_c1, iv)
             end do
          end do
       end do
#endif
    else
#if $D == 2
       do j = lo(2), hi(2)
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          do i = lo(1), hi(1)
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             box_c%cc(i, j, ivc) = box_p%cc(i_c1, j_c1, iv)
          end do
       end do
#elif $D == 3
       do k = lo(3), hi(3)
          k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
          do j = lo(2), hi(2)
             j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
             do i = lo(1), hi(1)
                i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
                box_c%cc(i, j, k, ivc) = box_p%cc(i_c1, j_c1, k_c1, iv)
             end do
          end do
       end do
#endif
    end if
  end subroutine a$D_prolong_copy

  !> Linear prolongation to children. We use 2-1-1 interpolation (2d) and
  !> 1-1-1-1 interpolation (3D), which do not require corner ghost cells.
  subroutine a$D_prolong_linear_from(boxes, id, iv, iv_to, add)
    type(box$D_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Box whose children we will fill
    integer, intent(in)           :: iv       !< Variable that is filled
    integer, intent(in), optional :: iv_to    !< Destination variable
    logical, intent(in), optional :: add      !< Add to old values
    integer                       :: i_c, c_id

    do i_c = 1, a$D_num_children
       c_id = boxes(id)%children(i_c)
       if (c_id == af_no_box) cycle
       call a$D_prolong_linear(boxes(id), boxes(c_id), iv, iv_to, add)
    end do
  end subroutine a$D_prolong_linear_from

  !> Prolongation to a child (from parent) using linear interpolation. We use
  !> 2-1-1 interpolation (2D) and 1-1-1-1 interpolation (3D) which do not need
  !> corner ghost cells.
  subroutine a$D_prolong_linear(box_p, box_c, iv, iv_to, add)
    type(box$D_t), intent(in)      :: box_p !< Parent box
    type(box$D_t), intent(inout)   :: box_c !< Child box
    integer, intent(in)           :: iv       !< Variable to fill
    integer, intent(in), optional :: iv_to    !< Destination variable
    logical, intent(in), optional :: add      !< Add to old values
    integer                       :: hnc, nc, ix_offset($D), ivc
    integer                       :: i, j, i_c, i_f, j_c, j_f
    real(dp)                      :: f0, flx, fhx, fly, fhy
    logical                       :: add_to
#if $D == 3
    real(dp)                      :: flz, fhz
    integer                       :: k, k_c, k_f
#endif

    nc        = box_c%n_cell
    hnc       = ishft(box_c%n_cell, -1)
    ix_offset = a$D_get_child_offset(box_c)
    add_to    = .false.; if (present(add)) add_to = add
    ivc       = iv; if (present(iv_to)) ivc = iv_to

    if (.not. add_to) then
#if $D == 2
       box_c%cc(1:nc, 1:nc, ivc) = 0
#elif $D == 3
       box_c%cc(1:nc, 1:nc, 1:nc, ivc) = 0
#endif
    end if

#if $D == 2
    do j = 1, hnc
       j_c = j + ix_offset(2)
       j_f = 2 * j - 1
       do i = 1, hnc
          i_c = i + ix_offset(1)
          i_f = 2 * i - 1

          f0 = 0.5_dp * box_p%cc(i_c, j_c, iv)
          flx = 0.25_dp * box_p%cc(i_c-1, j_c, iv)
          fhx = 0.25_dp * box_p%cc(i_c+1, j_c, iv)
          fly = 0.25_dp * box_p%cc(i_c, j_c-1, iv)
          fhy = 0.25_dp * box_p%cc(i_c, j_c+1, iv)

          box_c%cc(i_f,   j_f,   ivc) = f0 + flx + fly &
               + box_c%cc(i_f,   j_f,   ivc)
          box_c%cc(i_f+1, j_f,   ivc) = f0 + fhx + fly &
               + box_c%cc(i_f+1, j_f,   ivc)
          box_c%cc(i_f,   j_f+1, ivc) = f0 + flx + fhy &
               + box_c%cc(i_f,   j_f+1, ivc)
          box_c%cc(i_f+1, j_f+1, ivc) = f0 + fhx + fhy &
               + box_c%cc(i_f+1, j_f+1, ivc)
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

             f0  = 0.25_dp * box_p%cc(i_c,   j_c,   k_c,   iv)
             flx = 0.25_dp * box_p%cc(i_c-1, j_c,   k_c,   iv)
             fhx = 0.25_dp * box_p%cc(i_c+1, j_c,   k_c,   iv)
             fly = 0.25_dp * box_p%cc(i_c,   j_c-1, k_c,   iv)
             fhy = 0.25_dp * box_p%cc(i_c,   j_c+1, k_c,   iv)
             flz = 0.25_dp * box_p%cc(i_c,   j_c,   k_c-1, iv)
             fhz = 0.25_dp * box_p%cc(i_c,   j_c,   k_c+1, iv)

             box_c%cc(i_f,   j_f,   k_f,   ivc) = f0 + flx + &
                  fly + flz + box_c%cc(i_f,   j_f,   k_f,   ivc)
             box_c%cc(i_f+1, j_f,   k_f,   ivc) = f0 + fhx + &
                  fly + flz + box_c%cc(i_f+1, j_f,   k_f,   ivc)
             box_c%cc(i_f,   j_f+1, k_f,   ivc) = f0 + flx + &
                  fhy + flz + box_c%cc(i_f,   j_f+1, k_f,   ivc)
             box_c%cc(i_f+1, j_f+1, k_f,   ivc) = f0 + fhx + &
                  fhy + flz + box_c%cc(i_f+1, j_f+1, k_f,   ivc)
             box_c%cc(i_f,   j_f,   k_f+1, ivc) = f0 + flx + &
                  fly + fhz + box_c%cc(i_f,   j_f,   k_f+1, ivc)
             box_c%cc(i_f+1, j_f,   k_f+1, ivc) = f0 + fhx + &
                  fly + fhz + box_c%cc(i_f+1, j_f,   k_f+1, ivc)
             box_c%cc(i_f,   j_f+1, k_f+1, ivc) = f0 + flx + &
                  fhy + fhz + box_c%cc(i_f,   j_f+1, k_f+1, ivc)
             box_c%cc(i_f+1, j_f+1, k_f+1, ivc) = f0 + fhx + &
                  fhy + fhz + box_c%cc(i_f+1, j_f+1, k_f+1, ivc)
          end do
       end do
    end do
#endif
  end subroutine a$D_prolong_linear

  !> Quadratic prolongation to children. We use stencils that do not require
  !> corner ghost cells.
  subroutine a$D_prolong_quadratic_from(boxes, id, iv, iv_to, add)
    type(box$D_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Box whose children we will fill
    integer, intent(in)           :: iv       !< Variable that is filled
    integer, intent(in), optional :: iv_to    !< Destination variable
    logical, intent(in), optional :: add      !< Add to old values
    integer                       :: i_c, c_id

    do i_c = 1, a$D_num_children
       c_id = boxes(id)%children(i_c)
       if (c_id == af_no_box) cycle
       call a$D_prolong_quadratic(boxes(id), boxes(c_id), iv, iv_to, add)
    end do
  end subroutine a$D_prolong_quadratic_from

  !> Prolongation to a child (from parent) using quadratic interpolation. The 3D
  !> version of this routine misses the cross-derivative term f_xyz, and the 2D
  !> version relies on corner ghost cells (which are filled with an
  !> extrapolation in a$D_gc_box).
  subroutine a$D_prolong_quadratic(box_p, box_c, iv, iv_to, add)
    type(box$D_t), intent(in)      :: box_p !< Parent box
    type(box$D_t), intent(inout)   :: box_c !< Child box
    integer, intent(in)           :: iv       !< Variable to fill
    integer, intent(in), optional :: iv_to    !< Destination variable
    logical, intent(in), optional :: add      !< Add to old values
    logical                      :: add_to
    integer                      :: hnc, ix_offset($D)
    integer                      :: i, j, nc, ivc
    integer                      :: i_c, i_f, j_c, j_f
    real(dp)                     :: f0, fx, fy, fxx, fyy, f2
#if $D == 2
    real(dp)                     :: fxy(2**$D)
#elif $D == 3
    real(dp)                     :: fz, fzz
    integer                      :: k, k_c, k_f
#endif

    nc        = box_c%n_cell
    hnc       = ishft(box_c%n_cell, -1)
    ix_offset = a$D_get_child_offset(box_c)
    add_to    = .false.; if (present(add)) add_to = add
    ivc       = iv; if (present(iv_to)) ivc = iv_to

    if (.not. add_to) then
#if $D == 2
       box_c%cc(1:nc, 1:nc, ivc) = 0
#elif $D == 3
       box_c%cc(1:nc, 1:nc, 1:nc, ivc) = 0
#endif
    end if

#if $D == 2
    do j = 1, hnc
       j_c = j + ix_offset(2)
       j_f = 2 * j - 1
       do i = 1, hnc
          i_c = i + ix_offset(1)
          i_f = 2 * i - 1

          f0 = box_p%cc(i_c, j_c, iv)
          fx = 0.125_dp * (box_p%cc(i_c+1, j_c, iv) - &
               box_p%cc(i_c-1, j_c, iv))
          fy = 0.125_dp * (box_p%cc(i_c, j_c+1, iv) - &
               box_p%cc(i_c, j_c-1, iv))
          fxx = 0.03125_dp * (box_p%cc(i_c-1, j_c, iv) - &
               2 * f0 + box_p%cc(i_c+1, j_c, iv))
          fyy = 0.03125_dp * (box_p%cc(i_c, j_c-1, iv) - &
               2 * f0 + box_p%cc(i_c, j_c+1, iv))
          f2 = fxx + fyy

          fxy(1) = 0.0625_dp * (box_p%cc(i_c-1, j_c-1, iv) + f0 - &
               box_p%cc(i_c-1, j_c, iv) - box_p%cc(i_c, j_c-1, iv))
          fxy(2) = 0.0625_dp * (box_p%cc(i_c+1, j_c-1, iv) + f0 - &
               box_p%cc(i_c+1, j_c, iv) - box_p%cc(i_c, j_c-1, iv))
          fxy(3) = 0.0625_dp * (box_p%cc(i_c-1, j_c+1, iv) + f0 - &
               box_p%cc(i_c-1, j_c, iv) - box_p%cc(i_c, j_c+1, iv))
          fxy(4) = 0.0625_dp * (box_p%cc(i_c+1, j_c+1, iv) + f0 - &
               box_p%cc(i_c+1, j_c, iv) - box_p%cc(i_c, j_c+1, iv))

          box_c%cc(i_f,   j_f,   ivc) = f0 - fx - fy + f2 + fxy(1) + &
               box_c%cc(i_f,   j_f,   ivc)
          box_c%cc(i_f+1, j_f,   ivc) = f0 + fx - fy + f2 + fxy(2) + &
               box_c%cc(i_f+1, j_f,   ivc)
          box_c%cc(i_f,   j_f+1, ivc) = f0 - fx + fy + f2 + fxy(3) + &
               box_c%cc(i_f,   j_f+1, ivc)
          box_c%cc(i_f+1, j_f+1, ivc) = f0 + fx + fy + f2 + fxy(4) + &
               box_c%cc(i_f+1, j_f+1, ivc)
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

             f0 = box_p%cc(i_c, j_c, k_c, iv)
             fx = 0.125_dp * (box_p%cc(i_c+1, j_c, k_c, iv) - &
                  box_p%cc(i_c-1, j_c, k_c, iv))
             fy = 0.125_dp * (box_p%cc(i_c, j_c+1, k_c, iv) - &
                  box_p%cc(i_c, j_c-1, k_c, iv))
             fz = 0.125_dp * (box_p%cc(i_c, j_c, k_c+1, iv) - &
                  box_p%cc(i_c, j_c, k_c-1, iv))
             fxx = 0.03125_dp * (box_p%cc(i_c-1, j_c, k_c, iv) - &
                  2 * f0 + box_p%cc(i_c+1, j_c, k_c, iv))
             fyy = 0.03125_dp * (box_p%cc(i_c, j_c-1, k_c, iv) - &
                  2 * f0 + box_p%cc(i_c, j_c+1, k_c, iv))
             fzz = 0.03125_dp * (box_p%cc(i_c, j_c, k_c-1, iv) - &
                  2 * f0 + box_p%cc(i_c, j_c, k_c+1, iv))
             f2 = fxx + fyy + fzz

             box_c%cc(i_f,   j_f,   k_f,   ivc) = f0 - fx - fy - fz + f2 + &
                  box_c%cc(i_f,   j_f,   k_f,   ivc)
             box_c%cc(i_f+1, j_f,   k_f,   ivc) = f0 + fx - fy - fz + f2 + &
                  box_c%cc(i_f+1, j_f,   k_f,   ivc)
             box_c%cc(i_f,   j_f+1, k_f,   ivc) = f0 - fx + fy - fz + f2 + &
                  box_c%cc(i_f,   j_f+1, k_f,   ivc)
             box_c%cc(i_f+1, j_f+1, k_f,   ivc) = f0 + fx + fy - fz + f2 + &
                  box_c%cc(i_f+1, j_f+1, k_f,   ivc)
             box_c%cc(i_f,   j_f,   k_f+1, ivc) = f0 - fx - fy + fz + f2 + &
                  box_c%cc(i_f,   j_f,   k_f+1, ivc)
             box_c%cc(i_f+1, j_f,   k_f+1, ivc) = f0 + fx - fy + fz + f2 + &
                  box_c%cc(i_f+1, j_f,   k_f+1, ivc)
             box_c%cc(i_f,   j_f+1, k_f+1, ivc) = f0 - fx + fy + fz + f2 + &
                  box_c%cc(i_f,   j_f+1, k_f+1, ivc)
             box_c%cc(i_f+1, j_f+1, k_f+1, ivc) = f0 + fx + fy + fz + f2 + &
                  box_c%cc(i_f+1, j_f+1, k_f+1, ivc)
          end do
       end do
    end do
#endif
  end subroutine a$D_prolong_quadratic

end module m_a$D_prolong
