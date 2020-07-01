!> This module contains the routines related to prolongation: going from
!> coarse to fine variables.
module m_af_prolong
  use m_af_types

  implicit none
  private

  public :: af_prolong_copy_from
  public :: af_prolong_copy
  public :: af_prolong_zeroth
  public :: af_prolong_linear_from
  public :: af_prolong_sparse
  public :: af_prolong_linear
  ! public :: af_prolong_quadratic_from
  ! public :: af_prolong_quadratic

  public :: af_prolong_limit_pos
  public :: af_prolong_limit
  public :: af_prolong_linear_cons

contains

  !> Zeroth-order prolongation to children.
  subroutine af_prolong_copy_from(boxes, id, iv, iv_to, add)
    type(box_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Box whose children we will fill
    integer, intent(in)           :: iv       !< Variable that is prolonged
    integer, intent(in), optional :: iv_to    !< Destination variable
    logical, intent(in), optional :: add      !< Add to old values
    integer                       :: i_c, c_id

    do i_c = 1, af_num_children
       c_id = boxes(id)%children(i_c)
       if (c_id == af_no_box) cycle
       call af_prolong_copy(boxes(id), boxes(c_id), iv, iv_to=iv_to, add=add)
    end do
  end subroutine af_prolong_copy_from

  !> Partial prolongation to a child (from parent) using injection (simply copy value)
  subroutine af_prolong_copy(box_p, box_c, iv, iv_to, low, high, add)
    type(box_t), intent(in)      :: box_p !< Parent box
    type(box_t), intent(inout)   :: box_c !< Child box
    integer, intent(in)           :: iv       !< Variable to fill
    integer, intent(in), optional :: iv_to    !< Destination variable
    integer, intent(in), optional :: low(NDIM) !< Min cell index at child
    integer, intent(in), optional :: high(NDIM) !< Max cell index at child
    logical, intent(in), optional :: add      !< Add to old values
    logical                       :: add_to
    integer                       :: nc, ix_offset(NDIM), ivc
    integer                       :: i, j, i_c1, j_c1, lo(NDIM), hi(NDIM)
#if NDIM == 3
    integer                       :: k, k_c1
#endif

    nc   = box_c%n_cell
    add_to = .false.; if (present(add)) add_to = add
    ivc = iv; if (present(iv_to)) ivc = iv_to
    lo   = 1; if (present(low)) lo = low
    hi   = nc; if (present(high)) hi = high

    ! Offset of child w.r.t. parent
    ix_offset = af_get_child_offset(box_c)

    if (add_to) then
#if NDIM == 2
       do j = lo(2), hi(2)
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          do i = lo(1), hi(1)
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             box_c%cc(i, j, ivc) = box_c%cc(i, j, ivc) + &
                  box_p%cc(i_c1, j_c1, iv)
          end do
       end do
#elif NDIM == 3
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
#if NDIM == 2
       do j = lo(2), hi(2)
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          do i = lo(1), hi(1)
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             box_c%cc(i, j, ivc) = box_p%cc(i_c1, j_c1, iv)
          end do
       end do
#elif NDIM == 3
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
  end subroutine af_prolong_copy

  !> Zeroth order prolongation
  subroutine af_prolong_zeroth(box_p, box_c, iv, iv_to, add)
    type(box_t), intent(in)       :: box_p !< Parent box
    type(box_t), intent(inout)    :: box_c !< Child box
    integer, intent(in)           :: iv    !< Variable to fill
    integer, intent(in), optional :: iv_to !< Destination variable
    logical, intent(in), optional :: add   !< Add to old values
    integer                       :: hnc, nc, ix_offset(NDIM), ivc
    integer                       :: i, j, i_c, i_f, j_c, j_f
    real(dp)                      :: f0
    logical                       :: add_to
#if NDIM == 3
    integer                       :: k, k_c, k_f
#endif

    nc        = box_c%n_cell
    hnc       = ishft(box_c%n_cell, -1)
    ix_offset = af_get_child_offset(box_c)
    add_to    = .false.; if (present(add)) add_to = add
    ivc       = iv; if (present(iv_to)) ivc = iv_to

    if (.not. add_to) then
#if NDIM == 2
       box_c%cc(1:nc, 1:nc, ivc) = 0
#elif NDIM == 3
       box_c%cc(1:nc, 1:nc, 1:nc, ivc) = 0
#endif
    end if

#if NDIM == 2
    do j = 1, hnc
       j_c = j + ix_offset(2)
       j_f = 2 * j - 1
       do i = 1, hnc
          i_c = i + ix_offset(1)
          i_f = 2 * i - 1

          f0 = box_p%cc(i_c, j_c, iv)
          box_c%cc(i_f:i_f+1, j_f:j_f+1, ivc) = &
               box_c%cc(i_f:i_f+1, j_f:j_f+1, ivc) + f0
       end do
    end do
#elif NDIM == 3
    do k = 1, hnc
       k_c = k + ix_offset(3)
       k_f = 2 * k - 1
       do j = 1, hnc
          j_c = j + ix_offset(2)
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = i + ix_offset(1)
             i_f = 2 * i - 1

             f0  = box_p%cc(i_c,   j_c,   k_c,   iv)
             box_c%cc(i_f:i_f+1, j_f:j_f+1, k_f:k_f+1, ivc) = &
               box_c%cc(i_f:i_f+1, j_f:j_f+1, k_f:k_f+1, ivc) + f0
          end do
       end do
    end do
#endif
  end subroutine af_prolong_zeroth

  !> Linear prolongation to children. We use 2-1-1 interpolation (2d) and
  !> 1-1-1-1 interpolation (3D), which do not require corner ghost cells.
  subroutine af_prolong_linear_from(boxes, id, iv, iv_to, add)
    type(box_t), intent(inout)  :: boxes(:) !< List of all boxes
    integer, intent(in)           :: id       !< Box whose children we will fill
    integer, intent(in)           :: iv       !< Variable that is filled
    integer, intent(in), optional :: iv_to    !< Destination variable
    logical, intent(in), optional :: add      !< Add to old values
    integer                       :: i_c, c_id

    do i_c = 1, af_num_children
       c_id = boxes(id)%children(i_c)
       if (c_id == af_no_box) cycle
       call af_prolong_linear(boxes(id), boxes(c_id), iv, iv_to, add)
    end do
  end subroutine af_prolong_linear_from

  !> Prolongation to a child (from parent) using linear interpolation. We use
  !> 2-1-1 interpolation (2D) and 1-1-1-1 interpolation (3D) which do not need
  !> corner ghost cells.
  subroutine af_prolong_sparse(box_p, box_c, iv, iv_to, add)
    type(box_t), intent(in)      :: box_p !< Parent box
    type(box_t), intent(inout)   :: box_c !< Child box
    integer, intent(in)           :: iv       !< Variable to fill
    integer, intent(in), optional :: iv_to    !< Destination variable
    logical, intent(in), optional :: add      !< Add to old values
    integer                       :: hnc, nc, ix_offset(NDIM), ivc
    integer                       :: i, j, i_c, i_f, j_c, j_f
    real(dp)                      :: f0, flx, fhx, fly, fhy
    logical                       :: add_to
#if NDIM == 3
    real(dp)                      :: flz, fhz
    integer                       :: k, k_c, k_f
#endif

    nc        = box_c%n_cell
    hnc       = ishft(box_c%n_cell, -1)
    ix_offset = af_get_child_offset(box_c)
    add_to    = .false.; if (present(add)) add_to = add
    ivc       = iv; if (present(iv_to)) ivc = iv_to

    if (.not. add_to) then
#if NDIM == 2
       box_c%cc(1:nc, 1:nc, ivc) = 0
#elif NDIM == 3
       box_c%cc(1:nc, 1:nc, 1:nc, ivc) = 0
#endif
    end if

#if NDIM == 2
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
#elif NDIM == 3
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
  end subroutine af_prolong_sparse

  !> Conservative prolongation using the gradient of the coarse cells, and
  !> limited to preserve positivity
  subroutine af_prolong_limit_pos(box_p, box_c, iv, iv_to, add)
    type(box_t), intent(in)     :: box_p !< Parent box
    type(box_t), intent(inout)  :: box_c !< Child box
    integer, intent(in)           :: iv    !< Variable to fill
    integer, intent(in), optional :: iv_to !< Destination variable
    logical, intent(in), optional :: add   !< Add to old values
    integer                       :: hnc, nc, ix_offset(NDIM), ivc
    integer                       :: i, j, i_c, i_f, j_c, j_f
    real(dp)                      :: f0, fx, fy, tmp
    logical                       :: add_to
#if NDIM == 3
    real(dp)                      :: fz
    integer                       :: k, k_c, k_f
#endif

    nc        = box_c%n_cell
    hnc       = ishft(box_c%n_cell, -1)
    ix_offset = af_get_child_offset(box_c)
    add_to    = .false.; if (present(add)) add_to = add
    ivc       = iv; if (present(iv_to)) ivc = iv_to

    if (.not. add_to) then
#if NDIM == 2
       box_c%cc(1:nc, 1:nc, ivc) = 0
#elif NDIM == 3
       box_c%cc(1:nc, 1:nc, 1:nc, ivc) = 0
#endif
    end if

#if NDIM == 2
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

          if (box_p%coord_t == af_cyl) then
             ! Ensure densities stay positive
             tmp = abs(0.5_dp * f0 / (1 + 0.25_dp * box_p%dr(1) / &
                  af_cyl_radius_cc(box_p, i_c)))
             if (abs(fx) > tmp) fx = sign(tmp, fx)
             if (abs(fy) > tmp) fy = sign(tmp, fy)

             ! Correction for cylindrical coordinates
             f0 = f0 - 0.25_dp * box_p%dr(1) * fx / &
                  af_cyl_radius_cc(box_p, i_c)
          else
             ! Ensure densities stay positive
             tmp = abs(0.5_dp * f0)
             if (abs(fx) > tmp) fx = sign(tmp, fx)
             if (abs(fy) > tmp) fy = sign(tmp, fy)
          end if

          box_c%cc(i_f,   j_f,   ivc) = f0 - fx - fy &
               + box_c%cc(i_f,   j_f,   ivc)
          box_c%cc(i_f+1, j_f,   ivc) = f0 + fx - fy &
               + box_c%cc(i_f+1, j_f,   ivc)
          box_c%cc(i_f,   j_f+1, ivc) = f0 - fx + fy &
               + box_c%cc(i_f,   j_f+1, ivc)
          box_c%cc(i_f+1, j_f+1, ivc) = f0 + fx + fy &
               + box_c%cc(i_f+1, j_f+1, ivc)
       end do
    end do
#elif NDIM == 3
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
             fx = 0.125_dp * (box_p%cc(i_c+1, j_c,   k_c,   iv) - &
                  box_p%cc(i_c-1, j_c,   k_c,   iv))
             fy = 0.125_dp * (box_p%cc(i_c,   j_c+1, k_c,   iv) - &
                  box_p%cc(i_c,   j_c-1, k_c,   iv))
             fz = 0.125_dp * (box_p%cc(i_c,   j_c,   k_c+1, iv) - &
                  box_p%cc(i_c,   j_c,   k_c-1, iv))

             ! Ensure densities stay positive
             tmp = abs(f0/3.0_dp)
             if (abs(fx) > tmp) fx = sign(tmp, fx)
             if (abs(fy) > tmp) fy = sign(tmp, fy)
             if (abs(fz) > tmp) fz = sign(tmp, fz)

             box_c%cc(i_f,   j_f,   k_f,   ivc) = f0 - fx - &
                  fy - fz + box_c%cc(i_f,   j_f,   k_f,   ivc)
             box_c%cc(i_f+1, j_f,   k_f,   ivc) = f0 + fx - &
                  fy - fz + box_c%cc(i_f+1, j_f,   k_f,   ivc)
             box_c%cc(i_f,   j_f+1, k_f,   ivc) = f0 - fx + &
                  fy - fz + box_c%cc(i_f,   j_f+1, k_f,   ivc)
             box_c%cc(i_f+1, j_f+1, k_f,   ivc) = f0 + fx + &
                  fy - fz + box_c%cc(i_f+1, j_f+1, k_f,   ivc)
             box_c%cc(i_f,   j_f,   k_f+1, ivc) = f0 - fx - &
                  fy + fz + box_c%cc(i_f,   j_f,   k_f+1, ivc)
             box_c%cc(i_f+1, j_f,   k_f+1, ivc) = f0 + fx - &
                  fy + fz + box_c%cc(i_f+1, j_f,   k_f+1, ivc)
             box_c%cc(i_f,   j_f+1, k_f+1, ivc) = f0 - fx + &
                  fy + fz + box_c%cc(i_f,   j_f+1, k_f+1, ivc)
             box_c%cc(i_f+1, j_f+1, k_f+1, ivc) = f0 + fx + &
                  fy + fz + box_c%cc(i_f+1, j_f+1, k_f+1, ivc)
          end do
       end do
    end do
#endif
  end subroutine af_prolong_limit_pos

  !> Conservative prolongation using the gradient from the coarse cells, taking
  !> the minimum of the slopes and zero if they differ
  subroutine af_prolong_limit(box_p, box_c, iv, iv_to, add)
    type(box_t), intent(in)     :: box_p !< Parent box
    type(box_t), intent(inout)  :: box_c !< Child box
    integer, intent(in)           :: iv    !< Variable to fill
    integer, intent(in), optional :: iv_to !< Destination variable
    logical, intent(in), optional :: add   !< Add to old values
    integer                       :: hnc, nc, ix_offset(NDIM), ivc
    integer                       :: i, j, i_c, i_f, j_c, j_f
    real(dp)                      :: f0, fx, fy
    logical                       :: add_to
#if NDIM == 3
    real(dp)                      :: fz
    integer                       :: k, k_c, k_f
#endif

    nc        = box_c%n_cell
    hnc       = ishft(box_c%n_cell, -1)
    ix_offset = af_get_child_offset(box_c)
    add_to    = .false.; if (present(add)) add_to = add
    ivc       = iv; if (present(iv_to)) ivc = iv_to

    if (.not. add_to) then
#if NDIM == 2
       box_c%cc(1:nc, 1:nc, ivc) = 0
#elif NDIM == 3
       box_c%cc(1:nc, 1:nc, 1:nc, ivc) = 0
#endif
    end if

#if NDIM == 2
    do j = 1, hnc
       j_c = j + ix_offset(2)
       j_f = 2 * j - 1
       do i = 1, hnc
          i_c = i + ix_offset(1)
          i_f = 2 * i - 1

          f0 = box_p%cc(i_c, j_c, iv)
          fx = 0.25_dp * limit_slope( &
               box_p%cc(i_c, j_c, iv) - box_p%cc(i_c-1, j_c, iv), &
               box_p%cc(i_c+1, j_c, iv) - box_p%cc(i_c, j_c, iv))
          fy = 0.25_dp * limit_slope( &
               box_p%cc(i_c, j_c, iv) - box_p%cc(i_c, j_c-1, iv), &
               box_p%cc(i_c, j_c+1, iv) - box_p%cc(i_c, j_c, iv))

          box_c%cc(i_f,   j_f,   ivc) = f0 - fx - fy &
               + box_c%cc(i_f,   j_f,   ivc)
          box_c%cc(i_f+1, j_f,   ivc) = f0 + fx - fy &
               + box_c%cc(i_f+1, j_f,   ivc)
          box_c%cc(i_f,   j_f+1, ivc) = f0 - fx + fy &
               + box_c%cc(i_f,   j_f+1, ivc)
          box_c%cc(i_f+1, j_f+1, ivc) = f0 + fx + fy &
               + box_c%cc(i_f+1, j_f+1, ivc)
       end do
    end do
#elif NDIM == 3
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
             fx = 0.25_dp * limit_slope( &
                  box_p%cc(i_c, j_c, k_c, iv) - box_p%cc(i_c-1, j_c, k_c, iv), &
                  box_p%cc(i_c+1, j_c, k_c, iv) - box_p%cc(i_c, j_c, k_c, iv))
             fy = 0.25_dp * limit_slope( &
                  box_p%cc(i_c, j_c, k_c, iv) - box_p%cc(i_c, j_c-1, k_c, iv), &
                  box_p%cc(i_c, j_c+1, k_c, iv) - box_p%cc(i_c, j_c, k_c, iv))
             fz = 0.25_dp * limit_slope( &
                  box_p%cc(i_c, j_c, k_c, iv) - box_p%cc(i_c, j_c, k_c-1, iv), &
                  box_p%cc(i_c, j_c, k_c+1, iv) - box_p%cc(i_c, j_c, k_c, iv))

             box_c%cc(i_f,   j_f,   k_f,   ivc) = f0 - fx - &
                  fy - fz + box_c%cc(i_f,   j_f,   k_f,   ivc)
             box_c%cc(i_f+1, j_f,   k_f,   ivc) = f0 + fx - &
                  fy - fz + box_c%cc(i_f+1, j_f,   k_f,   ivc)
             box_c%cc(i_f,   j_f+1, k_f,   ivc) = f0 - fx + &
                  fy - fz + box_c%cc(i_f,   j_f+1, k_f,   ivc)
             box_c%cc(i_f+1, j_f+1, k_f,   ivc) = f0 + fx + &
                  fy - fz + box_c%cc(i_f+1, j_f+1, k_f,   ivc)
             box_c%cc(i_f,   j_f,   k_f+1, ivc) = f0 - fx - &
                  fy + fz + box_c%cc(i_f,   j_f,   k_f+1, ivc)
             box_c%cc(i_f+1, j_f,   k_f+1, ivc) = f0 + fx - &
                  fy + fz + box_c%cc(i_f+1, j_f,   k_f+1, ivc)
             box_c%cc(i_f,   j_f+1, k_f+1, ivc) = f0 - fx + &
                  fy + fz + box_c%cc(i_f,   j_f+1, k_f+1, ivc)
             box_c%cc(i_f+1, j_f+1, k_f+1, ivc) = f0 + fx + &
                  fy + fz + box_c%cc(i_f+1, j_f+1, k_f+1, ivc)
          end do
       end do
    end do
#endif
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

  end subroutine af_prolong_limit

  ! Prolong with a linear (unlimited) slope in the coarse cells, which can
  ! result in negative densities. This procedure is conservative.
  subroutine af_prolong_linear_cons(box_p, box_c, iv, iv_to, add)
    type(box_t), intent(in)     :: box_p !< Parent box
    type(box_t), intent(inout)  :: box_c !< Child box
    integer, intent(in)           :: iv    !< Variable to fill
    integer, intent(in), optional :: iv_to !< Destination variable
    logical, intent(in), optional :: add   !< Add to old values
    integer                       :: hnc, nc, ix_offset(NDIM), ivc
    integer                       :: i, j, i_c, i_f, j_c, j_f
    real(dp)                      :: f0, fx, fy
    logical                       :: add_to
#if NDIM == 3
    real(dp)                      :: fz
    integer                       :: k, k_c, k_f
#endif

    nc        = box_c%n_cell
    hnc       = ishft(box_c%n_cell, -1)
    ix_offset = af_get_child_offset(box_c)
    add_to    = .false.; if (present(add)) add_to = add
    ivc       = iv; if (present(iv_to)) ivc = iv_to

    if (.not. add_to) then
#if NDIM == 2
       box_c%cc(1:nc, 1:nc, ivc) = 0
#elif NDIM == 3
       box_c%cc(1:nc, 1:nc, 1:nc, ivc) = 0
#endif
    end if

#if NDIM == 2
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

          if (box_p%coord_t == af_cyl) then
             ! Conservative prolongation for cylindrical coords
             f0 = f0 - 0.25_dp * box_p%dr(1) * fx / &
                  af_cyl_radius_cc(box_p, i_c)
          end if

          box_c%cc(i_f,   j_f,   ivc) = f0 - fx - fy &
               + box_c%cc(i_f,   j_f,   ivc)
          box_c%cc(i_f+1, j_f,   ivc) = f0 + fx - fy &
               + box_c%cc(i_f+1, j_f,   ivc)
          box_c%cc(i_f,   j_f+1, ivc) = f0 - fx + fy &
               + box_c%cc(i_f,   j_f+1, ivc)
          box_c%cc(i_f+1, j_f+1, ivc) = f0 + fx + fy &
               + box_c%cc(i_f+1, j_f+1, ivc)
       end do
    end do
#elif NDIM == 3
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
             fx = 0.125_dp * (box_p%cc(i_c+1, j_c,   k_c,   iv) - &
                  box_p%cc(i_c-1, j_c,   k_c,   iv))
             fy = 0.125_dp * (box_p%cc(i_c,   j_c+1, k_c,   iv) - &
                  box_p%cc(i_c,   j_c-1, k_c,   iv))
             fz = 0.125_dp * (box_p%cc(i_c,   j_c,   k_c+1, iv) - &
                  box_p%cc(i_c,   j_c,   k_c-1, iv))

             box_c%cc(i_f,   j_f,   k_f,   ivc) = f0 - fx - &
                  fy - fz + box_c%cc(i_f,   j_f,   k_f,   ivc)
             box_c%cc(i_f+1, j_f,   k_f,   ivc) = f0 + fx - &
                  fy - fz + box_c%cc(i_f+1, j_f,   k_f,   ivc)
             box_c%cc(i_f,   j_f+1, k_f,   ivc) = f0 - fx + &
                  fy - fz + box_c%cc(i_f,   j_f+1, k_f,   ivc)
             box_c%cc(i_f+1, j_f+1, k_f,   ivc) = f0 + fx + &
                  fy - fz + box_c%cc(i_f+1, j_f+1, k_f,   ivc)
             box_c%cc(i_f,   j_f,   k_f+1, ivc) = f0 - fx - &
                  fy + fz + box_c%cc(i_f,   j_f,   k_f+1, ivc)
             box_c%cc(i_f+1, j_f,   k_f+1, ivc) = f0 + fx - &
                  fy + fz + box_c%cc(i_f+1, j_f,   k_f+1, ivc)
             box_c%cc(i_f,   j_f+1, k_f+1, ivc) = f0 - fx + &
                  fy + fz + box_c%cc(i_f,   j_f+1, k_f+1, ivc)
             box_c%cc(i_f+1, j_f+1, k_f+1, ivc) = f0 + fx + &
                  fy + fz + box_c%cc(i_f+1, j_f+1, k_f+1, ivc)
          end do
       end do
    end do
#endif

  end subroutine af_prolong_linear_cons

  !> Bi/trilinear prolongation to a child (from parent)
  subroutine af_prolong_linear(box_p, box_c, iv, iv_to, add)
    type(box_t), intent(in)     :: box_p !< Parent box
    type(box_t), intent(inout)  :: box_c !< Child box
    integer, intent(in)           :: iv    !< Variable to fill
    integer, intent(in), optional :: iv_to !< Destination variable
    logical, intent(in), optional :: add   !< Add to old values
    integer                       :: hnc, nc, ix_offset(NDIM), ivc
    integer                       :: i, j, i_c, i_f, j_c, j_f
    logical                       :: add_to
#if NDIM == 2
    real(dp)                      :: f0, flx, fhx, fly, fhy
    real(dp)                      :: fll, fhl, flh, fhh
    real(dp), parameter           :: f1  = 1/16.0_dp, f3=3/16.0_dp, f9=9/16.0_dp
#elif NDIM == 3
    real(dp)                      :: f000, f00l, f0l0, f0ll, fl00, fl0l, fll0
    real(dp)                      :: flll, f00h, f0h0, f0hh, fh00, fh0h, fhh0
    real(dp)                      :: fhhh, f0lh, f0hl, fl0h, fh0l, flh0, fhl0
    real(dp)                      :: fllh, flhl, fhll, fhhl, fhlh, flhh
    real(dp), parameter           :: f1  = 1/64.0_dp, f3=3/64.0_dp, f9=9/64.0_dp
    real(dp), parameter           :: f27 = 27/64.0_dp
    integer                       :: k, k_c, k_f
#endif

    nc        = box_c%n_cell
    hnc       = ishft(box_c%n_cell, -1)
    ix_offset = af_get_child_offset(box_c)
    add_to    = .false.; if (present(add)) add_to = add
    ivc       = iv; if (present(iv_to)) ivc = iv_to

    if (.not. add_to) then
#if NDIM == 2
       box_c%cc(1:nc, 1:nc, ivc) = 0
#elif NDIM == 3
       box_c%cc(1:nc, 1:nc, 1:nc, ivc) = 0
#endif
    end if

#if NDIM == 2
    do j = 1, hnc
       j_c = j + ix_offset(2)
       j_f = 2 * j - 1
       do i = 1, hnc
          i_c = i + ix_offset(1)
          i_f = 2 * i - 1

          f0 = f9 * box_p%cc(i_c, j_c, iv)
          flx = f3 * box_p%cc(i_c-1, j_c, iv)
          fhx = f3 * box_p%cc(i_c+1, j_c, iv)
          fly = f3 * box_p%cc(i_c, j_c-1, iv)
          fhy = f3 * box_p%cc(i_c, j_c+1, iv)
          fll = f1 * box_p%cc(i_c-1, j_c-1, iv)
          fhl = f1 * box_p%cc(i_c+1, j_c-1, iv)
          flh = f1 * box_p%cc(i_c-1, j_c+1, iv)
          fhh = f1 * box_p%cc(i_c+1, j_c+1, iv)

          box_c%cc(i_f,   j_f,   ivc) = f0 + flx + fly + fll &
               + box_c%cc(i_f,   j_f,   ivc)
          box_c%cc(i_f+1, j_f,   ivc) = f0 + fhx + fly + fhl &
               + box_c%cc(i_f+1, j_f,   ivc)
          box_c%cc(i_f,   j_f+1, ivc) = f0 + flx + fhy + flh &
               + box_c%cc(i_f,   j_f+1, ivc)
          box_c%cc(i_f+1, j_f+1, ivc) = f0 + fhx + fhy + fhh &
               + box_c%cc(i_f+1, j_f+1, ivc)
       end do
    end do
#elif NDIM == 3
    do k = 1, hnc
       k_c = k + ix_offset(3)
       k_f = 2 * k - 1
       do j = 1, hnc
          j_c = j + ix_offset(2)
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = i + ix_offset(1)
             i_f = 2 * i - 1

             f000 = f27 * box_p%cc(i_c,  j_c,   k_c,   iv)

             f00l = f9 * box_p%cc(i_c,   j_c,   k_c-1, iv)
             f0l0 = f9 * box_p%cc(i_c,   j_c-1, k_c,   iv)
             f0ll = f3 * box_p%cc(i_c,   j_c-1, k_c-1, iv)
             fl00 = f9 * box_p%cc(i_c-1, j_c,   k_c,   iv)
             fl0l = f3 * box_p%cc(i_c-1, j_c,   k_c-1, iv)
             fll0 = f3 * box_p%cc(i_c-1, j_c-1, k_c,   iv)
             flll = f1 * box_p%cc(i_c-1, j_c-1, k_c-1, iv)

             f00h = f9 * box_p%cc(i_c,   j_c,   k_c+1, iv)
             f0h0 = f9 * box_p%cc(i_c,   j_c+1, k_c,   iv)
             f0hh = f3 * box_p%cc(i_c,   j_c+1, k_c+1, iv)
             fh00 = f9 * box_p%cc(i_c+1, j_c,   k_c,   iv)
             fh0h = f3 * box_p%cc(i_c+1, j_c,   k_c+1, iv)
             fhh0 = f3 * box_p%cc(i_c+1, j_c+1, k_c,   iv)
             fhhh = f1 * box_p%cc(i_c+1, j_c+1, k_c+1, iv)

             fl0h = f3 * box_p%cc(i_c-1, j_c,   k_c+1, iv)
             fh0l = f3 * box_p%cc(i_c+1, j_c,   k_c-1, iv)
             flh0 = f3 * box_p%cc(i_c-1, j_c+1, k_c,   iv)
             fhl0 = f3 * box_p%cc(i_c+1, j_c-1, k_c,   iv)
             f0lh = f3 * box_p%cc(i_c,   j_c-1, k_c+1, iv)
             f0hl = f3 * box_p%cc(i_c,   j_c+1, k_c-1, iv)

             fllh = f1 * box_p%cc(i_c-1, j_c-1, k_c+1, iv)
             flhl = f1 * box_p%cc(i_c-1, j_c+1, k_c-1, iv)
             fhll = f1 * box_p%cc(i_c+1, j_c-1, k_c-1, iv)

             fhhl = f1 * box_p%cc(i_c+1, j_c+1, k_c-1, iv)
             fhlh = f1 * box_p%cc(i_c+1, j_c-1, k_c+1, iv)
             flhh = f1 * box_p%cc(i_c-1, j_c+1, k_c+1, iv)

             box_c%cc(i_f,   j_f,   k_f,   ivc) = f000 + fl00 + &
                  f0l0 + f00l + fll0 + fl0l + f0ll + flll + &
                  box_c%cc(i_f,   j_f,   k_f,   ivc)
             box_c%cc(i_f+1, j_f,   k_f,   ivc) = f000 + fh00 + &
                  f0l0 + f00l + fhl0 + fh0l + f0ll + fhll + &
                  box_c%cc(i_f+1, j_f,   k_f,   ivc)
             box_c%cc(i_f,   j_f+1, k_f,   ivc) = f000 + fl00 + &
                  f0h0 + f00l + flh0 + fl0l + f0hl + flhl + &
                  box_c%cc(i_f,   j_f+1, k_f,   ivc)
             box_c%cc(i_f+1, j_f+1, k_f,   ivc) = f000 + fh00 + &
                  f0h0 + f00l + fhh0 + fh0l + f0hl + fhhl + &
                  box_c%cc(i_f+1, j_f+1, k_f,   ivc)
             box_c%cc(i_f,   j_f,   k_f+1, ivc) = f000 + fl00 + &
                  f0l0 + f00h + fll0 + fl0h + f0lh + fllh + &
                  box_c%cc(i_f,   j_f,   k_f+1, ivc)
             box_c%cc(i_f+1, j_f,   k_f+1, ivc) = f000 + fh00 + &
                  f0l0 + f00h + fhl0 + fh0h + f0lh + fhlh + &
                  box_c%cc(i_f+1, j_f,   k_f+1, ivc)
             box_c%cc(i_f,   j_f+1, k_f+1, ivc) = f000 + fl00 + &
                  f0h0 + f00h + flh0 + fl0h + f0hh + flhh + &
                  box_c%cc(i_f,   j_f+1, k_f+1, ivc)
             box_c%cc(i_f+1, j_f+1, k_f+1, ivc) = f000 + fh00 + &
                  f0h0 + f00h + fhh0 + fh0h + f0hh + fhhh + &
                  box_c%cc(i_f+1, j_f+1, k_f+1, ivc)
          end do
       end do
    end do
#endif
  end subroutine af_prolong_linear

  ! !> Quadratic prolongation to children. We use stencils that do not require
  ! !> corner ghost cells.
  ! subroutine af_prolong_quadratic_from(boxes, id, iv, iv_to, add)
  !   type(box_t), intent(inout)  :: boxes(:) !< List of all boxes
  !   integer, intent(in)           :: id       !< Box whose children we will fill
  !   integer, intent(in)           :: iv       !< Variable that is filled
  !   integer, intent(in), optional :: iv_to    !< Destination variable
  !   logical, intent(in), optional :: add      !< Add to old values
  !   integer                       :: i_c, c_id

  !   do i_c = 1, af_num_children
  !      c_id = boxes(id)%children(i_c)
  !      if (c_id == af_no_box) cycle
  !      call af_prolong_quadratic(boxes(id), boxes(c_id), iv, iv_to, add)
  !   end do
  ! end subroutine af_prolong_quadratic_from

!   !> Prolongation to a child (from parent) using quadratic interpolation (using
!   !> a local Taylor approximation)
!   !> \todo Mixed derivatives in 3D version
!   subroutine af_prolong_quadratic(box_p, box_c, iv, iv_to, add)
!     type(box_t), intent(in)      :: box_p !< Parent box
!     type(box_t), intent(inout)   :: box_c !< Child box
!     integer, intent(in)           :: iv       !< Variable to fill
!     integer, intent(in), optional :: iv_to    !< Destination variable
!     logical, intent(in), optional :: add      !< Add to old values
!     logical                      :: add_to
!     integer                      :: hnc, ix_offset(NDIM)
!     integer                      :: i, j, nc, ivc
!     integer                      :: i_c, i_f, j_c, j_f
!     real(dp)                     :: f0, fx, fy, fxx, fyy, f2
! #if NDIM == 2
!     real(dp)                     :: fxy(2**NDIM)
! #elif NDIM == 3
!     real(dp)                     :: fz, fzz
!     integer                      :: k, k_c, k_f
! #endif

!     nc        = box_c%n_cell
!     hnc       = ishft(box_c%n_cell, -1)
!     ix_offset = af_get_child_offset(box_c)
!     add_to    = .false.; if (present(add)) add_to = add
!     ivc       = iv; if (present(iv_to)) ivc = iv_to

!     if (.not. add_to) then
! #if NDIM == 2
!        box_c%cc(1:nc, 1:nc, ivc) = 0
! #elif NDIM == 3
!        box_c%cc(1:nc, 1:nc, 1:nc, ivc) = 0
! #endif
!     end if

! #if NDIM == 2
!     do j = 1, hnc
!        j_c = j + ix_offset(2)
!        j_f = 2 * j - 1
!        do i = 1, hnc
!           i_c = i + ix_offset(1)
!           i_f = 2 * i - 1

!           f0 = box_p%cc(i_c, j_c, iv)
!           fx = 0.125_dp * (box_p%cc(i_c+1, j_c, iv) - &
!                box_p%cc(i_c-1, j_c, iv))
!           fy = 0.125_dp * (box_p%cc(i_c, j_c+1, iv) - &
!                box_p%cc(i_c, j_c-1, iv))
!           fxx = 0.03125_dp * (box_p%cc(i_c-1, j_c, iv) - &
!                2 * f0 + box_p%cc(i_c+1, j_c, iv))
!           fyy = 0.03125_dp * (box_p%cc(i_c, j_c-1, iv) - &
!                2 * f0 + box_p%cc(i_c, j_c+1, iv))
!           f2 = fxx + fyy

!           fxy(1) = 0.0625_dp * (box_p%cc(i_c-1, j_c-1, iv) + f0 - &
!                box_p%cc(i_c-1, j_c, iv) - box_p%cc(i_c, j_c-1, iv))
!           fxy(2) = 0.0625_dp * (box_p%cc(i_c+1, j_c-1, iv) + f0 - &
!                box_p%cc(i_c+1, j_c, iv) - box_p%cc(i_c, j_c-1, iv))
!           fxy(3) = 0.0625_dp * (box_p%cc(i_c-1, j_c+1, iv) + f0 - &
!                box_p%cc(i_c-1, j_c, iv) - box_p%cc(i_c, j_c+1, iv))
!           fxy(4) = 0.0625_dp * (box_p%cc(i_c+1, j_c+1, iv) + f0 - &
!                box_p%cc(i_c+1, j_c, iv) - box_p%cc(i_c, j_c+1, iv))

!           box_c%cc(i_f,   j_f,   ivc) = f0 - fx - fy + f2 + fxy(1) + &
!                box_c%cc(i_f,   j_f,   ivc)
!           box_c%cc(i_f+1, j_f,   ivc) = f0 + fx - fy + f2 + fxy(2) + &
!                box_c%cc(i_f+1, j_f,   ivc)
!           box_c%cc(i_f,   j_f+1, ivc) = f0 - fx + fy + f2 + fxy(3) + &
!                box_c%cc(i_f,   j_f+1, ivc)
!           box_c%cc(i_f+1, j_f+1, ivc) = f0 + fx + fy + f2 + fxy(4) + &
!                box_c%cc(i_f+1, j_f+1, ivc)
!        end do
!     end do
! #elif NDIM == 3
!     do k = 1, hnc
!        k_c = k + ix_offset(3)
!        k_f = 2 * k - 1
!        do j = 1, hnc
!           j_c = j + ix_offset(2)
!           j_f = 2 * j - 1
!           do i = 1, hnc
!              i_c = i + ix_offset(1)
!              i_f = 2 * i - 1

!              f0 = box_p%cc(i_c, j_c, k_c, iv)
!              fx = 0.125_dp * (box_p%cc(i_c+1, j_c, k_c, iv) - &
!                   box_p%cc(i_c-1, j_c, k_c, iv))
!              fy = 0.125_dp * (box_p%cc(i_c, j_c+1, k_c, iv) - &
!                   box_p%cc(i_c, j_c-1, k_c, iv))
!              fz = 0.125_dp * (box_p%cc(i_c, j_c, k_c+1, iv) - &
!                   box_p%cc(i_c, j_c, k_c-1, iv))
!              fxx = 0.03125_dp * (box_p%cc(i_c-1, j_c, k_c, iv) - &
!                   2 * f0 + box_p%cc(i_c+1, j_c, k_c, iv))
!              fyy = 0.03125_dp * (box_p%cc(i_c, j_c-1, k_c, iv) - &
!                   2 * f0 + box_p%cc(i_c, j_c+1, k_c, iv))
!              fzz = 0.03125_dp * (box_p%cc(i_c, j_c, k_c-1, iv) - &
!                   2 * f0 + box_p%cc(i_c, j_c, k_c+1, iv))
!              f2 = fxx + fyy + fzz

!              box_c%cc(i_f,   j_f,   k_f,   ivc) = f0 - fx - fy - fz + f2 + &
!                   box_c%cc(i_f,   j_f,   k_f,   ivc)
!              box_c%cc(i_f+1, j_f,   k_f,   ivc) = f0 + fx - fy - fz + f2 + &
!                   box_c%cc(i_f+1, j_f,   k_f,   ivc)
!              box_c%cc(i_f,   j_f+1, k_f,   ivc) = f0 - fx + fy - fz + f2 + &
!                   box_c%cc(i_f,   j_f+1, k_f,   ivc)
!              box_c%cc(i_f+1, j_f+1, k_f,   ivc) = f0 + fx + fy - fz + f2 + &
!                   box_c%cc(i_f+1, j_f+1, k_f,   ivc)
!              box_c%cc(i_f,   j_f,   k_f+1, ivc) = f0 - fx - fy + fz + f2 + &
!                   box_c%cc(i_f,   j_f,   k_f+1, ivc)
!              box_c%cc(i_f+1, j_f,   k_f+1, ivc) = f0 + fx - fy + fz + f2 + &
!                   box_c%cc(i_f+1, j_f,   k_f+1, ivc)
!              box_c%cc(i_f,   j_f+1, k_f+1, ivc) = f0 - fx + fy + fz + f2 + &
!                   box_c%cc(i_f,   j_f+1, k_f+1, ivc)
!              box_c%cc(i_f+1, j_f+1, k_f+1, ivc) = f0 + fx + fy + fz + f2 + &
!                   box_c%cc(i_f+1, j_f+1, k_f+1, ivc)
!           end do
!        end do
!     end do
! #endif
!   end subroutine af_prolong_quadratic

end module m_af_prolong
