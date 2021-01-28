! Multigrid operators for internal boundary conditions. A level set function
! defines the location of the interface(s).
module m_mg_lsf
  implicit none
  private

  ! Methods for normal Laplacian with internal boundary conditions (LSF)
  public :: mg_box_lpllsf
  public :: mg_box_gsrb_lpllsf
  public :: mg_box_corr_lpllsf
  public :: mg_box_rstr_lpllsf

contains

  !> For a point a, compute value and distance (between 0, 1) of a neighbor b.
  subroutine lsf_dist_val(lsf_val_bval_a, lsf_val_bval_b, dist, val)
    !> Level set function at a, value at a, boundary value at a
    real(dp), intent(in)  :: lsf_val_bval_a(3)
    !> Level set function at b, value at b, boundary value at b
    real(dp), intent(in)  :: lsf_val_bval_b(3)
    !> Distance to neighbor point (value between 0 and 1)
    real(dp), intent(out) :: dist
    !> Value at neighbor point
    real(dp), intent(out) :: val
    real(dp)              :: lsf_a, lsf_b, bval_a, bval_b

    lsf_a = lsf_val_bval_a(1)
    lsf_b = lsf_val_bval_b(1)

    if (lsf_a * lsf_b < 0) then
       ! There is a boundary between the points
       dist = lsf_a / (lsf_a - lsf_b)
       bval_a = lsf_val_bval_a(3)
       bval_b = lsf_val_bval_b(3)

       ! Interpolate between boundary values
       val  = bval_a * (1-dist) + bval_b * dist
    else
       ! Simply use the value at b
       dist = 1
       val  = lsf_val_bval_b(2)
    end if
  end subroutine lsf_dist_val

  subroutine mg_box_corr_lpllsf(box_p, box_c, mg)
    type(box_t), intent(inout) :: box_c !< Child box
    type(box_t), intent(in)    :: box_p !< Parent box
    type(mg_t), intent(in)     :: mg !< Multigrid options
    integer                     :: i_phi, i_corr, i_lsf, ix_offset(NDIM)
    integer                     :: nc, i, j, i_c1, i_c2, j_c1, j_c2
    real(dp)                    :: v_a(3), v_b(3), val(NDIM+1), dist(NDIM+1), c(NDIM+1)
#if NDIM == 3
    integer                     :: k, k_c1, k_c2
#endif

    nc        = box_c%n_cell
    ix_offset = af_get_child_offset(box_c)
    i_phi     = mg%i_phi
    i_corr    = mg%i_tmp
    i_lsf     = mg%i_lsf

    ! In these loops, we calculate the closest coarse index (_c1), and the
    ! one-but-closest (_c2). The fine cell lies in between.
#if NDIM == 2
    do j = 1, nc
       j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

          v_a(1:2) = box_c%cc(i, j, [i_lsf, i_corr])
          v_a(3) = 0.0_dp       ! Boundary value for correction is 0
          v_b(3) = 0.0_dp       ! Idem
          v_b(1:2) = box_p%cc(i_c1, j_c1, [i_lsf, i_corr])
          call lsf_dist_val(v_a, v_b, dist(1), val(1))
          v_b(1:2) = box_p%cc(i_c2, j_c1, [i_lsf, i_corr])
          call lsf_dist_val(v_a, v_b, dist(2), val(2))
          v_b(1:2) = box_p%cc(i_c1, j_c2, [i_lsf, i_corr])
          call lsf_dist_val(v_a, v_b, dist(3), val(3))

          ! This expresses general interpolation between 3 points (on the lines
          ! between the fine and the 3 coarse values).
          c(1) = 2 * dist(2) * dist(3)
          c(2) = dist(1) * dist(3)
          c(3) = dist(1) * dist(2)
          box_c%cc(i, j, i_phi) = box_c%cc(i, j, i_phi) + sum(c * val)/sum(c)
       end do
    end do
#elif NDIM == 3
    do k = 1, nc
       k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
       k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

             v_a(1:2) = box_c%cc(i, j, k, [i_lsf, i_corr])
             v_a(3) = 0.0_dp       ! Boundary value for correction is 0
             v_b(3) = 0.0_dp       ! Idem
             v_b(1:2) = box_p%cc(i_c1, j_c1, k_c1, [i_lsf, i_corr])
             call lsf_dist_val(v_a, v_b, dist(1), val(1))
             v_b(1:2) = box_p%cc(i_c2, j_c1, k_c1, [i_lsf, i_corr])
             call lsf_dist_val(v_a, v_b, dist(2), val(2))
             v_b(1:2) = box_p%cc(i_c1, j_c2, k_c1, [i_lsf, i_corr])
             call lsf_dist_val(v_a, v_b, dist(3), val(3))
             v_b(1:2) = box_p%cc(i_c1, j_c1, k_c2, [i_lsf, i_corr])
             call lsf_dist_val(v_a, v_b, dist(4), val(4))

             ! This expresses general interpolation between 4 points (on the lines
             ! between the fine and the 4 coarse values).
             c(1) = dist(2) * dist(3) * dist(4)
             c(2) = dist(1) * dist(3) * dist(4)
             c(3) = dist(1) * dist(2) * dist(4)
             c(4) = dist(1) * dist(2) * dist(3)
             box_c%cc(i, j, k, i_phi) = box_c%cc(i, j, k, i_phi) + &
                  sum(c * val)/sum(c)
          end do
       end do
    end do
#endif
  end subroutine mg_box_corr_lpllsf

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator
  subroutine mg_box_gsrb_lpllsf(box, redblack_cntr, mg)
    type(box_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    type(mg_t), intent(in)     :: mg !< Multigrid options
    integer                     :: i, i0, j, nc, i_phi, i_rhs, i_lsf, i_bval
    real(dp)                    :: dx2, dd(2*NDIM), val(2*NDIM), v_a(3), v_b(3)
#if NDIM == 3
    integer                     :: k
#endif

    dx2    = box%dr**2
    nc     = box%n_cell
    i_phi  = mg%i_phi
    i_rhs  = mg%i_rhs
    i_lsf  = mg%i_lsf
    i_bval = mg%i_bval

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
#if NDIM == 2
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          v_a = box%cc(i, j, [i_lsf, i_phi, i_bval])
          v_b = box%cc(i-1, j, [i_lsf, i_phi, i_bval])
          call lsf_dist_val(v_a, v_b, dd(1), val(1))
          v_b = box%cc(i+1, j, [i_lsf, i_phi, i_bval])
          call lsf_dist_val(v_a, v_b, dd(2), val(2))
          v_b = box%cc(i, j-1, [i_lsf, i_phi, i_bval])
          call lsf_dist_val(v_a, v_b, dd(3), val(3))
          v_b = box%cc(i, j+1, [i_lsf, i_phi, i_bval])
          call lsf_dist_val(v_a, v_b, dd(4), val(4))

          ! Solve for generalized Laplacian (see routine mg_box_lpllsf)
          box%cc(i, j, i_phi) = 1 / &
               (dd(1) * dd(2) + dd(3) * dd(4)) * ( &
               (dd(2) * val(1) + dd(1) * val(2)) * &
               dd(3) * dd(4) / (dd(1) + dd(2)) + &
               (dd(4) * val(3) + dd(3) * val(4)) * &
               dd(1) * dd(2) / (dd(3) + dd(4)) - &
               0.5_dp * product(dd) * dx2 * box%cc(i, j, i_rhs))
       end do
    end do
#elif NDIM == 3
    do k = 1, nc
       do j = 1, nc
          i0 = 2 - iand(ieor(redblack_cntr, k+j), 1)
          do i = i0, nc, 2
             v_a = box%cc(i, j, k, [i_lsf, i_phi, i_bval])
             v_b = box%cc(i-1, j, k, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(1), val(1))
             v_b = box%cc(i+1, j, k, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(2), val(2))
             v_b = box%cc(i, j-1, k, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(3), val(3))
             v_b = box%cc(i, j+1, k, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(4), val(4))
             v_b = box%cc(i, j, k-1, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(5), val(5))
             v_b = box%cc(i, j, k+1, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(6), val(6))

             ! Solve for generalized Laplacian (see routine mg_box_lpllsf)
             box%cc(i, j, k, i_phi) = 1 / (1/(dd(1)*dd(2)) + &
                  1/(dd(3)*dd(4)) + 1/(dd(5)*dd(6))) * ( &
                  (dd(2) * val(1) + dd(1) * val(2)) / &
                  ((dd(1) + dd(2)) * dd(1) * dd(2)) + &
                  (dd(4) * val(3) + dd(3) * val(4)) / &
                  ((dd(3) + dd(4)) * dd(3) * dd(4)) + &
                  (dd(6) * val(5) + dd(5) * val(6)) / &
                  ((dd(5) + dd(6)) * dd(5) * dd(6)) - &
                  0.5_dp * dx2 * box%cc(i, j, k, i_rhs))

          end do
       end do
    end do
#endif
  end subroutine mg_box_gsrb_lpllsf

  !> Perform Laplacian operator on a box
  subroutine mg_box_lpllsf(box, i_out, mg)
    type(box_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    type(mg_t), intent(in)     :: mg !< Multigrid options
    integer                     :: i, j, nc, i_phi, i_lsf, i_bval
    real(dp)                    :: idr2, dd(2*NDIM), val(2*NDIM)
    real(dp)                    :: f0, v_a(3), v_b(3)
#if NDIM == 3
    integer                     :: k
#endif

    nc        = box%n_cell
    idr2 = 1 / box%dr**2
    i_phi     = mg%i_phi
    i_lsf     = mg%i_lsf
    i_bval    = mg%i_bval

#if NDIM == 2
    do j = 1, nc
       do i = 1, nc
          v_a = box%cc(i, j, [i_lsf, i_phi, i_bval])
          v_b = box%cc(i-1, j, [i_lsf, i_phi, i_bval])
          call lsf_dist_val(v_a, v_b, dd(1), val(1))
          v_b = box%cc(i+1, j, [i_lsf, i_phi, i_bval])
          call lsf_dist_val(v_a, v_b, dd(2), val(2))
          v_b = box%cc(i, j-1, [i_lsf, i_phi, i_bval])
          call lsf_dist_val(v_a, v_b, dd(3), val(3))
          v_b = box%cc(i, j+1, [i_lsf, i_phi, i_bval])
          call lsf_dist_val(v_a, v_b, dd(4), val(4))

          ! Generalized Laplacian for neighbors at distance dd * dx
          f0 = box%cc(i, j, i_phi)
          box%cc(i, j, i_out) = 2 * idr2 * ( &
               (dd(2) * val(1) + dd(1) * val(2) - (dd(1)+dd(2)) * f0) / &
               ((dd(1) + dd(2)) * dd(1) * dd(2)) + &
               (dd(4) * val(3) + dd(3) * val(4) - (dd(3)+dd(4)) * f0) / &
               ((dd(3) + dd(4)) * dd(3) * dd(4)))
       end do
    end do
#elif NDIM == 3
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             v_a = box%cc(i, j, k, [i_lsf, i_phi, i_bval])
             v_b = box%cc(i-1, j, k, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(1), val(1))
             v_b = box%cc(i+1, j, k, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(2), val(2))
             v_b = box%cc(i, j-1, k, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(3), val(3))
             v_b = box%cc(i, j+1, k, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(4), val(4))
             v_b = box%cc(i, j, k-1, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(5), val(5))
             v_b = box%cc(i, j, k+1, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(6), val(6))

             ! Generalized Laplacian for neighbors at distance dd * dx
             f0 = box%cc(i, j, k, i_phi)
             box%cc(i, j, k, i_out) = 2 * idr2 * ( &
                  (dd(2) * val(1) + dd(1) * val(2) - (dd(1)+dd(2)) * f0) / &
                  ((dd(1) + dd(2)) * dd(1) * dd(2)) + &
                  (dd(4) * val(3) + dd(3) * val(4) - (dd(3)+dd(4)) * f0) / &
                  ((dd(3) + dd(4)) * dd(3) * dd(4)) + &
                  (dd(6) * val(5) + dd(5) * val(6) - (dd(5)+dd(6)) * f0) / &
                  ((dd(5) + dd(6)) * dd(5) * dd(6)))
          end do
       end do
    end do
#endif
  end subroutine mg_box_lpllsf

  !> Restriction of child box (box_c) to its parent (box_p)
  subroutine mg_box_rstr_lpllsf(box_c, box_p, iv, mg)
    type(box_t), intent(in)      :: box_c         !< Child box to restrict
    type(box_t), intent(inout)   :: box_p         !< Parent box to restrict to
    integer, intent(in)           :: iv            !< Variable to restrict
    type(mg_t), intent(in)       :: mg !< Multigrid options
    integer                       :: i, j, i_f, j_f, i_c, j_c
    integer                       :: hnc, ix_offset(NDIM), n_ch
#if NDIM == 2
    logical                       :: child_mask(2, 2)
#elif NDIM == 3
    logical                       :: child_mask(2, 2, 2)
    integer                       :: k, k_f, k_c
#endif

    hnc       = ishft(box_c%n_cell, -1) ! n_cell / 2
    ix_offset = af_get_child_offset(box_c)

#if NDIM == 2
    do j = 1, hnc
       j_c = ix_offset(2) + j
       j_f = 2 * j - 1
       do i = 1, hnc
          i_c = ix_offset(1) + i
          i_f = 2 * i - 1

          child_mask = (box_p%cc(i_c, j_c, mg%i_lsf) * &
               box_c%cc(i_f:i_f+1, j_f:j_f+1, mg%i_lsf) > 0)
          n_ch = count(child_mask)

          if (n_ch < af_num_children .and. n_ch > 0) then
             box_p%cc(i_c, j_c, iv) = 1 / n_ch * &
                  sum(box_c%cc(i_f:i_f+1, j_f:j_f+1, iv), mask=child_mask)
          else                  ! Take average of children
             box_p%cc(i_c, j_c, iv) = 0.25_dp * &
                  sum(box_c%cc(i_f:i_f+1, j_f:j_f+1, iv))
          end if
       end do
    end do
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

             child_mask = (box_p%cc(i_c, j_c, k_c, mg%i_lsf) * &
                  box_c%cc(i_f:i_f+1, j_f:j_f+1, k_f:k_f+1, mg%i_lsf) > 0)
             n_ch = count(child_mask)

             if (n_ch < af_num_children .and. n_ch > 0) then
                box_p%cc(i_c, j_c, k_c, iv) = 1 / n_ch * &
                     sum(box_c%cc(i_f:i_f+1, j_f:j_f+1, k_f:k_f+1, iv), &
                     mask=child_mask)
             else                  ! Take average of children
                box_p%cc(i_c, j_c, k_c, iv) = 0.125_dp * &
                     sum(box_c%cc(i_f:i_f+1, j_f:j_f+1, k_f:k_f+1, iv))
             end if
          end do
       end do
    end do
#endif
  end subroutine mg_box_rstr_lpllsf

end module m_mg_lsf
