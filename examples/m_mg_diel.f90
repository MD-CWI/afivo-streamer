!> This module contains operators for multigrid if there is a face-aligned jump
!> in the (dielectric) constant. The jump should be aligned to the coarsest grid
!> that is used.
module m_mg_diel
  use m_afivo_2d
  use m_mg_2d

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

contains

  subroutine gsrb_lpl_box_diel(box, redblack_cntr, mg)
    type(box2_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, i0, j, nc, i_phi, i_eps, i_rhs
    real(dp)                    :: dx2, u0, u(4), a0, a(4), c(4)

    dx2   = box%dr**2
    nc    = box%n_cell
    i_phi = mg%i_phi
    i_eps = mg%i_eps
    i_rhs = mg%i_rhs

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          u0 = box%cc(i, j, i_phi) ! value of phi at i,j
          a0 = box%cc(i, j, i_eps) ! value of eps at i,j
          u(1:2) = box%cc(i-1:i+2:2, j, i_phi) ! values at neighbors
          a(1:2) = box%cc(i-1:i+2:2, j, i_eps)
          u(3:4) = box%cc(i, j-1:j+2:2, i_phi)
          a(3:4) = box%cc(i, j-1:j+2:2, i_eps)
          c(:) = 2 * a0 * a(:) / (a0 + a(:))

          box%cc(i, j, i_phi) = &
               (sum(c(:) * u(:)) - dx2 * box%cc(i, j, i_rhs)) / sum(c)
       end do
    end do
  end subroutine gsrb_lpl_box_diel

  !> Perform Laplacian operator on a box
  subroutine lpl_box_diel(box, i_out, mg)
    type(box2_t), intent(inout) :: box   !< Box to operate on
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, j, nc, i_phi, i_eps
    real(dp)                    :: inv_dr_sq, a0, u0, u(4), a(4)

    nc        = box%n_cell
    inv_dr_sq = 1 / box%dr**2
    i_phi     = mg%i_phi
    i_eps     = mg%i_eps

    do j = 1, nc
       do i = 1, nc
          u0 = box%cc(i, j, i_phi)
          a0 = box%cc(i, j, i_eps)
          u(1:2) = box%cc(i-1:i+2:2, j, i_phi)
          u(3:4) = box%cc(i, j-1:j+2:2, i_phi)
          a(1:2) = box%cc(i-1:i+2:2, j, i_eps)
          a(3:4) = box%cc(i, j-1:j+2:2, i_eps)

          box%cc(i, j, i_out) = inv_dr_sq * 2 * &
               sum(a0*a(:)/(a0 + a(:)) * (u(:) - u0))
       end do
    end do
  end subroutine lpl_box_diel

  subroutine corr_lpl_box_diel(box_p, box_c, mg)
    type(box2_t), intent(inout) :: box_c
    type(box2_t), intent(in)    :: box_p
    type(mg2_t), intent(in)     :: mg
    integer                     :: ix_offset(2), i_phi, i_corr, i_eps
    integer                     :: nc, i, j, i_c1, i_c2, j_c1, j_c2
    real(dp) :: u(3), a(3)

    nc = box_c%n_cell
    ix_offset = a2_get_child_offset(box_c)
    i_phi = mg%i_phi
    i_corr = mg%i_tmp
    i_eps = mg%i_eps

    ! In these loops, we calculate the closest coarse index (_c1), and the
    ! one-but-closest (_c2). The fine cell lies in between.
    do j = 1, nc
       j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

          u(1) = box_p%cc(i_c1, j_c1, i_corr)
          u(2) = box_p%cc(i_c2, j_c1, i_corr)
          u(3) = box_p%cc(i_c1, j_c2, i_corr)
          a(1) = box_p%cc(i_c1, j_c1, i_eps)
          a(2) = box_p%cc(i_c2, j_c1, i_eps)
          a(3) = box_p%cc(i_c1, j_c2, i_eps)

          box_c%cc(i, j, i_phi) = box_c%cc(i, j, i_phi) + 0.5_dp * ( &
               (a(1)*u(1)+a(2)*u(2)) / (a(1) + a(2)) + &
               (a(1)*u(1)+a(3)*u(3)) / (a(1) + a(3)))
       end do
    end do
  end subroutine corr_lpl_box_diel

end module m_mg_diel
