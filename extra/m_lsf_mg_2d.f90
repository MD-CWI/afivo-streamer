
!> Multigrid operators for internal boundary conditions. A level set
!> function defines the location of the interface(s).
!\author Jannis Teunissen \copyright GPLv3

! The following replacements take place on this file (m_mgb_Xd.f90) to generate
! 2D and 3D versions:
! 1. 2 -> 2 or 3 (dimension of code)
! 2. preprocess file with cpp
! 3. cat -s (merge multiple blank lines)
module m_lsf_mg_2d
  use m_afivo_2d
  use m_mg_2d

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  public :: lsf2_box_corr_lpl
  public :: lsf2_box_lpl
  public :: lsf2_box_gsrb_lpl

contains

  subroutine lsf2_box_corr_lpl(box_p, box_c, mg)
    type(box2_t), intent(inout) :: box_c
    type(box2_t), intent(in)    :: box_p
    type(mg2_t), intent(in)     :: mg
    integer                     :: i_phi, i_corr, i_lsf, ix_offset(2)
    integer                     :: nc, i, j, i_c1, i_c2, j_c1, j_c2
    real(dp)                    :: val(3), dist(3), tmp

    nc        = box_c%n_cell
    ix_offset = a2_get_child_offset(box_c)
    i_phi     = mg%i_phi
    i_corr    = mg%i_tmp
    i_lsf     = mg%i_lsf

    ! In these loops, we calculate the closest coarse index (_c1), and the
    ! one-but-closest (_c2). The fine cell lies in between.
    do j = 1, nc
       j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

          call lsf_dist_val(box_c%cc(i, j, i_lsf), &
               box_p%cc(i_c1, j_c1, [i_corr, i_lsf]), 0.0_dp, dist(1), val(1))
          call lsf_dist_val(box_c%cc(i, j, i_lsf), &
               box_p%cc(i_c2, j_c1, [i_corr, i_lsf]), 0.0_dp, dist(2), val(2))
          call lsf_dist_val(box_c%cc(i, j, i_lsf), &
               box_p%cc(i_c1, j_c2, [i_corr, i_lsf]), 0.0_dp, dist(3), val(3))

          tmp = 1 / (dist(1) * dist(2) + dist(1) * dist(3) + 2 * dist(2) * dist(3))
          box_c%cc(i, j, i_phi) = box_c%cc(i, j, i_phi) + &
               2 * dist(2) * dist(3) * tmp * val(1) + &
               dist(1) * dist(3) * tmp * val(2) + &
               dist(1) * dist(2) * tmp * val(3)
       end do
    end do
  end subroutine lsf2_box_corr_lpl

  subroutine lsf_dist_val(lsf_a, v_b, b_value, dist, val)
    real(dp), intent(in)  :: lsf_a, v_b(2), b_value
    real(dp), intent(out) :: dist, val

    ! Determine whether there is a boundary
    if (lsf_a * v_b(2) < 0) then
       dist = lsf_a / (lsf_a - v_b(2))
       val  = b_value
    else
       dist = 1
       val  = v_b(1)
    end if
  end subroutine lsf_dist_val

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator
  subroutine lsf2_box_gsrb_lpl(box, redblack_cntr, mg)
    type(box2_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, i0, j, nc, i_phi, i_rhs, i_lsf
    real(dp)                    :: dx2, dist(4), val(4)

    dx2   = box%dr**2
    nc    = box%n_cell
    i_phi = mg%i_phi
    i_rhs = mg%i_rhs
    i_lsf = mg%i_lsf

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          call lsf_dist_val(box%cc(i, j, i_lsf), &
               box%cc(i-1, j, [i_phi, i_lsf]), 0.0_dp, dist(1), val(1))
          call lsf_dist_val(box%cc(i, j, i_lsf), &
               box%cc(i+1, j, [i_phi, i_lsf]), 0.0_dp, dist(2), val(2))
          call lsf_dist_val(box%cc(i, j, i_lsf), &
               box%cc(i, j-1, [i_phi, i_lsf]), 0.0_dp, dist(3), val(3))
          call lsf_dist_val(box%cc(i, j, i_lsf), &
               box%cc(i, j+1, [i_phi, i_lsf]), 0.0_dp, dist(4), val(4))

          box%cc(i, j, i_phi) = 0.5_dp / &
               (dist(1) * dist(2) + dist(3) * dist(4)) * ( &
               (dist(2) * val(1) + dist(1) * val(2)) * &
               dist(3) * dist(4) / (0.5 * (dist(1) + dist(2))) + &
               (dist(4) * val(3) + dist(3) * val(4)) * &
               dist(1) * dist(2) / (0.5 * (dist(3) + dist(4))) - &
               product(dist) * dx2 * box%cc(i, j, i_rhs))
       end do
    end do
  end subroutine lsf2_box_gsrb_lpl

  !> Perform Laplacian operator on a box
  subroutine lsf2_box_lpl(box, i_out, mg)
    type(box2_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, j, nc, i_phi, i_lsf
    real(dp)                    :: inv_dr_sq, dist(4), val(4), f0

    nc        = box%n_cell
    inv_dr_sq = 1 / box%dr**2
    i_phi     = mg%i_phi
    i_lsf     = mg%i_lsf

    do j = 1, nc
       do i = 1, nc
          call lsf_dist_val(box%cc(i, j, i_lsf), &
               box%cc(i-1, j, [i_phi, i_lsf]), 0.0_dp, dist(1), val(1))
          call lsf_dist_val(box%cc(i, j, i_lsf), &
               box%cc(i+1, j, [i_phi, i_lsf]), 0.0_dp, dist(2), val(2))
          call lsf_dist_val(box%cc(i, j, i_lsf), &
               box%cc(i, j-1, [i_phi, i_lsf]), 0.0_dp, dist(3), val(3))
          call lsf_dist_val(box%cc(i, j, i_lsf), &
               box%cc(i, j+1, [i_phi, i_lsf]), 0.0_dp, dist(4), val(4))

          f0 = box%cc(i, j, i_phi)
          box%cc(i, j, i_out) = inv_dr_sq * ( &
               (dist(2) * val(1) + dist(1) * val(2) - (dist(1)+dist(2)) * f0) / &
               (0.5_dp * (dist(1) + dist(2)) * dist(1) * dist(2)) + &
               (dist(4) * val(3) + dist(3) * val(4) - (dist(3)+dist(4)) * f0) / &
               (0.5_dp * (dist(3) + dist(4)) * dist(3) * dist(4)))
       end do
    end do
  end subroutine lsf2_box_lpl

end module m_lsf_mg_2d
