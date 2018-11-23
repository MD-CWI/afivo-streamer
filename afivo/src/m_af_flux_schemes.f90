!> Module containing a couple simple flux schemes for scalar advection/diffusion
!> problems
!>
!> @todo Explain multidimensional index for flux, velocity
module m_af_flux_schemes

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  public :: flux_diff_1d, flux_diff_2d, flux_diff_3d
  public :: flux_koren_1d, flux_koren_2d, flux_koren_3d

contains

  !> Compute diffusive flux (second order)
  subroutine flux_diff_1d(cc, dc, inv_dr, nc, ngc)
    integer, intent(in)   :: nc               !< Number of cells
    integer, intent(in)   :: ngc              !< Number of ghost cells
    real(dp), intent(in)  :: cc(1-ngc:nc+ngc) !< Cell-centered values
    !> Input: diffusion coeff. at interfaces, output: fluxes
    real(dp), intent(inout)  :: dc(1:nc+1)
    real(dp), intent(in)  :: inv_dr           !< Inverse grid spacing
    integer               :: i

    do i = 1, nc+1
       dc(i) = -dc(i) * (cc(i) - cc(i-1)) * inv_dr
    end do
  end subroutine flux_diff_1d

  !> Compute diffusive flux (second order)
  subroutine flux_diff_2d(cc, dc, inv_dr, nc, ngc)
    integer, intent(in)     :: nc                      !< Number of cells
    integer, intent(in)     :: ngc                     !< Number of ghost cells
    !> Cell-centered values
    real(dp), intent(in)    :: cc(1-ngc:nc+ngc, 1-ngc:nc+ngc)
    !> Input: diffusion coeff. at interfaces, output: fluxes
    real(dp), intent(inout)    :: dc(1:nc+1, 1:nc+1, 2)
    real(dp), intent(in)    :: inv_dr           !< Inverse grid spacing
    real(dp)                :: cc_1d(1-ngc:nc+ngc), dc_1d(1:nc+1)
    integer                 :: n

    do n = 1, nc
       ! x-fluxes
       call flux_diff_1d(cc(:, n), dc(:, n, 1), inv_dr, nc, ngc)

       ! y-fluxes (use temporary variables for efficiency)
       cc_1d = cc(n, :)
       dc_1d = dc(n, :, 2)
       call flux_diff_1d(cc_1d, dc_1d, inv_dr, nc, ngc)
       dc(n, :, 2) = dc_1d     ! Copy result
    end do
  end subroutine flux_diff_2d

  !> Compute diffusive flux (second order)
  subroutine flux_diff_3d(cc, dc, inv_dr, nc, ngc)
    !> Number of cells
    integer, intent(in)     :: nc
    !> Number of ghost cells
    integer, intent(in)     :: ngc
    !> Cell-centered values
    real(dp), intent(in)    :: cc(1-ngc:nc+ngc, 1-ngc:nc+ngc, 1-ngc:nc+ngc)
    !> Input: diffusion coeff. at interfaces, output: fluxes
    real(dp), intent(inout)    :: dc(1:nc+1, 1:nc+1, 1:nc+1, 3)
    !> Inverse grid spacing
    real(dp), intent(in)    :: inv_dr
    real(dp)                :: cc_1d(1-ngc:nc+ngc), dc_1d(1:nc+1)
    integer                 :: n, m

    do m = 1, nc
       do n = 1, nc
          ! x-fluxes
          call flux_diff_1d(cc(:, n, m), dc(:, n, m, 1), &
               inv_dr, nc, ngc)

          ! y-fluxes (use temporary variables for efficiency)
          cc_1d = cc(n, :, m)
          dc_1d = dc(n, :, m, 2)
          call flux_diff_1d(cc_1d, dc_1d, inv_dr, nc, ngc)
          dc(n, :, m, 2) = dc_1d ! Copy result

          ! z-fluxes (use temporary variables for efficiency)
          cc_1d = cc(n, m, :)
          dc_1d = dc(n, m, :, 3)
          call flux_diff_1d(cc_1d, dc_1d, inv_dr, nc, ngc)
          dc(n, m, :, 3) = dc_1d ! Copy result
       end do
    end do
  end subroutine flux_diff_3d

  !> Compute flux according to Koren limiter
  subroutine flux_koren_1d(cc, v, nc, ngc)
    integer, intent(in)   :: nc               !< Number of cells
    integer, intent(in)   :: ngc              !< Number of ghost cells
    real(dp), intent(in)  :: cc(1-ngc:nc+ngc) !< Cell-centered values
    !> Input: velocities at interfaces, output: fluxes
    real(dp), intent(inout)  :: v(1:nc+1)
    real(dp)              :: gradp, gradc, gradn
    integer               :: n

    do n = 1, nc+1
       gradc = cc(n) - cc(n-1)  ! Current gradient
       if (v(n) < 0.0_dp) then
          gradn = cc(n+1) - cc(n) ! Next gradient
          v(n) = v(n) * (cc(n) - koren_mlim(gradc, gradn))
       else                     ! v(n) > 0
          gradp = cc(n-1) - cc(n-2) ! Previous gradient
          v(n) = v(n) * (cc(n-1) + koren_mlim(gradc, gradp))
       end if
    end do

  end subroutine flux_koren_1d

  !> Compute flux according to Koren limiter
  subroutine flux_koren_2d(cc, v, nc, ngc)
    integer, intent(in)     :: nc                      !< Number of cells
    integer, intent(in)     :: ngc                     !< Number of ghost cells
    !> Cell-centered values
    real(dp), intent(in)    :: cc(1-ngc:nc+ngc, 1-ngc:nc+ngc)
    !> Input: velocities at interfaces, output: fluxes
    real(dp), intent(inout)    :: v(1:nc+1, 1:nc+1, 2)
    real(dp)                :: cc_1d(1-ngc:nc+ngc), v_1d(1:nc+1)
    integer                 :: n

    do n = 1, nc
       ! x-fluxes
       call flux_koren_1d(cc(:, n), v(:, n, 1), nc, ngc)

       ! y-fluxes (use temporary variables for efficiency)
       cc_1d = cc(n, :)
       v_1d  = v(n, :, 2)
       call flux_koren_1d(cc_1d, v_1d, nc, ngc)
       v(n, :, 2) = v_1d     ! Copy result
    end do
  end subroutine flux_koren_2d

  !> Compute flux according to Koren limiter
  subroutine flux_koren_3d(cc, v, nc, ngc)
    !> Number of cells
    integer, intent(in)     :: nc
    !> Number of ghost cells
    integer, intent(in)     :: ngc
    !> Cell-centered values
    real(dp), intent(in)    :: cc(1-ngc:nc+ngc, 1-ngc:nc+ngc, 1-ngc:nc+ngc)
    !> Input: velocities at interfaces, output: fluxes
    real(dp), intent(inout)    :: v(1:nc+1, 1:nc+1, 1:nc+1, 3)
    real(dp)                :: cc_1d(1-ngc:nc+ngc), v_1d(1:nc+1)
    integer                 :: n, m

    do m = 1, nc
       do n = 1, nc
          ! x-fluxes
          call flux_koren_1d(cc(:, n, m), v(:, n, m, 1), &
               nc, ngc)

          ! y-fluxes (use temporary variables for efficiency)
          cc_1d = cc(n, :, m)
          v_1d  = v(n, :, m, 2)
          call flux_koren_1d(cc_1d, v_1d, nc, ngc)
          v(n, :, m, 2) = v_1d ! Copy result

          ! z-fluxes (use temporary variables for efficiency)
          cc_1d = cc(n, m, :)
          v_1d  = v(n, m, :, 3)
          call flux_koren_1d(cc_1d, v_1d, nc, ngc)
          v(n, m, :, 3) = v_1d ! Copy result
       end do
    end do
  end subroutine flux_koren_3d

  !> Modified implementation of Koren limiter, to avoid division and the min/max
  !> functions, which can be problematic / expensive. In most literature, you
  !> have r = a / b (ratio of gradients). Then the limiter phi(r) is multiplied
  !> with b. With this implementation, you get phi(r) * b
  elemental function koren_mlim(a, b) result(bphi)
    real(dp), intent(in) :: a  !< Density gradient (numerator)
    real(dp), intent(in) :: b  !< Density gradient (denominator)
    real(dp), parameter  :: sixth = 1/6.0_dp
    real(dp)             :: bphi, aa, ab

    aa = a * a
    ab = a * b

    if (ab <= 0) then
       ! a and b have different sign or one of them is zero, so r is either 0,
       ! inf or negative (special case a == b == 0 is ignored)
       bphi = 0
    else if (aa <= 0.25_dp * ab) then
       ! 0 < a/b <= 1/4, limiter has value a/b
       bphi = a
    else if (aa <= 2.5_dp * ab) then
       ! 1/4 < a/b <= 2.5, limiter has value (1+2*a/b)/6
       bphi = sixth * (b + 2*a)
    else
       ! (1+2*a/b)/6 >= 1, limiter has value 1
       bphi = b
    end if
  end function koren_mlim

end module m_af_flux_schemes
