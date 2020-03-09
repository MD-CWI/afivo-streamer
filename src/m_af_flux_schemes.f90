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
  public :: flux_upwind_1d, flux_upwind_2d, flux_upwind_3d
  public :: flux_kurganovTadmor_1d, reconstruct_lr_1d

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
    real(dp), intent(in)    :: inv_dr(2)           !< Inverse grid spacing
    real(dp)                :: cc_1d(1-ngc:nc+ngc), dc_1d(1:nc+1)
    integer                 :: n

    do n = 1, nc
       ! x-fluxes
       call flux_diff_1d(cc(:, n), dc(:, n, 1), inv_dr(1), nc, ngc)

       ! y-fluxes (use temporary variables for efficiency)
       cc_1d = cc(n, :)
       dc_1d = dc(n, :, 2)
       call flux_diff_1d(cc_1d, dc_1d, inv_dr(2), nc, ngc)
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
    real(dp), intent(in)    :: inv_dr(3)
    real(dp)                :: cc_1d(1-ngc:nc+ngc), dc_1d(1:nc+1)
    integer                 :: n, m

    do m = 1, nc
       do n = 1, nc
          ! x-fluxes
          call flux_diff_1d(cc(:, n, m), dc(:, n, m, 1), &
               inv_dr(1), nc, ngc)

          ! y-fluxes (use temporary variables for efficiency)
          cc_1d = cc(n, :, m)
          dc_1d = dc(n, :, m, 2)
          call flux_diff_1d(cc_1d, dc_1d, inv_dr(2), nc, ngc)
          dc(n, :, m, 2) = dc_1d ! Copy result

          ! z-fluxes (use temporary variables for efficiency)
          cc_1d = cc(n, m, :)
          dc_1d = dc(n, m, :, 3)
          call flux_diff_1d(cc_1d, dc_1d, inv_dr(3), nc, ngc)
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

  subroutine reconstruct_lr_1d(nc, ngc, n_vars, cc, u_lr)
    integer, intent(in)     :: nc               !< Number of cells
    integer, intent(in)     :: ngc              !< Number of ghost cells
    integer, intent(in)     :: n_vars            !< Number of variables
    real(dp), intent(in)    :: cc(1-ngc:nc+ngc, n_vars) !< Cell-centered values
    !> Reconstructed (left, right) values at every interface
    real(dp), intent(inout) :: u_lr(1:nc+1, 2, n_vars)
    real(dp)                :: grad1(n_vars), grad2(n_vars), grad3(n_vars)
    integer                 :: n

    do n=1, nc+1
      grad1 = cc(n-1, :) - cc(n-2, :)
      grad2 = cc(n, :) - cc(n-1, :)
      grad3 = cc(n+1, :) - cc(n, :)
      u_lr(n,1,:) = cc(n-1,:) + 0.5_dp*koren_mlim(grad1, grad2) !left
      u_lr(n,2,:) = cc(n,:) - 0.5_dp*koren_mlim(grad2, grad3) !right
      !vint(n) = vint(n)*(cc(n-1) + 0.5_dp*vanLeer_mlim(grad1, grad2)) !left
      !vout(n) = vout(n)*(cc(n) - 0.5_dp*vanLeer_mlim(grad2, grad3)) !right
      !u_lr(n, 1, :) = cc(n-1, :) + 0.5_dp*minmod_mlim(grad1, grad2) !left
      !u_lr(n, 2, :) = cc(n, :) - 0.5_dp*minmod_mlim(grad2, grad3) !right
    end do
  end subroutine reconstruct_lr_1d

  !> Compute flux for the KT scheme
  subroutine flux_kurganovTadmor_1d( nf, n_vars, flux_lr, u_lr, w_lr, flux )
    !integer, intent(in) :: nc
    integer, intent(in) :: nf
    integer, intent(in) :: n_vars
    real(dp), intent(in) :: flux_lr(nf, 2, n_vars), u_lr(nf,2,n_vars)
    real(dp), intent(in) :: w_lr(nf)
    real(dp), intent(inout) :: flux(nf, n_vars)

    flux = 0.5_dp*(flux_lr(:,1,:) + flux_lr(:,2,:)) - &
                spread(w_lr, 2,n_vars)*(u_lr(:,2,:) - u_lr(:,1,:))
  end subroutine flux_kurganovTadmor_1d

  ! Jannis: Maybe this is not such a great idea, since the overhead of the
  ! function calls can become significant in this line-by-line algorithm.
  !
  ! subroutine flux_generic(nc, n_gc, n_vars, cc_line, flux, flux_dim, &
  !      to_primitive, to_conservative, &
  !      max_wavespeed, flux_from_primitives)
  !   integer, intent(in)     :: nc
  !   integer, intent(in)     :: n_gc
  !   integer, intent(in)     :: n_vars
  !   real(dp), intent(in)    :: cc_line(1-n_gc:nc+n_gc, n_vars)
  !   real(dp), intent(out)    :: flux(1:nc+1, n_vars)
  !   integer, intent(in)     :: flux_dim
  !   procedure(subr_prim_cons) :: to_primitive, to_conservative
  !   procedure(subr_max_wavespeed) :: max_wavespeed
  !   procedure(subr_flux_from_prim) :: flux_from_primitives

  !   real(dp) :: prim_vars(1-n_gc:nc+n_gc, n_vars)
  !   real(dp) :: u_lr(1:nc+1, 2, n_vars)
  !   real(dp) :: w_lr(1:nc+1)
  !   real(dp) :: flux_lr(1:nc+1, 2, n_vars)
  !   integer  :: i, j

  !   if (n_gc /= 2) error stop "n_gc /= 2"

  !   !We need primitives at the ghost cells as they're used for reconstruction
  !   call to_primitive(nc+2*n_gc, n_vars, cc_line, prim_vars)
  !   call reconstruct_lr_1d(nc, n_gc, n_vars, prim_vars, u_lr)
  !   call get_max_wavespeed_lr_1d(nc+1, n_vars, 1, u_lr, w_lr)
  !   call get_fluxes_lr_1d(nc+1, n_vars, 1, u_lr, flux_lr)
  !   call to_conservatives(nc+1, n_vars, u_lr, cons_vars)
  !   call flux_kurganovTadmor_1d(nc+1, n_vars, flux_lr, u_lr, w_lr, flux)
  !   tree%boxes(id)%fc(:, j, 1, :) = flux

  ! end subroutine flux_generic

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

  elemental function vanLeer_mlim(a, b) result(phi)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp) :: phi
    phi = (2.0_dp*max(0.0_dp, a*b))/(a + b + epsilon(1.0_dp))

  end function vanLeer_mlim


  elemental function minmod_mlim(a, b) result(phi)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp)             :: phi
    phi = 0.5_dp*(sign(1.0_dp, a) + sign(1.0_dp, b)) * &
          min(abs(a), abs(b))
  end function minmod_mlim

  !> Compute flux with first order upwind scheme
  subroutine flux_upwind_1d(cc, v, nc, ngc)
    integer, intent(in)   :: nc               !< Number of cells
    integer, intent(in)   :: ngc              !< Number of ghost cells
    real(dp), intent(in)  :: cc(1-ngc:nc+ngc) !< Cell-centered values
    !> Input: velocities at interfaces, output: fluxes
    real(dp), intent(inout)  :: v(1:nc+1)
    integer               :: n

    do n = 1, nc+1
       if (v(n) < 0.0_dp) then
          v(n) = v(n) * cc(n)
       else                     ! v(n) > 0
          v(n) = v(n) * cc(n-1)
       end if
    end do

  end subroutine flux_upwind_1d

  !> Compute flux with first order upwind scheme
  subroutine flux_upwind_2d(cc, v, nc, ngc)
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
       call flux_upwind_1d(cc(:, n), v(:, n, 1), nc, ngc)

       ! y-fluxes (use temporary variables for efficiency)
       cc_1d = cc(n, :)
       v_1d  = v(n, :, 2)
       call flux_upwind_1d(cc_1d, v_1d, nc, ngc)
       v(n, :, 2) = v_1d     ! Copy result
    end do
  end subroutine flux_upwind_2d

  !> Compute flux with first order upwind scheme
  subroutine flux_upwind_3d(cc, v, nc, ngc)
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
          call flux_upwind_1d(cc(:, n, m), v(:, n, m, 1), &
               nc, ngc)

          ! y-fluxes (use temporary variables for efficiency)
          cc_1d = cc(n, :, m)
          v_1d  = v(n, :, m, 2)
          call flux_upwind_1d(cc_1d, v_1d, nc, ngc)
          v(n, :, m, 2) = v_1d ! Copy result

          ! z-fluxes (use temporary variables for efficiency)
          cc_1d = cc(n, m, :)
          v_1d  = v(n, m, :, 3)
          call flux_upwind_1d(cc_1d, v_1d, nc, ngc)
          v(n, m, :, 3) = v_1d ! Copy result
       end do
    end do
  end subroutine flux_upwind_3d

end module m_af_flux_schemes
