#include "../afivo/src/cpp_macros_$Dd.h"
!> Module containing a couple simple flux schemes for scalar advection/diffusion
!> problems
!>
!> @todo Explain multidimensional index for flux, velocity
module m_flux_schemes_$Dd
  use m_streamer
  use m_units_constants
  use m_a$D_all
  
  implicit none
  private


  public :: flux_diff_$Dd
  public :: flux_koren_1d, flux_koren_2d, flux_koren_3d
  
  public :: gradient_$Dd

contains

  !> Compute diffusive flux (second order)
  !> Compute diffusive flux (second order)
  subroutine flux_diff_$Dd(box, dc, iv)
    type(box$D_t), intent(in)  :: box
    integer, intent(in)        :: iv
    real(dp), intent(inout)    :: dc(DTIMES(1:box%n_cell+1), $D)
    integer                    :: nc, d, IJK

    nc = box%n_cell

    do KJI_DO(1, nc+1)
       do d = 1, $D  
          dc(IJK, d) = - dc(IJK, d) * gradient_$Dd(box, 1, [IJK], 2*d-1, iv)
       end do
    end do; CLOSE_DO
  end subroutine flux_diff_$Dd

  !> Compute flux according to Koren limiter
  subroutine flux_koren_1d(cc, v, nc, inv_dr)
    integer, intent(in)     :: nc               !< Number of cells
    real(dp), intent(in)    :: cc(-1:nc+2) !< Cell-centered values
    !> Input: velocities at interfaces, output: fluxes
    real(dp), intent(inout) :: v(1:nc+1)
    real(dp), intent(in)    :: inv_dr
    real(dp)                :: gradp, gradc, gradn
    integer                 :: n


    
    do n = 1, nc+1
       gradc = cc(n) - cc(n-1)  ! Current gradient
       if (v(n) < 0.0_dp) then
         gradn = cc(n+1) - cc(n)  ! Next gradient
         v(n) = v(n) * (cc(n) - koren_mlim(gradc, gradn))
       else                    
         gradp = cc(n-1) - cc(n-2)  ! Previous gradient
         v(n) = v(n) * (cc(n-1) + koren_mlim(gradc, gradp))
       end if
    end do

  end subroutine flux_koren_1d

  !> Compute flux according to Koren limiter
  subroutine flux_koren_2d(cc, v, nc, inv_dr)
    integer, intent(in)     :: nc                      !< Number of cells
    real(dp), intent(in)    :: cc(-1:nc+2, -1:nc+2)  !< Cell-centered values
    !> Input: velocities at interfaces, output: fluxes
    real(dp), intent(inout) :: v(1:nc+1, 1:nc+1, 2)
    real(dp), intent(in)    :: inv_dr
    real(dp)                :: cc_1d(-1:nc+2), v_1d(1:nc+1)
    integer                 :: n

    do n = 1, nc
       ! x-fluxes
       call flux_koren_1d(cc(:, n), v(:, n, 1), nc, inv_dr)
       ! y-fluxes (use temporary variables for efficiency)
       cc_1d    = cc(n, :)
       v_1d     = v(n, :, 2)
       call flux_koren_1d(cc_1d, v_1d, nc, inv_dr)
       v(n, :, 2) = v_1d     ! Copy result
    end do
  end subroutine flux_koren_2d

  !> Compute flux according to Koren limiter
  subroutine flux_koren_3d(cc, v, nc, inv_dr)
    integer, intent(in)     :: nc                      !< Number of cells
    !> Cell-centered values
    real(dp), intent(in)    :: cc(-1:nc+2, -1:nc+2, -1:nc+2)  !< Cell-centered values
    !> Input: velocities at interfaces, output: fluxes
    real(dp), intent(inout) :: v(1:nc+1, 1:nc+1, 1:nc+1, 3)
    real(dp), intent(in)    :: inv_dr
    real(dp)                :: cc_1d(-1:nc+2), v_1d(1:nc+1)
    integer                 :: n, m

    do m = 1, nc
       do n = 1, nc
          ! x-fluxes
          call flux_koren_1d(cc(:, n, m), v(:, n, m, 1), nc, inv_dr)

          ! y-fluxes (use temporary variables for efficiency)
          cc_1d    = cc(n, :, m)
          v_1d     = v(n, :, m, 2)
          call flux_koren_1d(cc_1d, v_1d, nc, inv_dr)
          v(n, :, m, 2) = v_1d ! Copy result

          ! z-fluxes (use temporary variables for efficiency)
          cc_1d    = cc(n, m, :)
          v_1d     = v(n, m, :, 3)
          call flux_koren_1d(cc_1d, v_1d, nc, inv_dr)
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
  
  !> Return the gradient off cell centered variable iv in a direction dir from inside cell ix
  function gradient_$Dd(box, cc_num, ix, dir, iv) result(grad)
    type(box$D_t), intent(in)         :: box
    integer, intent(in)               :: cc_num 
    integer, intent(in)               :: ix($D), dir, iv
    real(dp)                          :: grad, inv_dr, s_C, harm_ep, fac, f, f_n, theta
    integer                           :: n_ix($D)
    
    n_ix      = a$D_neighb_dix(:, dir)
    inv_dr    = 1/box%dr
 

#if $D == 2
    
    if (cc_num == 2) then
               
      harm_ep = 2.0_dp*box%cc(ix(1), ix(2), i_eps)*box%cc(ix(1)+n_ix(1), ix(2)+n_ix(2), i_eps) / &
                (box%cc(ix(1), ix(2), i_eps)+box%cc(ix(1)+n_ix(1), ix(2)+n_ix(2), i_eps))           
      s_C     = box%fc(ix(1)+n_ix(1)*a2_neighb_high_01(dir), ix(2)+n_ix(2)*a2_neighb_high_01(dir) &
                , a2_neighb_dim(dir), sigma_rhs)/box%cc(ix(1)+n_ix(1), ix(2)+n_ix(2), i_eps)
      
      
      grad    = a2_neighb_high_pm(dir)*(inv_dr*(box%cc(ix(1)+n_ix(1), ix(2)+n_ix(2), iv) - &
                box%cc(ix(1), ix(2), iv)) - 0.5_dp*s_C) * &
                harm_ep/box%cc(ix(1), ix(2), i_eps) 
    else if (cc_num == 1) then
      !fac   = 0.5_dp * box%fc(ix(1)+n_ix(1)*a2_neighb_high_01(dir), ix(2)+n_ix(2)*a2_neighb_high_01(dir),&
      !        a2_neighb_dim(dir), sigma_rhs)/(box%fc(ix(1)+n_ix(1)*a2_neighb_high_01(dir), ix(2) + &
      !        n_ix(2)*a2_neighb_high_01(dir), a2_neighb_dim(dir), sigma_rhs) + &
      !        box%cc(ix(1)+n_ix(1), ix(2)+n_ix(2), iv)/inv_dr + epsilon(1.0_dp))
        
      grad  = inv_dr*a2_neighb_high_pm(dir)*(box%cc(ix(1)+n_ix(1), ix(2)+n_ix(2), iv)-box%cc(ix(1), ix(2), iv))

    end if



#elif $D == 3


      if (cc_num == 2) then
        harm_ep = 2.0_dp*box%cc(ix(1), ix(2), ix(3), i_eps)*box%cc(ix(1)+n_ix(1), ix(2)+n_ix(2), ix(3)+n_ix(3), i_eps)/&
                  (box%cc(ix(1), ix(2), ix(3), i_eps)+box%cc(ix(1)+n_ix(1), ix(2)+n_ix(2), ix(3)+n_ix(3), i_eps))           
        s_C     = box%fc(ix(1)+n_ix(1)*a3_neighb_high_01(dir), ix(2)+n_ix(2)*a3_neighb_high_01(dir),&
                  ix(3)+n_ix(3)*a3_neighb_high_01(dir), a3_neighb_dim(dir), sigma_rhs)
        grad    = (inv_dr*a3_neighb_high_pm(dir)*(box%cc(ix(1)+n_ix(1), ix(2)+n_ix(2), ix(3)+n_ix(3), iv) -&
                  box%cc(ix(1), ix(2), ix(3), iv)) - 0.5_dp*s_C/box%cc(ix(1)+n_ix(1), ix(2)+n_ix(2), ix(3)+n_ix(3), i_eps))&
                  * harm_ep/box%cc(ix(1), ix(2), ix(3), i_eps)
      else if (cc_num == 1) then
        fac   = 0.5_dp * box%fc(ix(1)+n_ix(1)*a3_neighb_high_01(dir), ix(2)+n_ix(2)*a3_neighb_high_01(dir),&
                ix(3)+n_ix(3)*a3_neighb_high_01(dir), a3_neighb_dim(dir), sigma_rhs)/ &
                (box%fc(ix(1)+n_ix(1)* a3_neighb_high_01(dir), ix(2)+n_ix(2)*a3_neighb_high_01(dir), &
                ix(3)+n_ix(3)*a3_neighb_high_01(dir), a3_neighb_dim(dir), sigma_rhs) + &
                box%cc(ix(1)+n_ix(1), ix(2)+n_ix(2), ix(3)+n_ix(3), iv)/inv_dr + epsilon(1.0_dp))
        
        grad  = inv_dr*a3_neighb_high_pm(dir)*(box%fc(ix(1)+n_ix(1)*a3_neighb_high_01(dir), ix(2)+&
                n_ix(2)*a3_neighb_high_01(dir), ix(3)+n_ix(3)*a3_neighb_high_01(dir), a3_neighb_dim(dir), sigma_rhs)*inv_dr &
                +box%cc(ix(1)+n_ix(1), ix(2)+n_ix(2), ix(3)+n_ix(3), iv)-box%cc(ix(1), ix(2), ix(3), iv)) / (1-fac)
      end if

#endif

  end function gradient_$Dd

end module m_flux_schemes_$Dd
