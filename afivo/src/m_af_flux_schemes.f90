#include "cpp_macros.h"
!> Module containing a couple flux schemes for solving hyperbolic problems
!> explicitly, as well as handling diffusion explicitly.
module m_af_flux_schemes
  use m_af_types

  implicit none
  private

  integer, parameter :: limiter_vanleer = 1
  integer, parameter :: limiter_koren = 2
  integer, parameter :: limiter_minmod = 3

  interface
     subroutine subr_prim_cons(n_values, n_vars, u)
       import
       integer, intent(in)     :: n_values, n_vars
       real(dp), intent(inout) :: u(n_values, n_vars)
     end subroutine subr_prim_cons

     subroutine subr_max_wavespeed(nf, n_var, flux_dim, u, w)
       import
       integer, intent(in)   :: nf    !< Number of cell faces
       integer, intent(in)   :: n_var !< Number of variables
       integer, intent(in)   :: flux_dim !< In which dimension fluxes are computed
       real(dp), intent(in)  :: u(nf, n_var)
       real(dp), intent(out) :: w(nf)
     end subroutine subr_max_wavespeed

     subroutine subr_flux_from_prim(nf, n_var, flux_dim, u, flux, box, line_ix)
       import
       integer, intent(in)     :: nf              !< Number of cell faces
       integer, intent(in)     :: n_var           !< Number of variables
       integer, intent(in)     :: flux_dim        !< In which dimension fluxes are computed
       real(dp), intent(in)    :: u(nf, n_var)    !< Primitive variables
       real(dp), intent(out)   :: flux(nf, n_var) !< Computed fluxes
       type(box_t), intent(in) :: box             !< Current box
       integer, intent(in)     :: line_ix(NDIM-1) !< Index of line for dim /= flux_dim
     end subroutine subr_flux_from_prim

     subroutine subr_source(box, dt, n_vars, i_cc, s_deriv, s_out)
       import
       type(box_t), intent(inout) :: box
       real(dp), intent(in)       :: dt
       integer, intent(in)        :: n_vars
       integer, intent(in)        :: i_cc(n_vars)
       integer, intent(in)        :: s_deriv
       integer, intent(in)        :: s_out
     end subroutine subr_source
  end interface

  public :: flux_diff_1d, flux_diff_2d, flux_diff_3d
  public :: flux_koren_1d, flux_koren_2d, flux_koren_3d
  public :: flux_upwind_1d, flux_upwind_2d, flux_upwind_3d
  public :: flux_kurganovTadmor_1d, reconstruct_lr_1d
  public :: flux_generic_box, flux_generic_tree
  public :: flux_update_densities

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

  subroutine reconstruct_lr_1d(nc, ngc, n_vars, cc, u_l, u_r, limiter)
    integer, intent(in)           :: nc                       !< Number of cells
    integer, intent(in)           :: ngc                      !< Number of ghost cells
    integer, intent(in)           :: n_vars                   !< Number of variables
    real(dp), intent(in)          :: cc(1-ngc:nc+ngc, n_vars) !< Cell-centered values
    !> Reconstructed (left, right) values at every interface
    real(dp), intent(inout)       :: u_l(1:nc+1, n_vars)
    real(dp), intent(inout)       :: u_r(1:nc+1, n_vars)
    integer, intent(in), optional :: limiter                  !< Which limiter to use
    real(dp)                      :: slopes(0:nc+1, n_vars)
    integer                       :: use_limiter

    use_limiter = limiter_vanleer
    if (present(limiter)) use_limiter = limiter

    ! Compute limited slopes at the cell centers
    select case (use_limiter)
    case (limiter_vanleer)
       slopes = vanLeer_mlim(cc(1:nc+2, :) - cc(0:nc+1, :), &
            cc(0:nc+1, :) - cc(-1:nc, :))
    case (limiter_koren)
       slopes = koren_mlim(cc(1:nc+2, :) - cc(0:nc+1, :), &
            cc(0:nc+1, :) - cc(-1:nc, :))
    case (limiter_minmod)
       slopes = minmod_mlim(cc(1:nc+2, :) - cc(0:nc+1, :), &
            cc(0:nc+1, :) - cc(-1:nc, :))
    case default
       error stop "unknown limiter"
    end select

    ! Reconstruct values on the faces
    u_l(1:nc+1, :) = cc(0:nc, :) + 0.5_dp * slopes(0:nc, :) ! left
    u_r(1:nc+1, :) = cc(1:nc+1, :) - 0.5_dp * slopes(1:nc+1, :) ! right
  end subroutine reconstruct_lr_1d

  !> Compute flux for the KT scheme
  subroutine flux_kurganovTadmor_1d(n_values, n_vars, flux_l, flux_r, &
       u_l, u_r, wmax, flux)
    integer, intent(in) :: n_values
    integer, intent(in) :: n_vars
    real(dp), intent(in) :: flux_l(n_values, n_vars)
    real(dp), intent(in) :: flux_r(n_values, n_vars)
    real(dp), intent(in) :: u_l(n_values, n_vars)
    real(dp), intent(in) :: u_r(n_values, n_vars)
    real(dp), intent(in) :: wmax(n_values)
    real(dp), intent(inout) :: flux(n_values, n_vars)

    flux = 0.5_dp * (flux_l + flux_r - spread(wmax, 2, n_vars) * (u_r - u_l))
  end subroutine flux_kurganovTadmor_1d

  subroutine flux_update_densities(tree, dt, n_vars, i_cc, i_flux, &
       s_deriv, s_prev, s_out, add_source_box)
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt             !< Time step
    integer, intent(in)       :: n_vars         !< Number of variables
    integer, intent(in)       :: i_cc(n_vars)   !< Cell-centered variables
    integer, intent(in)       :: i_flux(n_vars) !< Flux variables
    integer, intent(in)       :: s_deriv        !< State to compute derivatives from
    integer, intent(in)       :: s_prev         !< Previous time state
    integer, intent(in)       :: s_out          !< Output time state
    !> Method to include source terms
    procedure(subr_source), optional :: add_source_box
    integer                   :: lvl, n, id, IJK, nc
    real(dp)                  :: dt_dr(NDIM)
    real(dp)                  :: rfac(2, tree%n_cell)

    nc = tree%n_cell
    rfac = 0.0_dp ! Prevent warnings in 3D

    !$omp parallel private(lvl, n, id, IJK, dt_dr, rfac)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do n = 1, size(tree%lvls(lvl)%leaves)
          id    = tree%lvls(lvl)%leaves(n)
          dt_dr = dt/tree%boxes(id)%dr

          tree%boxes(id)%cc(DTIMES(1:nc), i_cc+s_out) = &
               tree%boxes(id)%cc(DTIMES(1:nc), i_cc+s_prev)

          if (present(add_source_box)) then
             call add_source_box(tree%boxes(id), dt, &
                  n_vars, i_cc, s_deriv, s_out)
          end if

          associate(cc => tree%boxes(id)%cc, fc => tree%boxes(id)%fc)
#if NDIM == 1
            do KJI_DO(1, nc)
               cc(i, i_cc+s_out) = cc(i, i_cc+s_out) + &
                    dt_dr(1) * &
                    (fc(i, 1, i_flux) - fc(i+1, 1, i_flux))
            end do; CLOSE_DO
#elif NDIM == 2
            if (tree%coord_t == af_cyl) then
               call af_cyl_flux_factors(tree%boxes(id), rfac)
               do KJI_DO(1, nc)
                  cc(i, j, i_cc+s_out) = cc(i, j, i_cc+s_out) + &
                       dt_dr(1) * (&
                       rfac(1, i) * fc(i, j, 1, i_flux) - &
                       rfac(2, i) * fc(i+1, j, 1, i_flux)) &
                       + dt_dr(2) * &
                       (fc(i, j, 2, i_flux) - fc(i, j+1, 2, i_flux))
               end do; CLOSE_DO
            else
               do KJI_DO(1, nc)
                  cc(i, j, i_cc+s_out) = cc(i, j, i_cc+s_out) + &
                       dt_dr(1) * &
                       (fc(i, j, 1, i_flux) - fc(i+1, j, 1, i_flux)) &
                       + dt_dr(2) * &
                       (fc(i, j, 2, i_flux) - fc(i, j+1, 2, i_flux))
               end do; CLOSE_DO
            end if
#elif NDIM == 3
            do KJI_DO(1, nc)
               cc(i, j, k, i_cc+s_out) = cc(i, j, k, i_cc+s_out) + &
                    dt_dr(1) * &
                    (fc(i, j, k, 1, i_flux) - fc(i+1, j, k, 1, i_flux)) + &
                    dt_dr(2) * &
                    (fc(i, j, k, 2, i_flux) - fc(i, j+1, k, 2, i_flux)) + &
                    dt_dr(3) * &
                    (fc(i, j, k, 3, i_flux) - fc(i, j, k+1, 3, i_flux))
            end do; CLOSE_DO
#endif
          end associate
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine flux_update_densities

  subroutine flux_generic_tree(tree, n_vars, i_cc, i_flux, wmax, &
       max_wavespeed, flux_from_primitives, to_primitive, to_conservative)
    use m_af_restrict
    use m_af_core
    type(af_t), intent(inout)      :: tree
    integer, intent(in)            :: n_vars         !< Number of variables
    integer, intent(in)            :: i_cc(n_vars)   !< Cell-centered variables
    integer, intent(in)            :: i_flux(n_vars) !< Flux variables
    real(dp), intent(out)          :: wmax(NDIM)     !< Maximum wave speed found
    !> Compute the maximum wave speed
    procedure(subr_max_wavespeed)  :: max_wavespeed
    !> Compute the flux from primitive variables
    procedure(subr_flux_from_prim) :: flux_from_primitives
    !> Convert conservative variables to primitive ones
    procedure(subr_prim_cons), optional :: to_primitive
    !> Convert primitive variables to conservative ones
    procedure(subr_prim_cons), optional :: to_conservative

    integer :: lvl, i

    ! Ensure ghost cells near refinement boundaries can be properly filled
    call af_restrict_ref_boundary(tree, i_cc)

    wmax(:) = 0

    !$omp parallel private(lvl, i) reduction(max:wmax)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          call flux_generic_box(tree, tree%lvls(lvl)%leaves(i), tree%n_cell, &
               n_vars, i_cc, i_flux, wmax, max_wavespeed, &
               flux_from_primitives, to_primitive, to_conservative)
       end do
       !$omp end do
    end do
    !$omp end parallel

    ! Compute coarse fluxes from the fine ones at refinement boundaries
    call af_consistent_fluxes(tree, i_flux)

  end subroutine flux_generic_tree

  !> Compute generic finite volume flux
  subroutine flux_generic_box(tree, id, nc, n_vars, i_cc, i_flux, wmax, &
       max_wavespeed, flux_from_primitives, to_primitive, to_conservative)
    use m_af_types
    use m_af_ghostcell
    type(af_t), intent(inout)      :: tree
    integer, intent(in)            :: id             !< Id of box
    integer, intent(in)            :: nc             !< Number of cells
    integer, intent(in)            :: n_vars         !< Number of variables
    integer, intent(in)            :: i_cc(n_vars)   !< Cell-centered variables
    integer, intent(in)            :: i_flux(n_vars) !< Flux variables
    real(dp), intent(inout)        :: wmax(NDIM)     !< Maximum wave speed found
    !> Compute the maximum wave speed
    procedure(subr_max_wavespeed)  :: max_wavespeed
    !> Compute the flux from primitive variables
    procedure(subr_flux_from_prim) :: flux_from_primitives
    !> Convert conservative variables to primitive ones
    procedure(subr_prim_cons), optional :: to_primitive
    !> Convert primitive variables to conservative ones
    procedure(subr_prim_cons), optional :: to_conservative


    real(dp) :: cc(DTIMES(-1:nc+2), n_vars)
    real(dp) :: cc_line(-1:nc+2, n_vars)
    real(dp) :: flux(nc+1, n_vars)
    real(dp) :: u_l(nc+1, n_vars), u_r(nc+1, n_vars)
    real(dp) :: w_l(nc+1), w_r(nc+1)
    real(dp) :: flux_l(nc+1, n_vars), flux_r(nc+1, n_vars)
    integer  :: flux_dim, line_ix(NDIM-1)
#if NDIM > 1
    integer  :: i
#endif
#if NDIM > 2
    integer  :: j
#endif

    ! Get two layers of ghost cell data
    call af_gc2_box(tree, id, i_cc, cc)

    ! Jannis: Below, there are function calls in the inner part of a loop. When
    ! I did some benchmarks, it was not significantly slower than using a buffer
    ! and fewer functions calls.

    do flux_dim = 1, NDIM
#if NDIM > 2
       do j = 1, nc
#endif
#if NDIM > 1
          do i = 1, nc
#endif
             ! Extract cell-centered values along a line
             select case (flux_dim)
#if NDIM == 1
             case (1)
                cc_line = cc(:, :)
#elif NDIM == 2
             case (1)
                cc_line = cc(:, i, :)
             case (2)
                cc_line = cc(i, :, :)
#elif NDIM == 3
             case (1)
                cc_line = cc(:, i, j, :)
             case (2)
                cc_line = cc(i, :, j, :)
             case (3)
                cc_line = cc(:, i, j, :)
#endif
             end select

#if NDIM == 2
             line_ix = [i]
#elif NDIM == 3
             line_ix = [i, j]
#endif

             if (present(to_primitive)) then
                call to_primitive(nc+4, n_vars, cc_line)
             end if

             ! Reconstruct to cell faces
             call reconstruct_lr_1d(nc, 2, n_vars, cc_line, u_l, u_r)

             call max_wavespeed(nc+1, n_vars, flux_dim, u_l, w_l)
             call max_wavespeed(nc+1, n_vars, flux_dim, u_r, w_r)
             call flux_from_primitives(nc+1, n_vars, flux_dim, u_l, flux_l, &
                  tree%boxes(id), line_ix)
             call flux_from_primitives(nc+1, n_vars, flux_dim, u_r, flux_r, &
                  tree%boxes(id), line_ix)

             if (present(to_conservative)) then
                call to_conservative(nc+1, n_vars, u_l)
                call to_conservative(nc+1, n_vars, u_r)
             end if

             w_l = max(w_l, w_r) ! Get maximum of left/right wave speed
             call flux_kurganovTadmor_1d(nc+1, n_vars, flux_l, flux_r, &
                  u_l, u_r, w_l, flux)

             ! Store maximum wave speed
             wmax(flux_dim) = max(wmax(flux_dim), maxval(w_l))

             ! Store the computed fluxes
             select case (flux_dim)
#if NDIM == 1
             case (1)
                tree%boxes(id)%fc(:, flux_dim, i_flux) = flux
#elif NDIM == 2
             case (1)
                tree%boxes(id)%fc(:, i, flux_dim, i_flux) = flux
             case (2)
                tree%boxes(id)%fc(i, :, flux_dim, i_flux) = flux
#elif NDIM == 3
             case (1)
                tree%boxes(id)%fc(:, i, j, flux_dim, i_flux) = flux
             case (2)
                tree%boxes(id)%fc(i, :, j, flux_dim, i_flux) = flux
             case (3)
                tree%boxes(id)%fc(:, i, j, flux_dim, i_flux) = flux
#endif
             end select
#if NDIM > 1
          end do
#endif
#if NDIM > 2
       end do
#endif
    end do

  end subroutine flux_generic_box

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
    real(dp) :: ab, phi

    ab = a * b
    if (ab > 0) then
       phi = 2 * ab / (a + b)
    else
       phi = 0
    end if
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
