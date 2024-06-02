!> Module containing a couple flux schemes for solving hyperbolic problems
!> explicitly, as well as handling diffusion explicitly.
module m_af_flux_schemes
#include "cpp_macros.h"
  use m_af_types
  use m_af_limiters

  implicit none
  private

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

     subroutine subr_flux(nf, n_var, flux_dim, u, flux, box, line_ix, s_deriv)
       import
       integer, intent(in)     :: nf              !< Number of cell faces
       integer, intent(in)     :: n_var           !< Number of variables
       integer, intent(in)     :: flux_dim        !< In which dimension fluxes are computed
       real(dp), intent(in)    :: u(nf, n_var)    !< Primitive variables
       real(dp), intent(out)   :: flux(nf, n_var) !< Computed fluxes
       type(box_t), intent(in) :: box             !< Current box
       integer, intent(in)     :: line_ix(NDIM-1) !< Index of line for dim /= flux_dim
       integer, intent(in)     :: s_deriv        !< State to compute derivatives from
     end subroutine subr_flux

     subroutine subr_flux_upwind(nf, n_var, flux_dim, u, flux, cfl_sum, &
          n_other_dt, other_dt, box, line_ix, s_deriv)
       import
       integer, intent(in)     :: nf              !< Number of cell faces
       integer, intent(in)     :: n_var           !< Number of variables
       integer, intent(in)     :: flux_dim        !< In which dimension fluxes are computed
       real(dp), intent(in)    :: u(nf, n_var)    !< Face values
       real(dp), intent(out)   :: flux(nf, n_var) !< Computed fluxes
       !> Terms per cell-center to be added to CFL sum, see flux_upwind_box
       real(dp), intent(out)   :: cfl_sum(nf-1)
       integer, intent(in)     :: n_other_dt !< Number of non-cfl time step restrictions
       real(dp), intent(inout) :: other_dt(n_other_dt) !< Non-cfl time step restrictions
       type(box_t), intent(in) :: box             !< Current box
       integer, intent(in)     :: line_ix(NDIM-1) !< Index of line for dim /= flux_dim
       integer, intent(in)     :: s_deriv         !< State to compute derivatives from
     end subroutine subr_flux_upwind

     subroutine subr_flux_modify(nf, n_var, flux_dim, flux, box, line_ix, s_deriv)
       import
       integer, intent(in)     :: nf              !< Number of cell faces
       integer, intent(in)     :: n_var           !< Number of variables
       integer, intent(in)     :: flux_dim        !< In which dimension fluxes are computed
       real(dp), intent(inout) :: flux(nf, n_var) !< Flux that will be modified
       type(box_t), intent(in) :: box             !< Current box
       integer, intent(in)     :: line_ix(NDIM-1) !< Index of line for dim /= flux_dim
       integer, intent(in)     :: s_deriv         !< State to compute derivatives from
     end subroutine subr_flux_modify

     subroutine subr_flux_dir(box, line_ix, s_deriv, n_var, flux_dim, direction_positive)
       import
       type(box_t), intent(in) :: box             !< Current box
       integer, intent(in)     :: line_ix(NDIM-1) !< Index of line for dim /= flux_dim
       integer, intent(in)     :: s_deriv         !< State to compute derivatives from
       integer, intent(in)     :: n_var           !< Number of variables
       integer, intent(in)     :: flux_dim        !< In which dimension fluxes are computed
       !> True means positive flow (to the "right"), false to the left
       logical, intent(out)    :: direction_positive(box%n_cell+1, n_var)
     end subroutine subr_flux_dir

     subroutine subr_line_modify(n_cc, n_var, cc_line, flux_dim, box, line_ix, s_deriv)
       import
       integer, intent(in)     :: n_cc                 !< Number of cell centers
       integer, intent(in)     :: n_var                !< Number of variables
       real(dp), intent(inout) :: cc_line(n_cc, n_var) !< Line values to modify
       integer, intent(in)     :: flux_dim             !< In which dimension fluxes are computed
       type(box_t), intent(in) :: box                  !< Current box
       integer, intent(in)     :: line_ix(NDIM-1)      !< Index of line for dim /= flux_dim
       integer, intent(in)     :: s_deriv              !< State to compute derivatives from
     end subroutine subr_line_modify

     subroutine subr_source(box, dt, n_vars, i_cc, s_deriv, s_out, n_dt, dt_lim, mask)
       import
       type(box_t), intent(inout) :: box
       real(dp), intent(in)       :: dt           !< Time step
       integer, intent(in)        :: n_vars       !< Number of cell-centered variables
       integer, intent(in)        :: i_cc(n_vars) !< Indices of all cell-centered variables
       integer, intent(in)        :: s_deriv      !< State to compute derivatives from
       integer, intent(in)        :: s_out        !< State to update
       !> Number of time steps restrictions
       integer, intent(in)        :: n_dt
       !> Maximal allowed time steps
       real(dp), intent(inout)    :: dt_lim(n_dt)
       !> Mask where to update solution
       logical, intent(in)        :: mask(DTIMES(box%n_cell))
     end subroutine subr_source

     subroutine subr_box_mask(box, mask)
       import
       type(box_t), intent(in) :: box
       logical, intent(out)    :: mask(DTIMES(box%n_cell))
     end subroutine subr_box_mask
  end interface

  public :: flux_diff_1d, flux_diff_2d, flux_diff_3d
  public :: flux_koren_1d, flux_koren_2d, flux_koren_3d
  public :: flux_upwind_1d, flux_upwind_2d, flux_upwind_3d
  public :: flux_kurganovTadmor_1d, reconstruct_lr_1d
  public :: flux_generic_box, flux_generic_tree
  public :: flux_upwind_box, flux_upwind_tree
  public :: flux_update_densities
  public :: flux_dummy_conversion
  public :: flux_dummy_source
  public :: flux_dummy_modify
  public :: flux_dummy_line_modify
  public :: flux_get_line_cc, flux_get_line_1cc
  public :: flux_get_line_fc, flux_get_line_1fc

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
          v(n) = v(n) * (cc(n) - 0.5_dp * af_limiter_koren(gradc, gradn))
       else                     ! v(n) > 0
          gradp = cc(n-1) - cc(n-2) ! Previous gradient
          v(n) = v(n) * (cc(n-1) + 0.5_dp * af_limiter_koren(gradc, gradp))
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

  !> Reconstruct "left" and "right" values at cell faces from cell-centered
  !> data. The left value is for right-going waves, and the right value for
  !> left-going waves.
  subroutine reconstruct_lr_1d(nc, ngc, n_vars, cc, u_l, u_r, limiter)
    integer, intent(in)     :: nc                       !< Number of cells
    integer, intent(in)     :: ngc                      !< Number of ghost cells
    integer, intent(in)     :: n_vars                   !< Number of variables
    real(dp), intent(in)    :: cc(1-ngc:nc+ngc, n_vars) !< Cell-centered values
    !> Reconstructed "left" values at every interface
    real(dp), intent(inout) :: u_l(1:nc+1, n_vars)
    !> Reconstructed "right" values at every interface
    real(dp), intent(inout) :: u_r(1:nc+1, n_vars)
    integer, intent(in)     :: limiter                  !< Which limiter to use
    real(dp)                :: slopes(0:nc+1, n_vars)

    associate (a=>cc(1:nc+2, :) - cc(0:nc+1, :), &
         b=>cc(0:nc+1, :) - cc(-1:nc, :))
      ! Compute limited slopes at the cell centers
      slopes = af_limiter_apply(a, b, limiter)

      ! Reconstruct values on the faces
      u_l(1:nc+1, :) = cc(0:nc, :) + 0.5_dp * slopes(0:nc, :) ! left

      if (af_limiter_symmetric(limiter)) then
         u_r(1:nc+1, :) = cc(1:nc+1, :) - 0.5_dp * slopes(1:nc+1, :) ! right
      else
         slopes = af_limiter_apply(b, a, limiter)
         u_r(1:nc+1, :) = cc(1:nc+1, :) - 0.5_dp * slopes(1:nc+1, :) ! right
      end if
    end associate
  end subroutine reconstruct_lr_1d

  !> Reconstruct upwind values at cell faces from cell-centered data.
  subroutine reconstruct_upwind_1d(nc, ngc, n_vars, cc, u_f, limiter, direction_positive)
    integer, intent(in)     :: nc                       !< Number of cells
    integer, intent(in)     :: ngc                      !< Number of ghost cells
    integer, intent(in)     :: n_vars                   !< Number of variables
    real(dp), intent(in)    :: cc(1-ngc:nc+ngc, n_vars) !< Cell-centered values
    !> Reconstructed values at every interface
    real(dp), intent(inout) :: u_f(1:nc+1, n_vars)
    integer, intent(in)     :: limiter                  !< Which limiter to use
    !> True means positive flow (to the "right"), false to the left
    logical, intent(in)     :: direction_positive(nc+1, n_vars)

    associate (a=>cc(1:nc+2, :) - cc(0:nc+1, :), &
         b=>cc(0:nc+1, :) - cc(-1:nc, :))
      where (direction_positive)
         u_f = cc(0:nc, :) + 0.5_dp * &
              af_limiter_apply(a(1:nc+1, :), b(1:nc+1, :), limiter)
      elsewhere
         u_f = cc(1:nc+1, :) - 0.5_dp * &
              af_limiter_apply(b(2:nc+2, :), a(2:nc+2, :), limiter)
      end where
    end associate
  end subroutine reconstruct_upwind_1d

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

  subroutine flux_update_densities(tree, dt, n_vars, i_cc, n_vars_flux, i_cc_flux, i_flux, &
       s_deriv, n_prev, s_prev, w_prev, s_out, add_source_box, n_dt, dt_lim, get_mask)
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt             !< Time step
    integer, intent(in)       :: n_vars         !< Number of cell-centered variables
    integer, intent(in)       :: i_cc(n_vars)   !< All cell-centered indices
    integer, intent(in)       :: n_vars_flux    !< Number of variables with fluxes
    integer, intent(in)       :: i_cc_flux(n_vars_flux) !< Cell-centered indices of flux variables
    integer, intent(in)       :: i_flux(n_vars_flux) !< Indices of fluxes
    integer, intent(in)       :: s_deriv        !< State to compute derivatives from
    integer, intent(in)       :: n_prev         !< Number of previous states
    integer, intent(in)       :: s_prev(n_prev) !< Previous states
    real(dp), intent(in)      :: w_prev(n_prev) !< Weights of previous states
    integer, intent(in)       :: s_out          !< Output time state
    !> Method to include source terms
    procedure(subr_source)    :: add_source_box
    !> Number of time steps restrictions
    integer, intent(in)       :: n_dt
    !> Maximal allowed time steps
    real(dp), intent(out)     :: dt_lim(n_dt)
    !> If present, only update where the mask is true
    procedure(subr_box_mask), optional :: get_mask
    integer                   :: lvl, n, id, IJK, nc, i_var, iv
    logical                   :: mask(DTIMES(tree%n_cell))
    real(dp)                  :: dt_dr(NDIM), my_dt(n_dt)
    real(dp)                  :: rfac(2, tree%n_cell)

    nc = tree%n_cell
    rfac = 0.0_dp ! Prevent warnings in 3D

    dt_lim = 1e100_dp

    !$omp parallel private(lvl, n, id, IJK, dt_dr, i_var, rfac, iv, mask, my_dt) &
    !$omp &reduction(min:dt_lim)
    my_dt = 1e100_dp

    do lvl = 1, tree%highest_lvl
       !$omp do
       do n = 1, size(tree%lvls(lvl)%leaves)
          id    = tree%lvls(lvl)%leaves(n)
          dt_dr = dt/tree%boxes(id)%dr

          associate(cc => tree%boxes(id)%cc, fc => tree%boxes(id)%fc)
            if (present(get_mask)) then
               call get_mask(tree%boxes(id), mask)
            else
               mask = .true.
            end if

            do i_var = 1, n_vars
               iv = i_cc(i_var)
               do KJI_DO(1, nc)
                  tree%boxes(id)%cc(IJK, iv+s_out) = sum(w_prev * &
                       tree%boxes(id)%cc(IJK, iv+s_prev))
               end do; CLOSE_DO
            end do

            call add_source_box(tree%boxes(id), dt, n_vars, i_cc, s_deriv, &
                 s_out, n_dt, my_dt, mask)
            dt_lim = min(dt_lim, my_dt)

#if NDIM == 1
            do KJI_DO(1, nc)
               if (mask(IJK)) then
                  cc(i, i_cc_flux+s_out) = cc(i, i_cc_flux+s_out) + &
                       dt_dr(1) * &
                       (fc(i, 1, i_flux) - fc(i+1, 1, i_flux))
               end if
            end do; CLOSE_DO
#elif NDIM == 2
            if (tree%coord_t == af_cyl) then
               call af_cyl_flux_factors(tree%boxes(id), rfac)
               do KJI_DO(1, nc)
                  if (mask(IJK)) then
                     cc(i, j, i_cc_flux+s_out) = cc(i, j, i_cc_flux+s_out) + &
                          dt_dr(1) * (&
                          rfac(1, i) * fc(i, j, 1, i_flux) - &
                          rfac(2, i) * fc(i+1, j, 1, i_flux)) &
                          + dt_dr(2) * &
                          (fc(i, j, 2, i_flux) - fc(i, j+1, 2, i_flux))
                  end if
               end do; CLOSE_DO
            else
               do KJI_DO(1, nc)
                  if (mask(IJK)) then
                     cc(i, j, i_cc_flux+s_out) = cc(i, j, i_cc_flux+s_out) + &
                          dt_dr(1) * &
                          (fc(i, j, 1, i_flux) - fc(i+1, j, 1, i_flux)) &
                          + dt_dr(2) * &
                          (fc(i, j, 2, i_flux) - fc(i, j+1, 2, i_flux))
                  end if
               end do; CLOSE_DO
            end if
#elif NDIM == 3
            do KJI_DO(1, nc)
               if (mask(IJK)) then
                  cc(i, j, k, i_cc_flux+s_out) = cc(i, j, k, i_cc_flux+s_out) + &
                       dt_dr(1) * &
                       (fc(i, j, k, 1, i_flux) - fc(i+1, j, k, 1, i_flux)) + &
                       dt_dr(2) * &
                       (fc(i, j, k, 2, i_flux) - fc(i, j+1, k, 2, i_flux)) + &
                       dt_dr(3) * &
                       (fc(i, j, k, 3, i_flux) - fc(i, j, k+1, 3, i_flux))
               end if
            end do; CLOSE_DO
#endif
          end associate
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine flux_update_densities

  !> Compute generic finite volume flux with a second order MUSCL scheme
  subroutine flux_generic_tree(tree, n_vars, i_cc, s_deriv, i_flux, dt_lim, &
       max_wavespeed, flux_from_primitives, flux_modify, line_modify, &
       to_primitive, to_conservative, limiter)
    use m_af_restrict
    use m_af_core
    type(af_t), intent(inout)     :: tree
    integer, intent(in)           :: n_vars         !< Number of variables
    integer, intent(in)           :: i_cc(n_vars)   !< Cell-centered variables
    integer, intent(in)           :: s_deriv        !< State to compute derivatives from
    integer, intent(in)           :: i_flux(n_vars) !< Flux variables
    !> Maximal time step, assuming a CFL number of 1.0
    real(dp), intent(out)         :: dt_lim
    !> Compute the maximum wave speed
    procedure(subr_max_wavespeed) :: max_wavespeed
    !> Compute the flux from primitive variables
    procedure(subr_flux)          :: flux_from_primitives
    !> Other flux contributions
    procedure(subr_flux_modify)   :: flux_modify
    !> Potentially modify line densities
    procedure(subr_line_modify)   :: line_modify
    !> Convert conservative variables to primitive ones
    procedure(subr_prim_cons)     :: to_primitive
    !> Convert primitive variables to conservative ones
    procedure(subr_prim_cons)     :: to_conservative
    !> Type of slope limiter to use for flux calculation
    integer, intent(in)           :: limiter
    real(dp)                      :: my_dt

    integer :: lvl, i

    ! Ensure ghost cells near refinement boundaries can be properly filled
    call af_restrict_ref_boundary(tree, i_cc+s_deriv)

    dt_lim = 1e100_dp
    my_dt = 1e100_dp

    !$omp parallel private(lvl, i, my_dt) reduction(min:dt_lim)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          call flux_generic_box(tree, tree%lvls(lvl)%leaves(i), tree%n_cell, &
               n_vars, i_cc, s_deriv, i_flux, my_dt, max_wavespeed, &
               flux_from_primitives, flux_modify, line_modify, &
               to_primitive, to_conservative, limiter)
          dt_lim = min(dt_lim, my_dt)
       end do
       !$omp end do
    end do
    !$omp end parallel

    ! Compute coarse fluxes from the fine ones at refinement boundaries
    call af_consistent_fluxes(tree, i_flux)

  end subroutine flux_generic_tree

  !> Compute generic finite volume flux with a second order MUSCL scheme
  subroutine flux_generic_box(tree, id, nc, n_vars, i_cc, s_deriv, i_flux, dt_lim, &
       max_wavespeed, flux_from_primitives, flux_modify, line_modify, &
       to_primitive, to_conservative, limiter)
    use m_af_types
    use m_af_ghostcell
    type(af_t), intent(inout)     :: tree
    integer, intent(in)           :: id             !< Id of box
    integer, intent(in)           :: nc             !< Number of cells
    integer, intent(in)           :: n_vars         !< Number of variables
    integer, intent(in)           :: i_cc(n_vars)   !< Cell-centered variables
    integer, intent(in)           :: s_deriv        !< State to compute derivatives from
    integer, intent(in)           :: i_flux(n_vars) !< Flux variables
    !> Maximal time step, assuming a CFL number of 1.0
    real(dp), intent(out)         :: dt_lim
    !> Compute the maximum wave speed
    procedure(subr_max_wavespeed) :: max_wavespeed
    !> Compute the flux from primitive variables on cell faces
    procedure(subr_flux)          :: flux_from_primitives
    !> Modify flux for other contributions
    procedure(subr_flux_modify)   :: flux_modify
    !> Potentially modify line densities
    procedure(subr_line_modify)   :: line_modify
    !> Convert conservative variables to primitive ones
    procedure(subr_prim_cons)     :: to_primitive
    !> Convert primitive variables to conservative ones
    procedure(subr_prim_cons)     :: to_conservative
    !> Type of slope limiter to use for flux calculation
    integer, intent(in)           :: limiter


    real(dp) :: cc(DTIMES(-1:nc+2), n_vars)
    real(dp) :: cc_line(-1:nc+2, n_vars)
    real(dp) :: flux(nc+1, n_vars)
    real(dp) :: u_l(nc+1, n_vars), u_r(nc+1, n_vars)
    real(dp) :: w_l(nc+1), w_r(nc+1)
    real(dp) :: flux_l(nc+1, n_vars), flux_r(nc+1, n_vars)
    real(dp) :: cfl_sum(DTIMES(nc)), cfl_sum_line(nc), inv_dr(NDIM)
    integer  :: flux_dim, line_ix(NDIM-1)
#if NDIM > 1
    integer  :: i
#endif
#if NDIM > 2
    integer  :: j
#endif

    inv_dr = 1/tree%boxes(id)%dr

    ! Get two layers of ghost cell data
    call af_gc2_box(tree, id, i_cc+s_deriv, cc)

    ! This will contain the sum of CFL-related conditions. For example, one can
    ! write dt * (vx/dx + vy/dy) < 1. This sum will then contain (vx/dx +
    ! vy/dy).
    cfl_sum = 0.0_dp

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
                cc_line = cc(i, j, :, :)
#endif
             end select

#if NDIM == 2
             line_ix = [i]
#elif NDIM == 3
             line_ix = [i, j]
#endif

             ! Optionally modify data, e.g. to take into account a boundary
             call line_modify(nc+4, n_vars, cc_line, flux_dim, &
                  tree%boxes(id), line_ix, s_deriv)

             ! Flux computation is based on primitive variables (this can be a
             ! dummy procedure)
             call to_primitive(nc+4, n_vars, cc_line)

             ! Reconstruct to cell faces, getting a left and right state at
             ! every face
             call reconstruct_lr_1d(nc, 2, n_vars, cc_line, u_l, u_r, limiter)

             ! Determine the maximal wave speed for the left and right state
             call max_wavespeed(nc+1, n_vars, flux_dim, u_l, w_l)
             call max_wavespeed(nc+1, n_vars, flux_dim, u_r, w_r)

             ! Compute the fluxes for the left and right state
             call flux_from_primitives(nc+1, n_vars, flux_dim, u_l, flux_l, &
                  tree%boxes(id), line_ix, s_deriv)
             call flux_from_primitives(nc+1, n_vars, flux_dim, u_r, flux_r, &
                  tree%boxes(id), line_ix, s_deriv)

             ! Convert back to conservative form
             call to_conservative(nc+1, n_vars, u_l)
             call to_conservative(nc+1, n_vars, u_r)

             ! Get maximum of left/right wave speed
             w_l = max(w_l, w_r)

             ! Get maximal wave speed for each cell center
             cfl_sum_line = max(w_l(1:nc), w_l(2:nc+1)) * inv_dr(NDIM)

             ! Combine left and right fluxes to obtain a single flux
             call flux_kurganovTadmor_1d(nc+1, n_vars, flux_l, flux_r, &
                  u_l, u_r, w_l, flux)

             ! Potentially add other flux components
             call flux_modify(nc+1, n_vars, flux_dim, flux, &
                  tree%boxes(id), line_ix, s_deriv)

             ! Store the computed fluxes
             select case (flux_dim)
#if NDIM == 1
             case (1)
                tree%boxes(id)%fc(:, flux_dim, i_flux) = flux
                cfl_sum = cfl_sum + cfl_sum_line
#elif NDIM == 2
             case (1)
                tree%boxes(id)%fc(:, i, flux_dim, i_flux) = flux
                cfl_sum(:, i) = cfl_sum(:, i) + cfl_sum_line
             case (2)
                tree%boxes(id)%fc(i, :, flux_dim, i_flux) = flux
                cfl_sum(i, :) = cfl_sum(i, :) + cfl_sum_line
#elif NDIM == 3
             case (1)
                tree%boxes(id)%fc(:, i, j, flux_dim, i_flux) = flux
                cfl_sum(:, i, j) = cfl_sum(:, i, j) + cfl_sum_line
             case (2)
                tree%boxes(id)%fc(i, :, j, flux_dim, i_flux) = flux
                cfl_sum(i, :, j) = cfl_sum(i, :, j) + cfl_sum_line
             case (3)
                tree%boxes(id)%fc(i, j, :, flux_dim, i_flux) = flux
                cfl_sum(i, j, :) = cfl_sum(i, j, :) + cfl_sum_line
#endif
             end select
#if NDIM > 1
          end do
#endif
#if NDIM > 2
       end do
#endif
    end do

    ! Determine maximal time step
    dt_lim = 1/max(maxval(cfl_sum), 1e-100_dp)

  end subroutine flux_generic_box

  !> Compute upwind fluxes
  subroutine flux_upwind_tree(tree, n_vars, i_cc, s_deriv, i_flux, n_dt, dt_lim, &
       flux_upwind, flux_direction, line_modify, limiter)
    use m_af_restrict
    use m_af_core
    type(af_t), intent(inout)   :: tree
    integer, intent(in)         :: n_vars         !< Number of variables
    integer, intent(in)         :: i_cc(n_vars)   !< Cell-centered variables
    integer, intent(in)         :: s_deriv        !< State to compute derivatives from
    integer, intent(in)         :: i_flux(n_vars) !< Flux variables
    !> Number of time steps restrictions (first is CFL)
    integer, intent(in)         :: n_dt
    !> Maximal allowed time steps, assuming a CFL number of 1.0
    real(dp), intent(out)       :: dt_lim(n_dt)
    !> Method to compute fluxes
    procedure(subr_flux_upwind) :: flux_upwind
    !> Method to get direction of flux (positive or negative)
    procedure(subr_flux_dir)    :: flux_direction
    !> Potentially modify line densities
    procedure(subr_line_modify) :: line_modify
    !> Type of slope limiter to use for flux calculation
    integer, intent(in)         :: limiter
    integer                     :: lvl, i
    real(dp)                    :: my_dt(n_dt)

    ! Ensure ghost cells near refinement boundaries can be properly filled
    call af_restrict_ref_boundary(tree, i_cc+s_deriv)

    dt_lim = 1e100_dp
    my_dt = 1e100_dp

    !$omp parallel private(lvl, i, my_dt) reduction(min:dt_lim)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          call flux_upwind_box(tree, tree%lvls(lvl)%leaves(i), tree%n_cell, &
               n_vars, i_cc, s_deriv, i_flux, n_dt, my_dt, flux_upwind, &
               flux_direction, line_modify, limiter)
          dt_lim = min(dt_lim, my_dt)
       end do
       !$omp end do
    end do
    !$omp end parallel

    ! Compute coarse fluxes from the fine ones at refinement boundaries
    call af_consistent_fluxes(tree, i_flux)

  end subroutine flux_upwind_tree

  !> Compute generic finite volume flux
  subroutine flux_upwind_box(tree, id, nc, n_vars, i_cc, s_deriv, i_flux, &
       n_dt, dt_lim, flux_upwind, flux_direction, line_modify, limiter)
    use m_af_types
    use m_af_ghostcell
    type(af_t), intent(inout)   :: tree
    integer, intent(in)         :: id             !< Id of box
    integer, intent(in)         :: nc             !< Number of cells
    integer, intent(in)         :: n_vars         !< Number of variables
    integer, intent(in)         :: i_cc(n_vars)   !< Cell-centered variables
    integer, intent(in)         :: s_deriv        !< State to compute derivatives from
    integer, intent(in)         :: i_flux(n_vars) !< Flux variables
    !> Number of time steps restrictions (first is CFL)
    integer, intent(in)         :: n_dt
    real(dp), intent(inout)     :: dt_lim(n_dt) !< Time step restrictions
    !> Method to compute fluxes
    procedure(subr_flux_upwind) :: flux_upwind
    !> Method to get direction of flux (positive or negative)
    procedure(subr_flux_dir)    :: flux_direction
    !> Potentially modify line densities
    procedure(subr_line_modify) :: line_modify
    !> Type of slope limiter to use for flux calculation
    integer, intent(in)         :: limiter

    real(dp) :: cc(DTIMES(-1:nc+2), n_vars)
    real(dp) :: cc_line(-1:nc+2, n_vars)
    real(dp) :: flux(nc+1, n_vars)
    real(dp) :: u_l(nc+1, n_vars)
    real(dp) :: cfl_sum(DTIMES(nc)), cfl_sum_line(nc)
    real(dp) :: other_dt(n_dt-1)
    logical  :: direction_positive(nc+1, n_vars)
    integer  :: flux_dim, line_ix(NDIM-1)
#if NDIM > 1
    integer  :: i
#endif
#if NDIM > 2
    integer  :: j
#endif

    ! Get two layers of ghost cell data
    call af_gc2_box(tree, id, i_cc+s_deriv, cc)

    ! This will contain the sum of CFL-related conditions. For example, one can
    ! write dt * (vx/dx + vy/dy) < 1. This sum will then contain (vx/dx +
    ! vy/dy).
    cfl_sum = 0.0_dp
    dt_lim(2:) = 1e100_dp
    other_dt = 1e100_dp

    associate (fc => tree%boxes(id)%fc)
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
                  cc_line = cc(i, j, :, :)
#endif
               end select

#if NDIM == 2
               line_ix = [i]
#elif NDIM == 3
               line_ix = [i, j]
#endif
               call flux_direction(tree%boxes(id), line_ix, s_deriv, &
                    n_vars, flux_dim, direction_positive)

               ! Optionally modify data, e.g. to take into account a boundary
               call line_modify(nc+4, n_vars, cc_line, flux_dim, &
                    tree%boxes(id), line_ix, s_deriv)

               ! Reconstruct to cell faces
               call reconstruct_upwind_1d(nc, 2, n_vars, cc_line, u_l, &
                    limiter, direction_positive)

               call flux_upwind(nc+1, n_vars, flux_dim, u_l, flux, cfl_sum_line, &
                    n_dt-1, other_dt, tree%boxes(id), line_ix, s_deriv)
               dt_lim(2:) = min(dt_lim(2:), other_dt)

               ! Store the computed fluxes
               select case (flux_dim)
#if NDIM == 1
               case (1)
                  fc(:, flux_dim, i_flux) = flux
                  cfl_sum = cfl_sum + cfl_sum_line
#elif NDIM == 2
               case (1)
                  fc(:, i, flux_dim, i_flux) = flux
                  cfl_sum(:, i) = cfl_sum(:, i) + cfl_sum_line
               case (2)
                  fc(i, :, flux_dim, i_flux) = flux
                  cfl_sum(i, :) = cfl_sum(i, :) + cfl_sum_line
#elif NDIM == 3
               case (1)
                  fc(:, i, j, flux_dim, i_flux) = flux
                  cfl_sum(:, i, j) = cfl_sum(:, i, j) + cfl_sum_line
               case (2)
                  fc(i, :, j, flux_dim, i_flux) = flux
                  cfl_sum(i, :, j) = cfl_sum(i, :, j) + cfl_sum_line
               case (3)
                  fc(i, j, :, flux_dim, i_flux) = flux
                  cfl_sum(i, j, :) = cfl_sum(i, j, :) + cfl_sum_line
#endif
               end select
#if NDIM > 1
            end do
#endif
#if NDIM > 2
         end do
#endif
      end do
    end associate

    ! Determine maximal CFL time step
    dt_lim(1) = 1/maxval(cfl_sum)

  end subroutine flux_upwind_box

  !> Extract cell-centered data along a line in a box, including a single layer
  !> of ghost cells. This is convenient to get extra variables in a flux
  !> computation.
  subroutine flux_get_line_cc(box, ivs, flux_dim, line_ix, cc_line)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: ivs(:)          !< Indices of the variables
    integer, intent(in)     :: flux_dim        !< Dimension of flux computation
    integer, intent(in)     :: line_ix(NDIM-1) !< Index of line
    real(dp), intent(inout) :: cc_line(box%n_cell+2, size(ivs))

    select case (flux_dim)
#if NDIM == 1
    case (1)
       cc_line = box%cc(:, ivs)
#elif NDIM == 2
    case (1)
       cc_line = box%cc(:, line_ix(1), ivs)
    case (2)
       cc_line = box%cc(line_ix(1), :, ivs)
#elif NDIM == 3
    case (1)
       cc_line = box%cc(:, line_ix(1), line_ix(2), ivs)
    case (2)
       cc_line = box%cc(line_ix(1), :, line_ix(2), ivs)
    case (3)
       cc_line = box%cc(line_ix(1), line_ix(2), :, ivs)
#endif
    end select
  end subroutine flux_get_line_cc

  !> Extract cell-centered data along a line in a box, including a single layer
  !> of ghost cells. This is convenient to get extra variables in a flux
  !> computation.
  subroutine flux_get_line_1cc(box, iv, flux_dim, line_ix, cc_line)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: iv              !< Index of the variable
    integer, intent(in)     :: flux_dim        !< Dimension of flux computation
    integer, intent(in)     :: line_ix(NDIM-1) !< Index of line
    real(dp), intent(inout) :: cc_line(box%n_cell+2)

    select case (flux_dim)
#if NDIM == 1
    case (1)
       cc_line = box%cc(:, iv)
#elif NDIM == 2
    case (1)
       cc_line = box%cc(:, line_ix(1), iv)
    case (2)
       cc_line = box%cc(line_ix(1), :, iv)
#elif NDIM == 3
    case (1)
       cc_line = box%cc(:, line_ix(1), line_ix(2), iv)
    case (2)
       cc_line = box%cc(line_ix(1), :, line_ix(2), iv)
    case (3)
       cc_line = box%cc(line_ix(1), line_ix(2), :, iv)
#endif
    end select
  end subroutine flux_get_line_1cc

  !> Extract face-centered data along a line in a box. This is convenient to get
  !> extra variables in a flux computation.
  subroutine flux_get_line_fc(box, ivs, flux_dim, line_ix, fc_line)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: ivs(:)          !< Indices of the variables
    integer, intent(in)     :: flux_dim        !< Dimension of flux computation
    integer, intent(in)     :: line_ix(NDIM-1) !< Index of line
    real(dp), intent(inout) :: fc_line(box%n_cell+1, size(ivs))

    select case (flux_dim)
#if NDIM == 1
    case (1)
       fc_line = box%fc(:, flux_dim, ivs)
#elif NDIM == 2
    case (1)
       fc_line = box%fc(:, line_ix(1), flux_dim, ivs)
    case (2)
       fc_line = box%fc(line_ix(1), :, flux_dim, ivs)
#elif NDIM == 3
    case (1)
       fc_line = box%fc(:, line_ix(1), line_ix(2), flux_dim, ivs)
    case (2)
       fc_line = box%fc(line_ix(1), :, line_ix(2), flux_dim, ivs)
    case (3)
       fc_line = box%fc(line_ix(1), line_ix(2), :, flux_dim, ivs)
#endif
    end select
  end subroutine flux_get_line_fc

  !> Extract face-centered data along a line in a box. This is convenient to get
  !> extra variables in a flux computation.
  subroutine flux_get_line_1fc(box, iv, flux_dim, line_ix, fc_line)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: iv              !< Index of the variable
    integer, intent(in)     :: flux_dim        !< Dimension of flux computation
    integer, intent(in)     :: line_ix(NDIM-1) !< Index of line
    real(dp), intent(inout) :: fc_line(box%n_cell+1)

    select case (flux_dim)
#if NDIM == 1
    case (1)
       fc_line = box%fc(:, flux_dim, iv)
#elif NDIM == 2
    case (1)
       fc_line = box%fc(:, line_ix(1), flux_dim, iv)
    case (2)
       fc_line = box%fc(line_ix(1), :, flux_dim, iv)
#elif NDIM == 3
    case (1)
       fc_line = box%fc(:, line_ix(1), line_ix(2), flux_dim, iv)
    case (2)
       fc_line = box%fc(line_ix(1), :, line_ix(2), flux_dim, iv)
    case (3)
       fc_line = box%fc(line_ix(1), line_ix(2), :, flux_dim, iv)
#endif
    end select
  end subroutine flux_get_line_1fc

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

  !> Dummy conversion between primitive and conservative variables
  subroutine flux_dummy_conversion(n_values, n_vars, u)
    integer, intent(in)     :: n_values, n_vars
    real(dp), intent(inout) :: u(n_values, n_vars)
  end subroutine flux_dummy_conversion

  subroutine flux_dummy_source(box, dt, n_vars, i_cc, s_deriv, s_out, n_dt, dt_lim, mask)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: dt
    integer, intent(in)        :: n_vars
    integer, intent(in)        :: i_cc(n_vars)
    integer, intent(in)        :: s_deriv
    integer, intent(in)        :: s_out
    integer, intent(in)        :: n_dt
    real(dp), intent(inout)    :: dt_lim(n_dt)
    logical, intent(in)        :: mask(DTIMES(box%n_cell))
  end subroutine flux_dummy_source

  subroutine flux_dummy_modify(nf, n_var, flux_dim, flux, box, line_ix, s_deriv)
    integer, intent(in)     :: nf              !< Number of cell faces
    integer, intent(in)     :: n_var           !< Number of variables
    integer, intent(in)     :: flux_dim        !< In which dimension fluxes are computed
    real(dp), intent(inout) :: flux(nf, n_var) !< Flux that will be modified
    type(box_t), intent(in) :: box             !< Current box
    integer, intent(in)     :: line_ix(NDIM-1) !< Index of line for dim /= flux_dim
    integer, intent(in)     :: s_deriv         !< State to compute derivatives from
  end subroutine flux_dummy_modify

  subroutine flux_dummy_line_modify(n_cc, n_var, cc_line, flux_dim, box, &
       line_ix, s_deriv)
    integer, intent(in)     :: n_cc                 !< Number of cell centers
    integer, intent(in)     :: n_var                !< Number of variables
    real(dp), intent(inout) :: cc_line(n_cc, n_var) !< Flux that will be modified
    integer, intent(in)     :: flux_dim             !< In which dimension fluxes are computed
    type(box_t), intent(in) :: box                  !< Current box
    integer, intent(in)     :: line_ix(NDIM-1)      !< Index of line for dim /= flux_dim
    integer, intent(in)     :: s_deriv              !< State to compute derivatives from
  end subroutine flux_dummy_line_modify

end module m_af_flux_schemes
