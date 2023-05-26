#include "../src/cpp_macros.h"
!> \example compressible_flow_wall.f90
!>
!> Example with compressible flow interacting with a solid wall
program compressible_flow_wall
  use m_af_all
  implicit none

  integer, parameter :: n_vars          = 2+NDIM
  integer, parameter :: box_size        = 8

  ! Indices defining the order of the flux variables
  integer, parameter :: i_rho       = 1
  integer, parameter :: i_mom(NDIM) = [2, 3]
  integer, parameter :: i_e         = 4
  integer            :: i_lsf

  real(dp), parameter :: k_boltzmann = 1.380649e-23_dp ! J/K
  real(dp), parameter :: dalton = 1.66053906660e-27_dp ! kg

  real(dp)           :: euler_gamma     = 1.4_dp
  real(dp)           :: cfl_number      = 0.5_dp
  integer            :: time_integrator = af_heuns_method
  real(dp)           :: rho_0, u_0, p_0, T_0, sound_speed, mach_number
  real(dp)           :: molecular_mass, domain_size(NDIM)
  real(dp)           :: wedge_angle
  integer            :: variables(n_vars)
  integer            :: cc_rho, cc_mom(2), cc_e
  integer            :: fluxes(n_vars)
  integer            :: grid_size(NDIM)
  real(dp)           :: u_init(1, n_vars)
  type(af_t)         :: tree
  real(dp)           :: dt, dt_lim, time, end_time, dt_output
  character(len=100) :: fname
  integer            :: n, output_cnt
  logical            :: use_amr

  character(len=10) :: var_names(n_vars)

  ! This test case with supersonic flow past a wedge is described in
  ! https://doi.org/10.3929/ethz-a-005575219 (master thesis)

  domain_size = [0.45_dp, 0.3_dp] ! m
  grid_size   = [12, 8] * box_size
  end_time       = 3e-3_dp       ! s
  euler_gamma    = 1.4_dp
  molecular_mass = 28.9_dp ! daltons
  p_0            = 101.35e3_dp        ! Pa
  T_0            = 288.9_dp           ! K
  wedge_angle    = 15.0_dp * acos(-1.0_dp) / 180.0_dp

  ! gas density in kg/m3
  rho_0          = p_0 * molecular_mass * dalton / (k_boltzmann * T_0)
  sound_speed    = sqrt(euler_gamma * p_0 / rho_0) ! m/s
  mach_number    = 2.5_dp
  u_0            = mach_number * sound_speed       ! m/s

  u_init(1, :) = [rho_0, u_0, 0.0_dp, p_0]
  call to_conservative(1, n_vars, u_init)

  var_names(i_rho) = "rho"
  var_names(i_mom) = ["mom_x", "mom_y"]
  var_names(i_e)   = "E"

  do n = 1, n_vars
     call af_add_cc_variable(tree, var_names(n), ix=variables(n), n_copies=2)
     call af_add_fc_variable(tree, "flux", ix=fluxes(n))
  end do

  call af_set_cc_methods(tree, variables(i_rho), bc_gas_rho)
  call af_set_cc_methods(tree, variables(i_mom(1)), bc_gas_mom_x)
  call af_set_cc_methods(tree, variables(i_mom(2)), bc_gas_mom_y)
  call af_set_cc_methods(tree, variables(i_e), bc_gas_e)

  call af_add_cc_variable(tree, "lsf", ix=i_lsf)
  call af_set_cc_methods(tree, i_lsf, funcval=set_lsf_box)

  cc_rho = variables(i_rho)
  cc_mom = variables(i_mom)
  cc_e   = variables(i_e)

  use_amr = .false.

  call af_init(tree, box_size, domain_size, grid_size)
  call af_loop_box(tree, set_init_conds)
  call af_gc_tree(tree, variables)

  time             = 0.0_dp
  dt               = 0.0_dp ! Start from zero time step
  output_cnt       = 0
  dt_output        = end_time / 100

  do
     if (output_cnt * dt_output <= time) then
        write(fname, "(A,I0)") "output/compressible_flow_wall_" &
             // DIMNAME // "_", output_cnt
        call af_write_silo(tree, trim(fname), output_cnt, time, &
             add_vars = write_primitives, add_names=["xVel","yVel","pres"])
        output_cnt = output_cnt + 1
     end if

     call af_advance(tree, dt, dt_lim, time, variables, &
          time_integrator, forward_euler)
     dt = dt_lim
     if (time > end_time) exit
  end do

  call af_destroy(tree)

contains

  subroutine forward_euler(tree, dt, dt_stiff, dt_lim, time, s_deriv, n_prev, &
       s_prev, w_prev, s_out, i_step, n_steps)
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt
    real(dp), intent(in)      :: dt_stiff       !< Time step for stiff terms
    real(dp), intent(inout)   :: dt_lim
    real(dp), intent(in)      :: time
    integer, intent(in)       :: s_deriv
    integer, intent(in)       :: n_prev         !< Number of previous states
    integer, intent(in)       :: s_prev(n_prev) !< Previous states
    real(dp), intent(in)      :: w_prev(n_prev) !< Weights of previous states
    integer, intent(in)       :: s_out
    integer, intent(in)       :: i_step, n_steps
    real(dp)                  :: wmax(NDIM)

    call flux_generic_tree(tree, n_vars, variables, s_deriv, fluxes, wmax, &
         max_wavespeed, get_fluxes, flux_dummy_modify, line_modify, &
         to_primitive, to_conservative, af_limiter_vanleer_t)
    call flux_update_densities(tree, dt, n_vars, variables, fluxes, &
         s_deriv, n_prev, s_prev, w_prev, s_out, flux_dummy_source, i_lsf)

    ! Compute new time step
    dt_lim = cfl_number / sum(wmax/af_lvl_dr(tree, tree%highest_lvl))
  end subroutine forward_euler

  subroutine set_init_conds(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK, nc

    nc = box%n_cell
    do KJI_DO(0, nc+1)
       box%cc(IJK, variables) = u_init(1, :)
    end do; CLOSE_DO
  end subroutine set_init_conds

  !> Convert variables to primitive form
  subroutine to_primitive(n_values, n_vars, u)
    integer, intent(in)     :: n_values, n_vars
    real(dp), intent(inout) :: u(n_values, n_vars)

    u(:, i_mom(1)) = u(:, i_mom(1))/u(:, i_rho)
    u(:, i_mom(2)) = u(:, i_mom(2))/u(:, i_rho)
    u(:, i_e) = (euler_gamma-1.0_dp) * (u(:, i_e) - &
         0.5_dp*u(:, i_rho)* sum(u(:, i_mom(:))**2, dim=2))
  end subroutine to_primitive

  !> Convert variables to conservative form
  subroutine to_conservative(n_values, n_vars, u)
    integer, intent(in)     :: n_values, n_vars
    real(dp), intent(inout) :: u(n_values, n_vars)
    real(dp)                :: kin_en(n_values)
    real(dp)                :: inv_fac
    integer                 :: i

    ! Compute kinetic energy (0.5 * rho * velocity^2)
    kin_en = 0.5_dp * u(:, i_rho) * sum(u(:, i_mom(:))**2, dim=2)

    ! Compute energy from pressure and kinetic energy
    inv_fac = 1/(euler_gamma - 1.0_dp)
    u(:, i_e) = u(:, i_e) * inv_fac + kin_en

    ! Compute momentum from density and velocity components
    do i = 1, NDIM
       u(:, i_mom(i)) = u(:, i_rho) * u(:, i_mom(i))
    end do
  end subroutine to_conservative

  !> Get maximal wave speed
  subroutine max_wavespeed(n_values, n_var, flux_dim, u, w)
    integer, intent(in)   :: n_values !< Number of cell faces
    integer, intent(in)   :: n_var    !< Number of variables
    integer, intent(in)   :: flux_dim !< In which dimension fluxes are computed
    real(dp), intent(in)  :: u(n_values, n_var) !< Primitive variables
    real(dp), intent(out) :: w(n_values) !< Maximum speed
    real(dp)              :: sound_speeds(n_values)

    sound_speeds = sqrt(euler_gamma * u(:, i_e) / u(:, i_rho))
    w = sound_speeds + abs(u(:, i_mom(flux_dim)))
  end subroutine max_wavespeed

  !> Compute fluxes from variables in primitive form
  subroutine get_fluxes(n_values, n_var, flux_dim, u, flux, box, line_ix, s_deriv)
    integer, intent(in)     :: n_values !< Number of cell faces
    integer, intent(in)     :: n_var    !< Number of variables
    integer, intent(in)     :: flux_dim !< In which dimension fluxes are computed
    real(dp), intent(in)    :: u(n_values, n_var)
    real(dp), intent(out)   :: flux(n_values, n_var)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: line_ix(NDIM-1)
    integer, intent(in)     :: s_deriv
    real(dp)                :: E(n_values), inv_fac
    integer                 :: i

    ! Compute left and right flux for conservative variables from the primitive
    ! reconstructed values.

    ! Density flux
    flux(:, i_rho) = u(:, i_rho) * u(:, i_mom(flux_dim))

    ! Momentum flux
    do i = 1, NDIM
       flux(:,  i_mom(i)) = u(:, i_rho) * &
            u(:, i_mom(i)) * u(:, i_mom(flux_dim))
    end do

    ! Add pressure term
    flux(:, i_mom(flux_dim)) = flux(:, i_mom(flux_dim)) + u(:, i_e)

    ! Compute energy
    inv_fac = 1/(euler_gamma-1.0_dp)
    E = u(:, i_e) * inv_fac + 0.5_dp * u(:, i_rho) * &
         sum(u(:, i_mom(:))**2, dim=2)

    ! Energy flux
    flux(:, i_e) = u(:, i_mom(flux_dim)) * (E + u(:, i_e))

  end subroutine get_fluxes

  !> Modify line data to approximate a slip boundary condition. However, in case
  !> of a staircase pattern, this becomes more of a no-slip boundary condition.
  subroutine line_modify(n_cc, n_var, cc_line, flux_dim, box, line_ix, s_deriv)
    integer, intent(in)     :: n_cc                 !< Number of cell centers
    integer, intent(in)     :: n_var                !< Number of variables
    real(dp), intent(inout) :: cc_line(-1:n_cc-2, n_var) !< Line values to modify
    integer, intent(in)     :: flux_dim             !< In which dimension fluxes are computed
    type(box_t), intent(in) :: box                  !< Current box
    integer, intent(in)     :: line_ix(NDIM-1)      !< Index of line for dim /= flux_dim
    integer, intent(in)     :: s_deriv              !< State to compute derivatives from

    real(dp)                :: lsf(0:box%n_cell+2, 1)
    integer                 :: i

    ! Get level set function along the line of the flux computation
    call flux_get_line_cc(box, [i_lsf], flux_dim, line_ix, lsf)

    if (all(lsf > 0)) return    ! no boundary

    do i = 0, box%n_cell
       if (lsf(i, 1) * lsf(i+1, 1) <= 0) then
          ! There is an interface
          if (lsf(i, 1) > 0) then
             cc_line(i+1, :) = cc_line(i, :)
             cc_line(i+2, :) = cc_line(i-1, :)
             cc_line(i+1, i_mom(flux_dim)) = -cc_line(i, i_mom(flux_dim))
             cc_line(i+2, i_mom(flux_dim)) = -cc_line(i-1, i_mom(flux_dim))
          else
             cc_line(i, :) = cc_line(i+1, :)
             cc_line(i-1, :) = cc_line(i+2, :)
             cc_line(i, i_mom(flux_dim)) = -cc_line(i+1, i_mom(flux_dim))
             cc_line(i-1, i_mom(flux_dim)) = -cc_line(i+2, i_mom(flux_dim))
          end if
       end if
    end do

  end subroutine line_modify

  !> Write primitive variables to output
  subroutine write_primitives(box, new_vars, n_var)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: n_var
    real(dp)                :: new_vars(DTIMES(0:box%n_cell+1),n_var)

    ! X Velocity
    new_vars(DTIMES(:), 1) = box%cc(DTIMES(:), cc_mom(1))/box%cc(DTIMES(:), cc_rho)
    ! Y Velocity
    new_vars(DTIMES(:), 2) = box%cc(DTIMES(:), cc_mom(2))/box%cc(DTIMES(:), cc_rho)
    ! Pressure
    new_vars(DTIMES(:), 3) = (euler_gamma-1.0_dp)*(box%cc(DTIMES(:), cc_e) - &
         sum(box%cc(DTIMES(:), cc_mom(:))**2, dim=NDIM+1) &
         /(2.0_dp*box%cc(DTIMES(:), cc_rho)))
  end subroutine write_primitives

  !> Boundary conditions for gas density
  subroutine bc_gas_rho(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    if (nb == af_neighb_lowx) then
       bc_type = af_bc_dirichlet
       bc_val = u_init(1, i_rho)
    else
       call af_bc_neumann_zero(box, nb, iv, coords, bc_val, bc_type)
    end if
  end subroutine bc_gas_rho

  !> Boundary conditions for x-momentum
  subroutine bc_gas_mom_x(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    if (nb == af_neighb_lowx) then
       bc_type = af_bc_dirichlet
       bc_val = u_init(1, i_mom(1))
    else
       call af_bc_neumann_zero(box, nb, iv, coords, bc_val, bc_type)
    end if
  end subroutine bc_gas_mom_x

  !> Boundary conditions for y-momentum
  subroutine bc_gas_mom_y(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    if (nb == af_neighb_lowx) then
       bc_type = af_bc_dirichlet
       bc_val = u_init(1, i_mom(2))
    else if (nb == af_neighb_lowy) then
       bc_type = af_bc_dirichlet
       bc_val = 0
    else if (nb == af_neighb_highy) then
       bc_type = af_bc_dirichlet
       bc_val = 0
    else
       call af_bc_neumann_zero(box, nb, iv, coords, bc_val, bc_type)
    end if
  end subroutine bc_gas_mom_y

  !> Boundary conditions for energy
  subroutine bc_gas_e(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    if (nb == af_neighb_lowx) then
       bc_type = af_bc_dirichlet
       bc_val = u_init(1, i_e)
    else
       call af_bc_neumann_zero(box, nb, iv, coords, bc_val, bc_type)
    end if
  end subroutine bc_gas_e

  subroutine set_lsf_box(box, iv)
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: iv
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM), norm_dr

    nc = box%n_cell
    norm_dr = norm2(box%dr)

    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, iv) = get_lsf(rr)
    end do; CLOSE_DO
  end subroutine set_lsf_box

  !> Level set function
  real(dp) function get_lsf(rr)
    real(dp), intent(in) :: rr(NDIM)

    ! Distance from a point to a line
    ! https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
    get_lsf = cos(wedge_angle) * (rr(2) - 0.0_dp) - &
         sin(wedge_angle) * (rr(1) - 0.15_dp)

  end function get_lsf

end program compressible_flow_wall
