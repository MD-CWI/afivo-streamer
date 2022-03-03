#include "../src/cpp_macros.h"

program KT_euler
  use m_af_all
  implicit none

  integer, parameter :: n_vars = 2+NDIM
  integer, parameter :: ncells = 8
  integer, parameter :: coord_type = af_xyz
  real(dp), parameter :: euler_gamma = 1.4_dp
  integer, parameter :: time_integrator = af_heuns_method

  ! Indices defining the order of the flux variables
  integer, parameter :: i_rho = 1
  integer, parameter :: i_mom(NDIM) = [2, 3]
  integer, parameter :: i_e = 4

  integer            :: variables(n_vars)
  integer            :: cc_rho, cc_mom(2), cc_e
  integer            :: fluxes(n_vars)
  real(dp)           :: l_max(NDIM), l_min(NDIM)
  integer            :: grid(NDIM)
  logical            :: periodicBC(NDIM)
  real(dp)           :: u0(4, n_vars)
  type(af_t)         :: tree
  real(dp)           :: dt, dt_lim, time, end_time, dt_output
  real(dp)           :: dt_amr, time_prev_refine
  character(len=100) :: fname
  character(len=20)  :: carg
  integer            :: n, output_cnt
  real(dp)           :: rho_max
  logical            :: benchmark
  logical            :: use_amr

  ! AMR stuff
  type(ref_info_t) :: refine_info
  integer          :: refine_steps

  character(len=10) :: var_names(n_vars)

  var_names(i_rho) = "rho"
  var_names(i_mom) = ["mom_x", "mom_y"]
  var_names(i_e)   = "E"

  print *, "Running Euler 2D with KT scheme"
  print *, "Number of threads", af_get_max_threads()

  do n = 1, n_vars
     call af_add_cc_variable(tree, var_names(n), ix=variables(n), n_copies=2)
     call af_add_fc_variable(tree, "flux", ix=fluxes(n))
     call af_set_cc_methods(tree, variables(n), af_bc_neumann_zero)
  end do

  cc_rho = variables(i_rho)
  cc_mom = variables(i_mom)
  cc_e   = variables(i_e)

  carg = "sod"
  if (command_argument_count() > 0) then
     call get_command_argument(1, carg)
  end if
  print *, "Initial condition type: ", trim(carg)

  select case (carg)
  case ("first")
     u0(:, i_e)      = [1.0_dp, 0.4_dp, 0.0439_dp, 0.15_dp]
     u0(:, i_rho)    = [1.0_dp, 0.5197_dp, 0.1072_dp, 0.2579_dp]
     u0(:, i_mom(1)) = [0.0_dp, -0.7259_dp, -0.7259_dp, 0.0_dp]
     u0(:, i_mom(2)) = [0.0_dp, 0.0_dp, -1.4045_dp, -1.4045_dp]
  case ("sixth")
     u0(:, i_e)      = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
     u0(:, i_rho)    = [1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]
     u0(:, i_mom(1)) = [0.75_dp, 0.75_dp, -0.75_dp, -0.75_dp]
     u0(:, i_mom(2)) = [-0.5_dp, 0.5_dp, 0.5_dp, -0.5_dp]
  case ("sod")
     ! 1D Sod shock test case
     u0(:, i_rho)    = [0.125_dp, 1.0_dp, 1.0_dp, 0.125_dp]
     u0(:, i_e)      = [0.1_dp, 1.0_dp, 1.0_dp, 0.1_dp]
     u0(:, i_mom(1)) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
     u0(:, i_mom(2)) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
  case default
     error stop "Unknown initial condition"
  end select

  call to_conservative(4, n_vars, u0)

  use_amr = .false.
  if (command_argument_count() > 1) then
     call get_command_argument(2, carg)
     read(carg, *) use_amr
  end if
  print *, "Use AMR: ", use_amr

  benchmark = .false.
  if (command_argument_count() > 2) then
     call get_command_argument(3, carg)
     read(carg, *) benchmark
  end if
  print *, "Benchmark: ", benchmark

  grid(:)       = 50*ncells
  l_max(:)      = 1.0_dp
  l_min(:)      = 0.0_dp
  periodicBC(:) = .false.

  call af_init(tree, ncells, l_max-l_min, grid, &
       periodic=periodicBC, r_min=l_min, &
       coord=coord_type)

  if (use_amr) then
     do
        refine_steps = refine_steps + 1
        !Settng init conds for each refinement is needed as we use that data as
        !refinement criterion
        call af_loop_box(tree, set_init_conds)
        call af_gc_tree(tree, variables)
        call af_adjust_refinement(tree, ref_rout, refine_info, 1)
        if (refine_info%n_add == 0) exit
     end do
     call af_restrict_tree(tree, variables)
     call af_gc_tree(tree, variables)
  else
     call af_loop_box(tree, set_init_conds)
     call af_gc_tree(tree, variables)
  end if

  time             = 0.0_dp
  time_prev_refine = time
  end_time         = 0.2_dp
  dt               = 0.0_dp ! Start from zero time step
  output_cnt       = 0
  dt_output        = end_time / 20
  dt_amr           = end_time / 100

  do
     if (.not. benchmark .and. output_cnt * dt_output <= time) then
        write(fname, "(A,I0)") "output/euler_" // DIMNAME // "_", output_cnt
        call af_write_silo(tree, trim(fname), output_cnt, time, &
             add_vars = write_primitives, add_names=["xVel","yVel","pres"])
        output_cnt = output_cnt + 1
     end if

     call af_advance(tree, dt, dt_lim, time, variables, &
          time_integrator, forward_euler)

     dt = 0.8_dp * dt_lim

     if (use_amr .and. time > time_prev_refine + dt_amr) then
        call af_gc_tree(tree, variables)
        call af_adjust_refinement(tree, ref_rout, refine_info, 1)
        time_prev_refine = time
     end if

     call af_tree_maxabs_cc(tree, variables(i_rho), rho_max)
     if (rho_max > 10.0_dp) &
          error stop "solution diverging!"

     if (time > end_time) exit
  end do

  call af_destroy(tree)

contains

  !> [forward_euler_gasd]
  subroutine forward_euler(tree, dt, dt_lim, time, s_deriv, s_prev, s_out, &
       i_step, n_steps)
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt
    real(dp), intent(inout)   :: dt_lim
    real(dp), intent(in)      :: time
    integer, intent(in)       :: s_deriv
    integer, intent(in)       :: s_prev
    integer, intent(in)       :: s_out
    integer, intent(in)       :: i_step, n_steps
    real(dp)                  :: wmax(NDIM)

    call flux_generic_tree(tree, n_vars, variables+s_deriv, fluxes, wmax, &
         max_wavespeed, get_fluxes, to_primitive, to_conservative)
    call flux_update_densities(tree, dt, n_vars, variables, fluxes, &
         s_deriv, s_prev, s_out, flux_dummy_source)

    ! Compute new time step
    dt_lim = 1.0_dp / sum(wmax/af_lvl_dr(tree, tree%highest_lvl))
  end subroutine forward_euler
  !> [forward_euler_gasd]

  subroutine set_init_conds(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM)

    nc = box%n_cell

    do KJI_DO(0, nc+1)
       rr = af_r_cc(box, [IJK])
       if (rr(1) > 0.5_dp .and. rr(2) > 0.5_dp) then
          box%cc(IJK, variables) = u0(1, :)
       elseif (rr(1) <= 0.5_dp .and. rr(2) >= 0.5_dp) then
          box%cc(IJK, variables) = u0(2, :)
       elseif (rr(1) <= 0.5_dp .and. rr(2) <= 0.5_dp) then
          box%cc(IJK, variables) = u0(3, :)
       else
          box%cc(IJK, variables) = u0(4, :)
       end if
    end do; CLOSE_DO
  end subroutine set_init_conds

  subroutine to_primitive(n_values, n_vars, u)
    integer, intent(in)     :: n_values, n_vars
    real(dp), intent(inout) :: u(n_values, n_vars)

    u(:, i_mom(1)) = u(:, i_mom(1))/u(:, i_rho)
    u(:, i_mom(2)) = u(:, i_mom(2))/u(:, i_rho)
    u(:, i_e) = (euler_gamma-1.0_dp) * (u(:, i_e) - &
         0.5_dp*u(:, i_rho)* sum(u(:, i_mom(:))**2, dim=2))
  end subroutine to_primitive

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

  subroutine get_fluxes(n_values, n_var, flux_dim, u, flux, box, line_ix)
    integer, intent(in)     :: n_values !< Number of cell faces
    integer, intent(in)     :: n_var    !< Number of variables
    integer, intent(in)     :: flux_dim !< In which dimension fluxes are computed
    real(dp), intent(in)    :: u(n_values, n_var)
    real(dp), intent(out)   :: flux(n_values, n_var)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: line_ix(NDIM-1)
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

  subroutine ref_rout( box, cell_flags )
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))
    real(dp)                :: diff, tol
    integer                 :: IJK, nc

    nc  = box%n_cell
    tol = 1.0e-6_dp
    do KJI_DO(1,nc)
       diff =  box%dr(1)**2*abs(box%cc(i+1,j,i_rho)+box%cc(i-1,j,i_rho) &
            -2*box%cc(i,j,i_rho)) + &
            box%dr(2)**2*abs(box%cc(i,j+1,i_rho)+box%cc(i,j-1,i_rho) &
            -2*box%cc(i,j,i_rho))
       if (diff > tol .and. box%lvl < 3) then
          cell_flags(IJK) = af_do_ref
       else if (diff < 0.1_dp*tol) then
          cell_flags(IJK) = af_rm_ref
       else
          cell_flags(IJK) = af_keep_ref
       end if
    end do; CLOSE_DO
  end subroutine ref_rout

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

end program KT_euler
