#include "../src/cpp_macros.h"
!> \example incompressible_flow.f90
!>
!> An incompressible flow example
program incompressible_flow
  use m_af_all

  implicit none

  integer, parameter          :: box_size        = 8
  integer, parameter          :: i_mom(NDIM)     = [1, 2]
  integer, parameter          :: n_vars          = 2
  real(dp), parameter         :: domain_len      = 2.0_dp
  character(len=2), parameter :: mom_names(NDIM) = ["ux", "uy"]

  ! x-velocity at top of domain
  real(dp) :: ux_top = 1.0_dp

  ! Viscosity
  real(dp) :: nu = 0.1_dp !1e-3_dp

  ! Source term x-velocity
  real(dp) :: source_ux = 0.0_dp

  type(mg_t) :: mg

  integer, parameter :: integrator = af_midpoint_method
  integer, parameter :: n_states = af_advance_num_steps(integrator)

  integer            :: global_iv_velocity(NDIM)
  integer            :: i, it
  integer            :: variables(NDIM)
  integer            :: fluxes(NDIM)
  type(af_t)         :: tree
  integer            :: output_cnt
  real(dp)           :: dt, dt_lim, dt_diff, time
  real(dp)           :: end_time
  real(dp)           :: dt_output
  character(len=100) :: fname, test_case

  print *, "Running incompressible_flow_" // DIMNAME // ""
  print *, "Number of threads", af_get_max_threads()

  do i = 1, NDIM
     call af_add_cc_variable(tree, mom_names(i), ix=variables(i), &
          n_copies=n_states)
     call af_add_fc_variable(tree, "flux", ix=fluxes(i))
  end do

  call af_add_cc_variable(tree, "phi", ix=mg%i_phi)
  call af_add_cc_variable(tree, "rhs", ix=mg%i_rhs)
  call af_add_cc_variable(tree, "tmp", ix=mg%i_tmp)
  mg%sides_bc => af_bc_neumann_zero
  mg%subtract_mean = .true.

  test_case = "channel_flow"

  if (command_argument_count() > 0) then
     call get_command_argument(1, test_case)
  end if

  select case (test_case)
  case ("cavity_flow")
     ux_top = 1.0_dp
     nu = 0.1_dp
     source_ux = 0.0_dp

     call af_set_cc_methods(tree, variables(1), bc_ux)
     call af_set_cc_methods(tree, variables(2), af_bc_dirichlet_zero)

     call af_init(tree, box_size, [DTIMES(domain_len)], [DTIMES(box_size)])
  case ("channel_flow")
     nu = 0.1_dp
     source_ux = 1.0_dp

     call af_set_cc_methods(tree, variables(1), af_bc_dirichlet_zero)
     call af_set_cc_methods(tree, variables(2), af_bc_dirichlet_zero)
     call af_init(tree, box_size, [DTIMES(domain_len)], [DTIMES(box_size)], &
          periodic=[.true., .false.])
  case default
     error stop "Invalid test case, choices: cavity_flow, channel_flow"
  end select

  call mg_init(tree, mg)

  it         = 0
  output_cnt = 0
  time       = 0
  dt_output  = 1.0_dp
  end_time   = 30.0_dp

  call af_refine_up_to_lvl(tree, 3)
  call af_loop_box(tree, set_initial_condition)

  call af_gc_tree(tree, i_mom)

  dt_diff = 0.5_dp * 0.25_dp * af_min_dr(tree)**2 / nu
  dt      = min(dt_diff, 1e-6_dp)

  ! Starting simulation
  do while (time < end_time)
     if (output_cnt * dt_output <= time) then
        output_cnt = output_cnt + 1
        write(fname, "(A,I0)") "output/incompressible_flow_" // DIMNAME // "_", output_cnt
        call af_gc_tree(tree, variables)
        call af_write_silo(tree, trim(fname), output_cnt, time)
        print *, it, dt, dt_lim
     end if

     call af_advance(tree, dt, dt_lim, time, variables, &
          af_midpoint_method, forward_euler)
     dt = 0.9_dp * min(dt_lim, dt_diff, 2.0_dp * dt)
     it = it + 1
     ! print *, it, dt, dt_lim
  end do

contains

  subroutine set_initial_condition(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK, nc

    nc = box%n_cell
    do KJI_DO(0,nc+1)
       box%cc(IJK, variables) = 0.0_dp
    end do; CLOSE_DO
  end subroutine set_initial_condition

  subroutine forward_euler(tree, dt, dt_lim, time, s_deriv, n_prev, s_prev, &
       w_prev, s_out, i_step, n_steps)
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt
    real(dp), intent(in)      :: time
    real(dp), intent(inout)   :: dt_lim
    integer, intent(in)       :: s_deriv
    integer, intent(in)       :: n_prev         !< Number of previous states
    integer, intent(in)       :: s_prev(n_prev) !< Previous states
    real(dp), intent(in)      :: w_prev(n_prev) !< Weights of previous states
    integer, intent(in)       :: s_out
    integer, intent(in)       :: i_step, n_steps
    real(dp)                  :: wmax(NDIM), tmp

    call flux_generic_tree(tree, n_vars, variables+s_deriv, fluxes, wmax, &
         max_wavespeed, get_flux, flux_dummy_conversion, flux_dummy_conversion)

    call flux_update_densities(tree, dt, n_vars, variables, fluxes, &
         s_deriv, n_prev, s_prev, w_prev, s_out, source_term)

    call remove_velocity_divergence(tree, variables+s_out)

    ! Compute maximal time step
    tmp = sum(wmax/af_lvl_dr(tree, tree%highest_lvl))
    dt_lim = 1.0_dp / max(tmp, 1.0e-10_dp)
  end subroutine forward_euler

  subroutine source_term(box, dt, n_vars, i_cc, s_deriv, s_out)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: dt
    integer, intent(in)        :: n_vars
    integer, intent(in)        :: i_cc(n_vars)
    integer, intent(in)        :: s_deriv
    integer, intent(in)        :: s_out

    integer :: nc
    nc = box%n_cell

    box%cc(DTIMES(1:nc), i_cc(i_mom(1)) + s_out) = &
         box%cc(DTIMES(1:nc), i_cc(i_mom(1)) + s_out) + dt * source_ux
  end subroutine source_term

  subroutine max_wavespeed(n_values, n_var, flux_dim, u, w)
    integer, intent(in)   :: n_values !< Number of cell faces
    integer, intent(in)   :: n_var    !< Number of variables
    integer, intent(in)   :: flux_dim !< In which dimension fluxes are computed
    real(dp), intent(in)  :: u(n_values, n_var) !< Primitive variables
    real(dp), intent(out) :: w(n_values) !< Maximum speed

    w = abs(u(:, i_mom(flux_dim)))
  end subroutine max_wavespeed

  subroutine get_flux(n_values, n_var, flux_dim, u, flux, box, line_ix, u_diff)
    integer, intent(in)     :: n_values !< Number of cell faces
    integer, intent(in)     :: n_var    !< Number of variables
    integer, intent(in)     :: flux_dim !< In which dimension fluxes are computed
    real(dp), intent(in)    :: u(n_values, n_var)
    real(dp), intent(out)   :: flux(n_values, n_var)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: line_ix(NDIM-1)
    real(dp), intent(in)    :: u_diff(n_values, n_var)
    real(dp)                :: inv_dr(NDIM)

    inv_dr = 1/box%dr

    do i = 1, NDIM
       ! Momentum flux
       flux(:,  i_mom(i)) = u(:, i_mom(i)) * u(:, i_mom(flux_dim))

       ! Viscosity
       flux(:,  i_mom(i)) = flux(:,  i_mom(i)) - &
            nu * u_diff(:, i_mom(i)) * inv_dr(i)
    end do

  end subroutine get_flux

  subroutine remove_velocity_divergence(tree, iv_velocity)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: iv_velocity(NDIM)
    logical                   :: use_4th_order = .false.
    real(dp)                  :: tmp
    logical, parameter        :: print_divergence = .false.

    global_iv_velocity = iv_velocity
    call af_restrict_ref_boundary(tree, iv_velocity)

    if (use_4th_order) then
       call af_loop_tree(tree, box_set_rhs_4th_order, .true.)
    else
       call af_gc_tree(tree, iv_velocity, corners=.false., leaves_only=.true.)
       call af_loop_tree(tree, box_set_rhs_2nd_order, .true.)
    end if

    if (print_divergence) then
       call af_tree_sum_cc(tree, mg%i_rhs, tmp)
       print *, "BEFORE", tmp
    end if

    ! call mg_fas_fmg(tree, mg, .false., .true.)
    ! call mg_fas_fmg(tree, mg, .true., .true.)
    call mg_fas_vcycle(tree, mg, .true.)
    call af_loop_box(tree, box_correct_velocity, .true.)

    if (print_divergence) then
       if (use_4th_order) then
          call af_loop_tree(tree, box_set_rhs_4th_order, .true.)
       else
          call af_gc_tree(tree, iv_velocity, corners=.false., leaves_only=.true.)
          call af_loop_tree(tree, box_set_rhs_2nd_order, .true.)
       end if
       call af_tree_sum_cc(tree, mg%i_rhs, tmp)
       print *, "AFTER", tmp
    end if
  end subroutine remove_velocity_divergence

  subroutine box_set_rhs_2nd_order(tree, id)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: id
    integer                   :: IJK, nc, iv(NDIM)
    real(dp)                  :: tmp(NDIM)

    associate (box => tree%boxes(id))
      tmp = 0.5_dp / box%dr
      iv = global_iv_velocity
      nc = box%n_cell

      do KJI_DO(1,nc)
         box%cc(IJK, mg%i_rhs) = &
              tmp(1) * (box%cc(i+1, j, iv(1)) - box%cc(i-1, j, iv(1))) + &
              tmp(2) * (box%cc(i, j+1, iv(2)) - box%cc(i, j-1, iv(2)))
      end do; CLOSE_DO
    end associate
  end subroutine box_set_rhs_2nd_order

  subroutine box_set_rhs_4th_order(tree, id)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: id
    integer                   :: IJK, nc, iv(NDIM)
    real(dp)                  :: tmp(NDIM)
    real(dp)                  :: cc(DTIMES(-1:tree%n_cell+2), n_vars)

    associate (box => tree%boxes(id))
      nc  = box%n_cell
      tmp = 1/(12 * box%dr)
      iv  = global_iv_velocity
      call af_gc2_box(tree, id, iv, cc)

      do KJI_DO(1,nc)
         box%cc(IJK, mg%i_rhs) = &
              tmp(1) * (8 * (cc(i+1, j, 1) - cc(i-1, j, 1)) - &
              (cc(i+2, j, 1) - cc(i-2, j, 1))) + &
              tmp(2) * (8 * (cc(i, j+1, 2) - cc(i, j-1, 2)) - &
              (cc(i, j+2, 2) - cc(i, j-2, 2)))
      end do; CLOSE_DO
    end associate
  end subroutine box_set_rhs_4th_order

  subroutine box_correct_velocity(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK, nc, iv(NDIM)
    real(dp)                   :: tmp(NDIM)

    iv = global_iv_velocity
    tmp = 0.5_dp / box%dr

    nc = box%n_cell
    do KJI_DO(1,nc)
       box%cc(IJK, iv(1)) = box%cc(IJK, iv(1)) - tmp(1) * &
            (box%cc(i+1, j, mg%i_phi) - box%cc(i-1, j, mg%i_phi))
       box%cc(IJK, iv(2)) = box%cc(IJK, iv(2)) - tmp(2) * &
            (box%cc(i, j+1, mg%i_phi) - box%cc(i, j-1, mg%i_phi))
    end do; CLOSE_DO
  end subroutine box_correct_velocity

  subroutine bc_ux(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    bc_type = af_bc_dirichlet
    if (nb == af_neighb_highy) then
       bc_val = ux_top
    else
       bc_val = 0.0_dp
    end if
  end subroutine bc_ux

end program incompressible_flow
