#include "../src/cpp_macros.h"
!> \example advection_v2.f90
!>
!> An advection example
program advection
  use m_af_all

  implicit none

  integer, parameter  :: box_size   = 8
  integer             :: i_phi
  integer             :: i_err
  integer             :: i_flux
  integer, parameter  :: sol_type   = 2
  real(dp), parameter :: domain_len = 2 * acos(-1.0_dp)

  type(af_t)         :: tree
  type(ref_info_t)   :: refine_info
  integer            :: refine_steps, output_cnt
  real(dp)           :: dt, dt_lim, time, time_prev_refine
  real(dp)           :: end_time, err, sum_err2
  real(dp)           :: sum_phi, sum_phi_t0
  real(dp)           :: dt_amr, dt_output
  real(dp)           :: velocity(NDIM)
  character(len=100) :: fname

  print *, "Running advection_" // DIMNAME // ""
  print *, "Number of threads", af_get_max_threads()

  ! Add variables to the mesh. This is a scalar advection example with second
  ! order time stepping, which is why there are two copies of phi.
  call af_add_cc_variable(tree, "phi", ix=i_phi, n_copies=2)
  call af_add_cc_variable(tree, "err", ix=i_err)
  call af_add_fc_variable(tree, "flux", ix=i_flux)

  call af_set_cc_methods(tree, i_phi, af_bc_neumann_zero, &
       prolong=af_prolong_limit)

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [DTIMES(domain_len)], &
       [DTIMES(box_size)], &
       periodic=[DTIMES(.true.)])

  output_cnt  = 0
  time        = 0
  dt_amr      = 0.01_dp
  dt_output   = 0.5_dp
  end_time    = 5.0_dp
  velocity(:) = -0.5_dp
  velocity(1) = 1.0_dp

  ! Set up the initial conditions
  refine_steps=0

  do
     refine_steps=refine_steps+1
     ! Set initial conditions on all boxes
     call af_loop_box(tree, set_initial_condition)

     ! Fill ghost cells for variables i_phi
     call af_gc_tree(tree, [i_phi])

     ! Adjust the refinement of a tree using refine_routine
     call af_adjust_refinement(tree, refine_routine, refine_info, 1)

     ! If no new boxes have been added, exit the loop
     if (refine_info%n_add == 0) exit
  end do

  call af_print_info(tree)

  ! Restrict the initial conditions
  call af_restrict_tree(tree, [i_phi])

  ! Fill ghost cells for variables i_phi on the sides of all boxes
  call af_gc_tree(tree, [i_phi])

  call af_tree_sum_cc(tree, i_phi, sum_phi_t0)

  dt = 0.5_dp / (sum(abs(velocity/af_lvl_dr(tree, tree%highest_lvl))) + &
       epsilon(1.0_dp))
  time_prev_refine = time

  ! Starting simulation
  do
     if (output_cnt * dt_output <= time) then
        output_cnt = output_cnt + 1
        write(fname, "(A,I0)") "output/advection_v2_" // DIMNAME // "_", output_cnt

        ! Call procedure set_error (see below) for each box in tree, with argument time
        call af_loop_box_arg(tree, set_error, [time])

        ! Write the cell centered data of tree to a vtk unstructured file fname.
        ! Only the leaves of the tree are used
        call af_write_silo(tree, trim(fname), output_cnt, time)

        ! Find maximum and minimum values of cc(..., i_err) and cc(..., i_phi).
        ! By default, only loop over leaves, and ghost cells are not used.
        call af_tree_maxabs_cc(tree, i_err, err)
        call af_tree_sum_cc(tree, i_phi, sum_phi)
        call af_tree_sum_cc(tree, i_err, sum_err2, power=2)
        write(*,"(3(A,Es16.8))")  &
             " max error:", err, &
             " mean error:", sqrt(sum_err2/af_total_volume(tree)), &
             " conservation error:  ", (sum_phi - sum_phi_t0) / sum_phi_t0
     end if

     if (time > end_time) exit

     call af_advance(tree, dt, dt_lim, time, [i_phi], &
          af_midpoint_method, forward_euler)
     dt = 0.8_dp * dt_lim

     if (time > time_prev_refine + dt_amr) then
        call af_gc_tree(tree, [i_phi])
        call af_adjust_refinement(tree, refine_routine, refine_info, 1)
        time_prev_refine = time
     end if
  end do

contains

  !> [refine_routine]
  !> Set refinement flags for box
  subroutine refine_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)      :: cell_flags(DTIMES(box%n_cell))
    real(dp)                  :: diff
    integer                   :: IJK, nc

    nc   = box%n_cell

    do KJI_DO(1,nc)
#if NDIM == 1
       diff = abs(box%dr(1) * (box%cc(i+1, i_phi) + &
            box%cc(i-1, i_phi) - 2 * box%cc(i, i_phi)))
#elif NDIM == 2
       diff = abs(box%dr(1) * (box%cc(i+1, j, i_phi) + &
            box%cc(i-1, j, i_phi) - 2 * box%cc(i, j, i_phi)) + &
            box%dr(2) * (box%cc(i, j+1, i_phi) + &
            box%cc(i, j-1, i_phi) - 2 * box%cc(i, j, i_phi)))
#elif NDIM == 3
       diff = abs(box%dr(1) * (box%cc(i+1, j, k, i_phi) + &
            box%cc(i-1, j, k, i_phi) - 2 * box%cc(i, j, k, i_phi)) + &
            box%dr(2) * (box%cc(i, j+1, k, i_phi) + &
            box%cc(i, j-1, k, i_phi) - 2 * box%cc(i, j, k, i_phi)) + &
            box%dr(3) * (box%cc(i, j, k+1, i_phi) + &
            box%cc(i, j, k-1, i_phi) - 2 * box%cc(i, j, k, i_phi)))
#endif

       if (box%lvl < 2 .or. diff > 2.0e-3_dp .and. box%lvl < 5) then
          cell_flags(IJK) = af_do_ref
       else if (diff < 0.1_dp * 2.0e-3_dp) then
          cell_flags(IJK) = af_rm_ref
       else
          cell_flags(IJK) = af_keep_ref
       end if
    end do; CLOSE_DO
  end subroutine refine_routine
  !> [refine_routine]

  !> This routine sets the initial conditions for each box
  subroutine set_initial_condition(box)
    type(box_t), intent(inout) :: box
    integer                      :: IJK, nc
    real(dp)                     :: rr(NDIM)

    nc = box%n_cell
    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, i_phi) = solution(rr, 0.0_dp)
    end do; CLOSE_DO
  end subroutine set_initial_condition

  !> This routine computes the error in i_phi
  subroutine set_error(box, time)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: time(:)
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM)

    nc = box%n_cell
    do KJI_DO(1,nc)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, i_err) = &
            box%cc(IJK, i_phi) - solution(rr, time(1))
    end do; CLOSE_DO
  end subroutine set_error

  !> This routine calculates the analytic solution in point rr
  function solution(rr, t) result(sol)
    real(dp), intent(in) :: rr(NDIM), t
    real(dp)             :: sol, rr_t(NDIM)

    rr_t = rr - velocity * t

    select case (sol_type)
    case (1)
#if NDIM > 1
       sol = sin(0.5_dp * rr_t(1))**8 * cos(0.5_dp * rr_t(2))**8
#else
       sol = sin(0.5_dp * rr_t(1))**8
#endif
    case (2)
       rr_t = modulo(rr_t, domain_len) / domain_len
       if (norm2(rr_t - 0.5_dp) < 0.1_dp) then
          sol = 1.0_dp
       else
          sol = 0.0_dp
       end if
    end select
  end function solution

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
    real(dp)                  :: wmax(NDIM)

    call flux_generic_tree(tree, 1, [i_phi+s_deriv], [i_flux], wmax, &
         max_wavespeed, get_flux, flux_dummy_conversion, flux_dummy_conversion)
    call flux_update_densities(tree, dt, 1, [i_phi], [i_flux], &
         s_deriv, n_prev, s_prev, w_prev, s_out, flux_dummy_source)

    ! Compute maximal time step
    dt_lim = 1.0_dp / sum(wmax/af_lvl_dr(tree, tree%highest_lvl))
  end subroutine forward_euler

  subroutine max_wavespeed(n_values, n_var, flux_dim, u, w)
    integer, intent(in)   :: n_values !< Number of cell faces
    integer, intent(in)   :: n_var    !< Number of variables
    integer, intent(in)   :: flux_dim !< In which dimension fluxes are computed
    real(dp), intent(in)  :: u(n_values, n_var) !< Primitive variables
    real(dp), intent(out) :: w(n_values) !< Maximum speed

    w = abs(velocity(flux_dim))
  end subroutine max_wavespeed

  subroutine get_flux(n_values, n_var, flux_dim, u, flux, box, line_ix)
    integer, intent(in)     :: n_values !< Number of cell faces
    integer, intent(in)     :: n_var    !< Number of variables
    integer, intent(in)     :: flux_dim !< In which dimension fluxes are computed
    real(dp), intent(in)    :: u(n_values, n_var)
    real(dp), intent(out)   :: flux(n_values, n_var)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: line_ix(NDIM-1)

    flux = u * velocity(flux_dim)
  end subroutine get_flux

end program
