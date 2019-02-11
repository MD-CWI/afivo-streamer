#include "../src/cpp_macros.h"
!> \example advection_Xd.f90
!>
!> An advection example using the Koren flux limiter. Time stepping is done with
!> the explicit trapezoidal rule.
program advection_Xd
  use m_af_all

  implicit none

  integer, parameter  :: box_size   = 8
  integer             :: i_phi
  integer             :: i_phi_old
  integer             :: i_err
  integer             :: i_flux
  integer, parameter  :: sol_type   = 1
  ! Set coord_type to af_cyl to test conservation in cyl. coordinates
  integer, parameter  :: coord_type = af_xyz
  real(dp), parameter :: domain_len = 2 * acos(-1.0_dp)

  type(af_t)         :: tree
  type(ref_info_t)   :: refine_info
  integer            :: refine_steps, time_steps, output_cnt
  integer            :: i, n, n_steps
  real(dp)           :: dt, time, end_time, p_err, n_err
  real(dp)           :: sum_phi, sum_phi_t0
  real(dp)           :: dt_adapt, dt_output
  real(dp)           :: velocity(NDIM), dr_min(NDIM)
  character(len=100) :: fname
  integer            :: count_rate, t_start, t_end

  print *, "Running advection_" // DIMNAME // ""
  print *, "Number of threads", af_get_max_threads()

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "old", ix=i_phi_old)
  call af_add_cc_variable(tree, "err", ix=i_err)
  call af_add_fc_variable(tree, "flux", ix=i_flux)

  call af_set_cc_methods(tree, i_phi, af_bc_neumann_zero, af_gc_interp_lim, &
       prolong=af_prolong_limit)

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [DTIMES(domain_len)], &
       [DTIMES(box_size)], &
       periodic=[DTIMES(.true.)], &
       coord=coord_type)

  output_cnt = 0
  time       = 0
  dt_adapt   = 0.01_dp
  dt_output  = 0.5_dp
  end_time   = 5.0_dp
  velocity(:) = 0.0_dp
  velocity(1) = 1.0_dp
  velocity(2) = -1.0_dp

  ! Set up the initial conditions
  call system_clock(t_start,count_rate)
  refine_steps=0

  do
     refine_steps=refine_steps+1
     ! We should only set the finest level, but this also works
     call af_loop_box(tree, set_initial_condition)

     ! Fill ghost cells for variables i_phi on the sides of all boxes, using
     ! af_gc_interp_lim on refinement boundaries: Interpolation between fine
     ! points and coarse neighbors to fill ghost cells near refinement
     ! boundaries. The ghost values are less than twice the coarse values. and
     ! af_bc_neumann_zero physical boundaries: fill ghost cells near physical
     ! boundaries using Neumann zero
     call af_gc_tree(tree, i_phi, af_gc_interp_lim, af_bc_neumann_zero)

     ! Adjust the refinement of a tree using refine_routine (see below) for grid
     ! refinement.
     ! Routine af_adjust_refinement sets the bit af_bit_new_children for each box
     ! that is refined.  On input, the tree should be balanced. On output,
     ! the tree is still balanced, and its refinement is updated (with at most
     ! one level per call).
     call af_adjust_refinement(tree, &           ! tree
          refine_routine, & ! Refinement function
          refine_info, &    ! Information about refinement
          1)                ! Buffer width (in cells)

     ! If no new boxes have been added, exit the loop
     if (refine_info%n_add == 0) exit
  end do
  call system_clock(t_end, count_rate)

  write(*, "(A,i0,A,Es10.3,A)") " Wall-clock time for ", &
       refine_steps, " refinement steps: ", &
       (t_end-t_start) / real(count_rate, dp), " seconds"

  call af_print_info(tree)

  ! Restrict the initial conditions Restrict the children of a box to the box
  ! (e.g., in NDIMD, average the values at the four children to get the value for
  ! the parent)
  call af_restrict_tree(tree, i_phi)

  ! Fill ghost cells for variables i_phi on the sides of all boxes, using
  ! af_gc_interp_lim on refinement boundaries and af_bc_neumann_zero on
  ! physical boundaries
  call af_gc_tree(tree, i_phi, af_gc_interp_lim, af_bc_neumann_zero)

  call system_clock(t_start, count_rate)
  time_steps = 0

  call af_tree_sum_cc(tree, i_phi, sum_phi_t0)

  ! Starting simulation
  do
     time_steps = time_steps + 1
     dr_min  = af_min_dr(tree)
     dt      = 0.5_dp / (sum(abs(velocity/dr_min)) + epsilon(1.0_dp))

     n_steps = ceiling(dt_adapt/dt)
     dt      = dt_adapt / n_steps

     if (output_cnt * dt_output <= time) then
        output_cnt = output_cnt + 1
        write(fname, "(A,I0)") "advection_" // DIMNAME // "_", output_cnt

        ! Call procedure set_error (see below) for each box in tree, with argument time
        call af_loop_box_arg(tree, set_error, [time])

        ! Write the cell centered data of tree to a vtk unstructured file fname.
        ! Only the leaves of the tree are used
        call af_write_silo(tree, trim(fname), output_cnt, time, dir="output")

        ! Find maximum and minimum values of cc(..., i_err) and cc(..., i_phi).
        ! By default, only loop over leaves, and ghost cells are not used.
        call af_tree_max_cc(tree, i_err, p_err)
        call af_tree_min_cc(tree, i_err, n_err)
        call af_tree_sum_cc(tree, i_phi, sum_phi)
        write(*,"(2(A,1x,Es12.4,2x))")  &
             " max error:", max(p_err, abs(n_err)), &
             "conservation error:  ", (sum_phi - sum_phi_t0) / sum_phi_t0
     end if

     if (time > end_time) exit

     do n = 1, n_steps
        ! Copy previous solution
        call af_tree_copy_cc(tree, i_phi, i_phi_old)

        ! Two forward Euler steps over dt
        do i = 1, 2
           ! Call procedure fluxes_koren for each id in tree, giving the list of boxes
           call af_loop_boxes(tree, fluxes_koren)

           ! Restrict fluxes from children to parents on refinement boundaries.
           call af_consistent_fluxes(tree, [i_phi])

           ! Call procedure update_solution (see below) for each box in tree, with argument dt
           call af_loop_box_arg(tree, update_solution, [dt])

           ! Restrict variables i_phi to all parent boxes
           call af_restrict_tree(tree, i_phi)
        end do

        ! Take average of phi_old and phi
        call af_loop_box(tree, average_phi)
        time = time + dt
     end do

     ! Fill ghost cells for variable i_phi
     call af_gc_tree(tree, i_phi, af_gc_interp_lim, af_bc_neumann_zero)

     !> [adjust_refinement]
     ! Adjust the refinement of a tree using refine_routine (see below) for grid
     ! refinement. On input, the tree should be balanced. On output, the tree is
     ! still balanced, and its refinement is updated (with at most one level per
     ! call).
     call af_adjust_refinement(tree, refine_routine, refine_info, 1)
     !> [adjust_refinement]
  end do

  call system_clock(t_end,count_rate)
  write(*, "(A,I0,A,Es10.3,A)") &
       " Wall-clock time after ",time_steps, &
       " time steps: ", (t_end-t_start) / real(count_rate, dp), &
       " seconds"

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call af_destroy(tree)

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
#if NDIM == 2
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

       if (box%lvl < 2 .or. diff > 2.0e-3_dp .and. box%lvl < 3) then
          cell_flags(IJK) = af_do_ref
       else if (diff < 0.1_dp * 0.1e-3_dp) then
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
    real(dp), intent(in)         :: time(:)
    integer                      :: IJK, nc
    real(dp)                     :: rr(NDIM)

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
       sol = sin(rr_t(1))**4 * cos(rr_t(2))**4
    case (2)
       rr_t = modulo(rr_t, domain_len) / domain_len
       if (norm2(rr_t - 0.5_dp) < 0.1_dp) then
          sol = 1.0_dp
       else
          sol = 0.0_dp
       end if
    end select
  end function solution

  !> This routine computes the x-fluxes and y-fluxes interior (advective part)
  !> with the Koren limiter
  subroutine fluxes_koren(boxes, id)
    use m_af_flux_schemes
    type(box_t), intent(inout) :: boxes(:)
    integer, intent(in)          :: id
    integer                      :: nc
    real(dp), allocatable        :: cc(DTIMES(:))
    real(dp), allocatable        :: v(DTIMES(:), :)

    nc     = boxes(id)%n_cell
    allocate(cc(DTIMES(-1:nc+2)))
    allocate(v(DTIMES(1:nc+1), NDIM))

    call af_gc_box(boxes, id, i_phi, af_gc_interp_lim, af_bc_neumann_zero)

    ! Get a second layer of ghost cell data (the 'normal' routines give just one
    ! layer of ghost cells). Use af_gc2_prolong_linear on refinement boundaries and
    ! af_bc2_neumann_zero on physical boundaries.
    call af_gc2_box(boxes, &      ! List of all the boxes
         id, &                     ! Id of box for which we set ghost cells
         i_phi, &                  ! Variable for which ghost cells are set
         af_gc2_prolong_linear, & ! Procedure called at refinement boundaries
         af_bc2_neumann_zero, &   ! Procedure called at physical boundaries
         cc, &                  ! The enlarged box with ghost cells
         nc)                       ! box%n_cell

#if NDIM == 2
    v(:, :, 1) = velocity(1)
    v(:, :, 2) = velocity(2)

    call flux_koren_2d(cc, v, nc, 2)
    boxes(id)%fc(:, :, :, i_phi) = v
#elif NDIM == 3
    v(:, :, :, 1) = velocity(1)
    v(:, :, :, 2) = velocity(2)
    v(:, :, :, 3) = velocity(3)

    call flux_koren_3d(cc, v, nc, 2)
    boxes(id)%fc(:, :, :, :, i_phi) = v
#endif

  end subroutine fluxes_koren

  !> This routine computes the update of the solution per box
  subroutine update_solution(box, dt)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)         :: dt(:)
    real(dp)                     :: inv_dr(NDIM)
    integer                      :: IJK, nc
#if NDIM == 2
    real(dp)                     :: rfac(2, box%n_cell)
#endif

    nc     = box%n_cell
    inv_dr = 1/box%dr

#if NDIM == 2
    if (coord_type == af_cyl) then
       call af_cyl_flux_factors(box, rfac)
       do j = 1, nc
          do i = 1, nc
             box%cc(i, j, i_phi) = box%cc(i, j, i_phi) + dt(1) * ( &
                  inv_dr(1) * (rfac(1, i) * box%fc(i, j, 1, i_phi) - &
                  rfac(2, i) * box%fc(i+1, j, 1, i_phi)) &
                  + inv_dr(2) * &
                  (box%fc(i, j, 2, i_phi) - box%fc(i, j+1, 2, i_phi)))
          end do
       end do
    else
       do j = 1, nc
          do i = 1, nc
             box%cc(i, j, i_phi) = box%cc(i, j, i_phi) + dt(1) * ( &
                  inv_dr(1) * &
                  (box%fc(i, j, 1, i_phi) - box%fc(i+1, j, 1, i_phi)) &
                  + inv_dr(2) * &
                  (box%fc(i, j, 2, i_phi) - box%fc(i, j+1, 2, i_phi)))
          end do
       end do
    end if
#elif NDIM == 3
    forall (i = 1:nc, j = 1:nc, k = 1:nc)
       box%cc(i, j, k, i_phi) = box%cc(i, j, k, i_phi) + dt(1) * ( &
            inv_dr(1) * &
            (box%fc(i, j, k, 1, i_phi) - box%fc(i+1, j, k, 1, i_phi)) + &
            inv_dr(2) * &
            (box%fc(i, j, k, 2, i_phi) - box%fc(i, j+1, k, 2, i_phi)) + &
            inv_dr(3) * &
            (box%fc(i, j, k, 3, i_phi) - box%fc(i, j, k+1, 3, i_phi)))
    end forall
#endif
  end subroutine update_solution

  !> This routine computes the update of the solution per box
  subroutine average_phi(box)
    type(box_t), intent(inout) :: box

    box%cc(DTIMES(:), i_phi) = 0.5_dp * (box%cc(DTIMES(:), i_phi) + &
         box%cc(DTIMES(:), i_phi_old))
  end subroutine average_phi

end program advection_Xd
