!> \example drift_diffusion_$Dd.f90
!>
!> A drift-diffusion example. Diffusion is implemented with centered
!> differences, and for the drift the Koren flux limiter is used. Time stepping
!> is done with the explicit trapezoidal rule.
program drift_diffusion_$Dd
  use m_a$D_types
  use m_a$D_core
  use m_a$D_output
  use m_a$D_utils
  use m_a$D_ghostcell
  use m_a$D_restrict

  implicit none

  integer, parameter  :: box_size   = 8
  integer, parameter  :: i_phi      = 1
  integer, parameter  :: i_phi_old  = 2
  integer, parameter  :: i_err      = 3
  integer, parameter  :: sol_type   = 2
  real(dp), parameter :: domain_len = 2 * acos(-1.0_dp)
  real(dp), parameter :: dr         = domain_len / box_size

  type(a$D_t)         :: tree
  type(ref_info_t)   :: refine_info
  integer            :: ix_list($D, 1)
  integer            :: nb_list(a$D_num_neighbors, 1)
  integer            :: refine_steps, time_steps, output_cnt
  integer            :: i, id, n, n_steps
  real(dp)           :: dt, time, end_time, p_err, n_err, sum_phi
  real(dp)           :: dt_adapt, dt_output
  real(dp)           :: diff_coeff, velocity($D), dr_min($D)
  character(len=100) :: fname
  integer            :: count_rate, t_start, t_end

  print *, "Running drift_diffusion_$Dd"
  print *, "Number of threads", af_get_max_threads()

  ! Initialize tree
  call a$D_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       n_var_cell=3, & ! Number of cell-centered variables
       n_var_face=2, & ! Number of face-centered variables
       dr=dr, &        ! Distance between cells on base level
       cc_names=["phi", "old", "err"]) ! Variable names

  ! Set up geometry
  id             = 1
  ix_list(:, id) = 1            ! Set index of box
  nb_list(:, id) = id           ! Box is periodic, so its own neighbor

  ! Create the base mesh, using the box indices and their neighbor information
  call a$D_set_base(tree, ix_list, nb_list)
  call a$D_print_info(tree)

  output_cnt = 0
  time       = 0
  dt_adapt   = 0.01_dp
  dt_output  = 0.05_dp
  end_time   = 1.5_dp
  diff_coeff = 0.0_dp
  velocity(:) = 0.0_dp
  velocity(1) = 1.0_dp
  velocity(2) = -1.0_dp

  ! Set up the initial conditions
  call system_clock(t_start,count_rate)
  refine_steps=0

  do
     refine_steps=refine_steps+1
     ! We should only set the finest level, but this also works
     call a$D_loop_box(tree, set_initial_condition)

     ! Fill ghost cells for variables i_phi on the sides of all boxes, using
     ! a$D_gc_interp_lim on refinement boundaries: Interpolation between fine
     ! points and coarse neighbors to fill ghost cells near refinement
     ! boundaries. The ghost values are less than twice the coarse values. and
     ! a$D_bc_neumann_zero physical boundaries: fill ghost cells near physical
     ! boundaries using Neumann zero
     call a$D_gc_tree(tree, i_phi, a$D_gc_interp_lim, a$D_bc_neumann_zero)

     ! Adjust the refinement of a tree using refine_routine (see below) for grid
     ! refinement.
     ! Routine a$D_adjust_refinement sets the bit af_bit_new_children for each box
     ! that is refined.  On input, the tree should be balanced. On output,
     ! the tree is still balanced, and its refinement is updated (with at most
     ! one level per call).
     call a$D_adjust_refinement(tree, &           ! tree
          refine_routine, & ! Refinement function
          refine_info, &    ! Information about refinement
          2)                ! Buffer width (in cells)

     ! If no new boxes have been added, exit the loop
     if (refine_info%n_add == 0) exit
  end do
  call system_clock(t_end, count_rate)

  write(*, "(A,i0,A,Es10.3,A)") " Wall-clock time for ", &
       refine_steps, " refinement steps: ", &
       (t_end-t_start) / real(count_rate, dp), " seconds"

  call a$D_print_info(tree)

  ! Restrict the initial conditions Restrict the children of a box to the box
  ! (e.g., in $DD, average the values at the four children to get the value for
  ! the parent)
  call a$D_restrict_tree(tree, i_phi)

  ! Fill ghost cells for variables i_phi on the sides of all boxes, using
  ! a$D_gc_interp_lim on refinement boundaries and a$D_bc_neumann_zero on
  ! physical boundaries
  call a$D_gc_tree(tree, i_phi, a$D_gc_interp_lim, a$D_bc_neumann_zero)

  call system_clock(t_start, count_rate)
  time_steps = 0

  ! Starting simulation
  do
     time_steps = time_steps + 1
     dr_min  = a$D_min_dr(tree)
     dt      = 0.5_dp / (2 * diff_coeff * sum(1/dr_min**2) + &
          sum(abs(velocity/dr_min)) + epsilon(1.0_dp))

     n_steps = ceiling(dt_adapt/dt)
     dt      = dt_adapt / n_steps

     if (output_cnt * dt_output <= time) then
        output_cnt = output_cnt + 1
        write(fname, "(A,I0)") "drift_diffusion_$Dd_", output_cnt

        ! Call procedure set_error (see below) for each box in tree, with argument time
        call a$D_loop_box_arg(tree, set_error, [time])

        ! Write the cell centered data of tree to a vtk unstructured file fname.
        ! Only the leaves of the tree are used
        call a$D_write_vtk(tree, trim(fname), output_cnt, time, &
             ixs_fc=[1], dir="output")

        ! Find maximum and minimum values of cc(..., i_err) and cc(..., i_phi).
        ! By default, only loop over leaves, and ghost cells are not used.
        call a$D_tree_max_cc(tree, i_err, p_err)
        call a$D_tree_min_cc(tree, i_err, n_err)
        call a$D_tree_sum_cc(tree, i_phi, sum_phi)
        write(*,"(2(A,1x,Es12.4,2x))")  &
             " max error:", max(p_err, abs(n_err)), &
             "sum phi:  ", sum_phi
     end if

     if (time > end_time) exit

     do n = 1, n_steps
        ! Copy previous solution
        call a$D_tree_copy_cc(tree, i_phi, i_phi_old)

        ! Two forward Euler steps over dt
        do i = 1, 2
           ! Call procedure fluxes_koren for each id in tree, giving the list of boxes
           call a$D_loop_boxes(tree, fluxes_koren)

           ! Restrict fluxes from children to parents on refinement boundaries.
           call a$D_consistent_fluxes(tree, [i_phi])

           ! Call procedure update_solution (see below) for each box in tree, with argument dt
           call a$D_loop_box_arg(tree, update_solution, [dt])

           ! Restrict variables i_phi to all parent boxes
           call a$D_restrict_tree(tree, i_phi)
        end do

        ! Take average of phi_old and phi
        call a$D_loop_box(tree, average_phi)
        time = time + dt
     end do

     ! Fill ghost cells for variable i_phi
     call a$D_gc_tree(tree, i_phi, a$D_gc_interp_lim, a$D_bc_neumann_zero)

     ! Adjust the refinement of a tree using refine_routine (see below) for grid
     ! refinement. On input, the tree should be balanced. On output, the tree is
     ! still balanced, and its refinement is updated (with at most one level per
     ! call).
     call a$D_adjust_refinement(tree, refine_routine, refine_info, 2)

     ! Prolongation of i_phi values to new children (see below)
     call prolong_to_new_children(tree, refine_info)
  end do

  call system_clock(t_end,count_rate)
  write(*, "(A,I0,A,Es10.3,A)") &
       " Wall-clock time after ",time_steps, &
       " time steps: ", (t_end-t_start) / real(count_rate, dp), &
       " seconds"

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call a$D_destroy(tree)

contains

  !> Set refinement flags for box
  subroutine refine_routine(box, cell_flags)
    type(box$D_t), intent(in) :: box
#if $D == 2
    integer, intent(out)      :: cell_flags(box%n_cell, box%n_cell)
#elif $D == 3
    integer, intent(out)      :: cell_flags(box%n_cell, box%n_cell, box%n_cell)
#endif
    real(dp)                  :: diff
#if $D == 2
    integer                   :: i, j, nc
#elif $D == 3
    integer                   :: i, j, k, nc
#endif

    nc   = box%n_cell

#if $D == 2
    do j = 1, nc
       do i = 1, nc
          diff = abs(box%cc(i+1, j, i_phi) + &
               box%cc(i-1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + &
               box%cc(i, j-1, i_phi) - &
               4 * box%cc(i, j, i_phi)) * box%dr

          if (box%lvl < 2 .or. diff > 0.5e-3_dp .and. box%lvl < 5) then
             cell_flags(i, j) = af_do_ref
          else if (diff < 0.1_dp * 0.1e-3_dp) then
             cell_flags(i, j) = af_rm_ref
          else
             cell_flags(i, j) = af_keep_ref
          end if
       end do
    end do
#elif $D == 3
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             diff = abs(box%cc(i+1, j, k, i_phi) + &
                  box%cc(i-1, j, k, i_phi) + &
                  box%cc(i, j+1, k, i_phi) + &
                  box%cc(i, j-1, k, i_phi) + &
                  box%cc(i, j, k+1, i_phi) + &
                  box%cc(i, j, k-1, i_phi) - &
                  6 * box%cc(i, j, k, i_phi)) * box%dr

             if (box%lvl < 2 .or. diff > 0.5e-3_dp .and. box%lvl < 5) then
                cell_flags(i, j, k) = af_do_ref
             else if (diff < 0.1_dp * 0.1e-3_dp) then
                cell_flags(i, j, k) = af_rm_ref
             else
                cell_flags(i, j, k) = af_keep_ref
             end if
          end do
       end do
    end do
#endif
  end subroutine refine_routine

  !> This routine sets the initial conditions for each box
  subroutine set_initial_condition(box)
    type(box$D_t), intent(inout) :: box
#if $D == 2
    integer                      :: i, j, nc
#elif $D == 3
    integer                      :: i, j, k, nc
#endif
    real(dp)                    :: rr($D)

    nc = box%n_cell
#if $D == 2
    do j = 0, nc+1
       do i = 0, nc+1
          rr = a$D_r_cc(box, [i,j])
          box%cc(i, j, i_phi) = solution(rr, 0.0_dp)
       end do
    end do
#elif $D == 3
    do k = 0, nc+1
       do j = 0, nc+1
          do i = 0, nc+1
             rr = a$D_r_cc(box, [i,j,k])
             box%cc(i, j, k, i_phi) = solution(rr, 0.0_dp)
          end do
       end do
    end do
#endif
  end subroutine set_initial_condition

  !> This routine computes the error in i_phi
  subroutine set_error(box, time)
    type(box$D_t), intent(inout) :: box
    real(dp), intent(in)         :: time(:)
#if $D == 2
    integer                      :: i, j, nc
#elif $D == 3
    integer                      :: i, j, k, nc
#endif
    real(dp)                     :: rr($D)

    nc = box%n_cell
#if $D == 2
    do j = 1, nc
       do i = 1, nc
          rr = a$D_r_cc(box, [i,j])
          box%cc(i, j, i_err) = &
               box%cc(i, j, i_phi) - solution(rr, time(1))
       end do
    end do
#elif $D == 3
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             rr = a$D_r_cc(box, [i,j,k])
             box%cc(i, j, k, i_err) = &
                  box%cc(i, j, k, i_phi) - solution(rr, time(1))
          end do
       end do
    end do
#endif
  end subroutine set_error

  !> This routine calculates the analytic solution in point rr
  function solution(rr, t) result(sol)
    real(dp), intent(in) :: rr($D), t
    real(dp)             :: sol, rr_t($D)

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

  !> Modified implementation of Koren limiter, to avoid division or the min/max
  !> functions, which can be problematic / expensive. In most literature, you
  !> have r = ga / gb (ratio of gradients). Then the limiter phi(r) is
  !> multiplied with gb. With this implementation, you get phi(r) * gb
  elemental function koren_mlim(ga, gb)
    real(dp), intent(in) :: ga  ! Density gradient (numerator)
    real(dp), intent(in) :: gb  ! Density gradient (denominator)
    real(dp), parameter  :: sixth = 1/6.0_dp
    real(dp)             :: koren_mlim, t1, t2

    t1 = ga * ga                ! Two temporary variables,
    t2 = ga * gb                ! so that we do not need sign()

    if (t2 <= 0) then
       ! ga and gb have different sign: local minimum/maximum
       koren_mlim = 0
    else if (t1 >= 2.5_dp * t2) then
       ! (1+2*ga/gb)/6 => 1, limiter has value 1
       koren_mlim = gb
    else if (t1 > 0.25_dp * t2) then
       ! 1 > ga/gb > 1/4, limiter has value (1+2*ga/gb)/6
       koren_mlim = sixth * (gb + 2*ga)
    else
       ! 0 < ga/gb < 1/4, limiter has value ga/gb
       koren_mlim = ga
    end if
  end function koren_mlim

  !> This routine computes the x-fluxes and y-fluxes interior (advective part)
  !> with the Koren limiter
  subroutine fluxes_koren(boxes, id)
    use m_a$D_prolong
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)          :: id
    real(dp)                     :: gradp, gradc, gradn
    real(dp)                     :: inv_dr
    integer                      :: dim, dix($D)
#if $D == 2
    real(dp)                     :: cc(-1:boxes(id)%n_cell+2, &
         -1:boxes(id)%n_cell+2)
    integer                      :: i, j, nc
#elif $D == 3
    real(dp)                     :: cc(-1:boxes(id)%n_cell+2, &
         -1:boxes(id)%n_cell+2, -1:boxes(id)%n_cell+2)
    integer                      :: i, j, k, nc
#endif

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    call a$D_gc_box(boxes, id, i_phi, a$D_gc_interp_lim, a$D_bc_neumann_zero)

    ! Get a second layer of ghost cell data (the 'normal' routines give just one
    ! layer of ghost cells). Use a$D_gc2_prolong_linear on refinement boundaries and
    ! a$D_bc2_neumann_zero on physical boundaries.
    call a$D_gc2_box(boxes, &      ! List of all the boxes
         id, &                     ! Id of box for which we set ghost cells
         i_phi, &                  ! Variable for which ghost cells are set
         a$D_gc2_prolong_linear, & ! Procedure called at refinement boundaries
         a$D_bc2_neumann_zero, &   ! Procedure called at physical boundaries
         cc, &                  ! The enlarged box with ghost cells
         nc)                       ! box%n_cell

    do dim = 1, $D
       dix(:) = 0
       dix(dim) = 1

#if $D == 2
       do j = 1, nc+dix(2)
          do i = 1, nc+dix(1)
             gradc = cc(i, j) - cc(i-dix(1), j-dix(2))

             if (velocity(dim) < 0.0_dp) then
                gradn = cc(i+1, j) - cc(i, j)
                boxes(id)%fc(i, j, dim, i_phi) = velocity(dim) * &
                     (cc(i, j) - koren_mlim(gradc, gradn))
             else                  ! velocity(dim) > 0
                gradp = cc(i-dix(1), j-dix(2)) - cc(i-2*dix(1), j-2*dix(2))
                boxes(id)%fc(i, j, dim, i_phi) = velocity(dim) * &
                     (cc(i-dix(1), j-dix(2)) + koren_mlim(gradc, gradp))
             end if

             ! Diffusive part with 2-nd order explicit method
             boxes(id)%fc(i, j, dim, i_phi) = &
                  boxes(id)%fc(i, j, dim, i_phi) - &
                  diff_coeff * gradc * inv_dr
          end do
       end do
#elif $D == 3
       do k = 1, nc+dix(3)
          do j = 1, nc+dix(2)
             do i = 1, nc+dix(1)
                gradc = cc(i, j, k) - cc(i-dix(1), j-dix(2), k-dix(3))

                if (velocity(dim) < 0.0_dp) then
                   gradn = cc(i+dix(1), j+dix(2), k+dix(3)) - cc(i, j, k)
                   boxes(id)%fc(i, j, k, dim, i_phi) = velocity(dim) * &
                        (cc(i, j, k) - koren_mlim(gradc, gradn))
                else                  ! velocity(dim) > 0
                   gradp = cc(i-dix(1), j-dix(2), k-dix(3)) - &
                        cc(i-2*dix(1), j-2*dix(2), k-2*dix(3))
                   boxes(id)%fc(i, j, k, dim, i_phi) = velocity(dim) * &
                        (cc(i-dix(1), j-dix(2), k-dix(3)) + &
                        koren_mlim(gradc, gradp))
                end if

                ! Diffusive part with 2-nd order explicit method
                boxes(id)%fc(i, j, k, dim, i_phi) = &
                     boxes(id)%fc(i, j, k, dim, i_phi) - &
                     diff_coeff * gradc * inv_dr
             end do
          end do
       end do
#endif
    end do
  end subroutine fluxes_koren

  !> This routine computes the update of the solution per box
  subroutine update_solution(box, dt)
    type(box$D_t), intent(inout) :: box
    real(dp), intent(in)         :: dt(:)
    real(dp)                     :: inv_dr
#if $D == 2
    integer                      :: i, j
#elif $D == 3
    integer                      :: i, j, k
#endif
    integer                      :: nc

    nc     = box%n_cell
    inv_dr = 1/box%dr

#if $D == 2
    forall (i = 1:nc, j = 1:nc)
       box%cc(i, j, i_phi) = box%cc(i, j, i_phi) + dt(1) * inv_dr * ( &
            box%fc(i, j, 1, i_phi) - box%fc(i+1, j, 1, i_phi) + &
            box%fc(i, j, 2, i_phi) - box%fc(i, j+1, 2, i_phi))
    end forall
#elif $D == 3
    forall (i = 1:nc, j = 1:nc, k = 1:nc)
       box%cc(i, j, k, i_phi) = box%cc(i, j, k, i_phi) + dt(1) * inv_dr * ( &
            box%fc(i, j, k, 1, i_phi) - box%fc(i+1, j, k, 1, i_phi) + &
            box%fc(i, j, k, 2, i_phi) - box%fc(i, j+1, k, 2, i_phi) + &
            box%fc(i, j, k, 3, i_phi) - box%fc(i, j, k+1, 3, i_phi))
    end forall
#endif
  end subroutine update_solution

  !> This routine computes the update of the solution per box
  subroutine average_phi(box)
    type(box$D_t), intent(inout) :: box

#if $D == 2
    box%cc(:, :, i_phi) = 0.5_dp * (box%cc(:, :, i_phi) + &
         box%cc(:, :, i_phi_old))
#elif $D == 3
    box%cc(:, :, :, i_phi) = 0.5_dp * (box%cc(:, :, :, i_phi) + &
         box%cc(:, :, :, i_phi_old))
#endif
  end subroutine average_phi

  ! Linear prolongation of i_phi values to new children
  subroutine prolong_to_new_children(tree, refine_info)
    use m_a$D_prolong
    type(a$D_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: refine_info
    integer                      :: lvl, i, id, p_id

    do lvl = 1, tree%highest_lvl
       !$omp parallel do private(id, p_id)
       do i = 1, size(refine_info%lvls(lvl)%add)
          id = refine_info%lvls(lvl)%add(i)
          p_id = tree%boxes(id)%parent

          ! Linear prolongation will not strictly conserve phi
          call a$D_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_phi)
       end do
       !$omp end parallel do

       ! Fill ghost cells for variables i_phi on the sides of a box, using
       ! a$D_gc_interp_lim on refinement boundaries and a$D_bc_neumann_zero
       ! on physical boundaries
       call a$D_gc_ids(tree%boxes, refine_info%lvls(lvl)%add, i_phi, &
            a$D_gc_interp_lim, a$D_bc_neumann_zero)
    end do
  end subroutine prolong_to_new_children

end program drift_diffusion_$Dd
