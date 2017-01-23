!> \example drift_diffusion_3d.f90
!>
!> A drift-diffusion example for m_a3_types
!> @TODO: document this
program drift_diffusion_3d
  use m_a3_types
  use m_a3_core
  use m_a3_output
  use m_a3_utils
  use m_a3_ghostcell
  use m_a3_restrict

  implicit none

  integer, parameter :: box_size  = 8
  integer, parameter :: i_phi     = 1
  integer, parameter :: i_phi_old = 2
  integer, parameter :: i_err     = 3

  real(dp), parameter :: domain_len = 2 * acos(-1.0_dp)
  real(dp), parameter :: dr = domain_len / box_size

  type(a3_t)         :: tree
  type(ref_info_t)   :: refine_info
  integer            :: ix_list(3, 1)
  integer            :: nb_list(6, 1)
  integer            :: refine_steps, time_steps, output_cnt
  integer            :: i, id, n, n_steps
  real(dp)           :: dt, time, end_time, p_err, n_err, sum_phi
  real(dp)           :: dt_adapt, dt_output
  real(dp)           :: diff_coeff, vel_x, vel_y, vel_z, dr_min(3)
  character(len=100) :: fname
  integer            :: count_rate, t_start, t_end
  integer            :: time_step_method = 2

  print *, "Running drift_diffusion_3d"
  print *, "Number of threads", af_get_max_threads()

  ! Initialize tree
  call a3_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       n_var_cell=3, & ! Number of cell-centered variables
       n_var_face=2, & ! Number of face-centered variables
       dr=dr, &        ! Distance between cells on base level
       cc_names=["phi", "old", "err"]) ! Variable names

  ! Set up geometry
  id             = 1
  ix_list(:, id) = [1,1,1] ! Set index of box
  nb_list(:, id) = id      ! Box is periodic, so its own neighbor

  ! Create the base mesh, using the box indices and their neighbor information
  call a3_set_base(tree, ix_list, nb_list)
  call a3_print_info(tree)

  output_cnt = 0
  time       = 0
  dt_adapt   = 0.01_dp
  dt_output  = 0.05_dp
  end_time   = 1.5_dp
  diff_coeff = 0.0_dp
  vel_x      = 1.0_dp
  vel_y      = 0.0_dp
  vel_z      = 1.0_dp

  ! Set up the initial conditions
  call system_clock(t_start,count_rate)
  refine_steps=0
  do
     refine_steps=refine_steps+1
     ! We should only set the finest level, but this also works
     call a3_loop_box(tree, set_initial_condition)

     ! Fill ghost cells for variables i_phi on the sides of all boxes, using
     ! a3_gc_interp_lim on refinement boundaries: Interpolation between fine
     ! points and coarse neighbors to fill ghost cells near refinement
     ! boundaries. The ghost values are less than twice the coarse values. and
     ! a3_bc_neumann_zero physical boundaries: fill ghost cells near physical
     ! boundaries using Neumann zero
     call a3_gc_tree(tree, i_phi, a3_gc_interp_lim, a3_bc_neumann_zero)

     ! Adjust the refinement of a tree using refine_routine (see below) for grid
     ! refinement. 
     ! Routine a3_adjust_refinement sets the bit af_bit_new_children for each box
     ! that is refined.  On input, the tree should be balanced. On output,
     ! the tree is still balanced, and its refinement is updated (with at most
     ! one level per call).
     call a3_adjust_refinement(tree, &           ! tree
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

  call a3_print_info(tree)

  ! Restrict the initial conditions
  ! Restrict the children of a box to the box (e.g., in 3D, average the values
  ! at the four children to get the value for the parent)
  call a3_restrict_tree(tree, i_phi)

  ! Fill ghost cells for variables i_phi on the sides of all boxes, using
  ! a3_gc_interp_lim on refinement boundaries and a3_bc_neumann_zero on physical boundaries
  call a3_gc_tree(tree, i_phi, a3_gc_interp_lim, a3_bc_neumann_zero)

  select case (time_step_method)
  case (1)
     print *,"Time stepping: Forward Euler"
  case (2)
     print*, "Time stepping: Two forward Euler steps over dt"
  end select

  call system_clock(t_start, count_rate)
  time_steps = 0

  ! Starting simulation
  do
     time_steps = time_steps + 1
     dr_min  = a3_min_dr(tree)
     dt      = 0.5_dp / (2 * diff_coeff * sum(1/dr_min**2) + &
          sum(abs([vel_x, vel_y, vel_z]) / dr_min ) + epsilon(1.0_dp))

     n_steps = ceiling(dt_adapt/dt)
     dt      = dt_adapt / n_steps

     if (output_cnt * dt_output <= time) then
        output_cnt = output_cnt + 1
        write(fname, "(A,I0)") "drift_diffusion_3d_", output_cnt

        ! Call procedure set_error (see below) for each box in tree, with argument time
        call a3_loop_box_arg(tree, set_error, [time])

        ! Write the cell centered data of tree to a vtk unstructured file fname.
        ! Only the leaves of the tree are used
        call a3_write_vtk(tree, trim(fname), output_cnt, time, &
             ixs_fc=[1], dir="output")

        ! Find maximum and minimum values of cc(..., i_err) and cc(..., i_phi).
        ! By default, only loop over leaves, and ghost cells are not used.
        call a3_tree_max_cc(tree, i_err, p_err)
        call a3_tree_min_cc(tree, i_err, n_err)
        call a3_tree_sum_cc(tree, i_phi, sum_phi)
        write(*,"(2(A,1x,Es12.4,2x))")  &
             " max error:", max(p_err, abs(n_err)), &
             "sum phi:  ", sum_phi
     end if

     if (time > end_time) exit

     select case (time_step_method)
     case (1)
        ! Forward Euler
        do n = 1, n_steps
           ! Call procedure fluxes_centdif for each id in tree, giving the list of boxes
           call a3_loop_boxes_arg(tree, fluxes_centdif, &
                [diff_coeff, vel_x, vel_y, vel_z])

           ! Restrict fluxes from children to parents on refinement boundaries.
           call a3_consistent_fluxes(tree, [1])

           ! Call procedure update_solution (see below) for each box in tree, with argument dt
           call a3_loop_box_arg(tree, update_solution, [dt])

           ! Restrict variables i_phi to all parent boxes
           call a3_restrict_tree(tree, i_phi)

           ! Fill ghost cells for variables i_phi on the sides of all boxes, using
           ! a3_gc_interp_lim on refinement boundaries and a3_bc_neumann_zero on physical boundaries
           call a3_gc_tree(tree, i_phi, a3_gc_interp_lim, a3_bc_neumann_zero)
           time = time + dt
        end do
     case (2)
        do n = 1, n_steps
           ! Copy previous solution
           call a3_tree_copy_cc(tree, i_phi, i_phi_old)

           ! Two forward Euler steps over dt
           do i = 1, 2
              ! Call procedure fluxes_koren for each id in tree, giving the list of boxes
              call a3_loop_boxes_arg(tree, fluxes_koren, &
                   [diff_coeff, vel_x, vel_y, vel_z])

              ! Restrict fluxes from children to parents on refinement boundaries.
              call a3_consistent_fluxes(tree, [1])

              ! Call procedure update_solution (see below) for each box in tree, with argument dt
              call a3_loop_box_arg(tree, update_solution, [dt])

              ! Restrict variables i_phi to all parent boxes
              call a3_restrict_tree(tree, i_phi)

              ! Fill ghost cells for variables i_phi on the sides of all boxes, using
              ! a3_gc_interp_lim on refinement boundaries and a3_bc_neumann_zero on physical boundaries
              call a3_gc_tree(tree, i_phi, a3_gc_interp_lim, a3_bc_neumann_zero)
           end do

           ! Take average of phi_old and phi
           call a3_loop_box(tree, average_phi)
           time = time + dt
        end do
     end select

     ! Adjust the refinement of a tree using refine_routine (see below) for grid
     ! refinement. 
     ! Routine a3_adjust_refinement sets the bit af_bit_new_children for each box
     ! that is refined.  On input, the tree should be balanced. On output,
     ! the tree is still balanced, and its refinement is updated (with at most
     ! one level per call).
     call a3_adjust_refinement(tree, refine_routine, refine_info, 2)

     ! Linear prolongation of i_phi values to new children (see below)
     call prolong_to_new_children(tree, refine_info)

     ! Fill ghost cells for variables i_phi on the sides of all boxes, using
     ! a3_gc_interp_lim on refinement boundaries and a3_bc_neumann_zero on physical boundaries
     call a3_gc_tree(tree, i_phi, a3_gc_interp_lim, a3_bc_neumann_zero)

     !> Reorder and resize the list of boxes
     call a3_tidy_up(tree, 0.8_dp, 10000)
  end do

  call system_clock(t_end,count_rate)
  write(*, "(A,I0,A,Es10.3,A)") &
       " Wall-clock time after ",time_steps, &
       " time steps: ", (t_end-t_start) / real(count_rate, dp), &
       " seconds"

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call a3_destroy(tree)

contains

  !> Set refinement flags for box
  subroutine refine_routine(box, cell_flags)
    type(box3_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(box%n_cell, box%n_cell, box%n_cell)
    real(dp)                 :: diff
    integer                  :: i, j, k, nc

    nc   = box%n_cell

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

             cell_flags(i, j, k) = af_keep_ref
             if (box%lvl < 2 .or. diff > 0.5e-3_dp .and. box%lvl < 6) then
                cell_flags(i, j, k) = af_do_ref
             else if (box%lvl > 4 .and. diff < 0.1_dp * 0.5e-3) then
                cell_flags(i, j, k) = af_rm_ref
             end if
          end do
       end do
    end do
  end subroutine refine_routine

  !> This routine sets the initial conditions for each box
  subroutine set_initial_condition(box)
    type(box3_t), intent(inout) :: box
    integer                     :: i, j, k, nc
    real(dp)                    :: xyz(3)

    nc = box%n_cell
    do k = 0, nc+1
       do j = 0, nc+1
          do i = 0, nc+1
             xyz = a3_r_cc(box, [i,j,k])
             box%cc(i, j, k, i_phi) = solution(xyz, 0.0_dp)
          end do
       end do
    end do
  end subroutine set_initial_condition

  !> This routine computes the error in i_phi
  subroutine set_error(box, time)
    type(box3_t), intent(inout) :: box
    real(dp), intent(in)        :: time(:)
    integer                     :: i, j, k, nc
    real(dp)                    :: xyz(3)

    nc = box%n_cell
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             xyz = a3_r_cc(box, [i,j,k])
             box%cc(i, j, k, i_err) = &
                  box%cc(i, j, k, i_phi) - solution(xyz, time(1))
          end do
       end do
    end do
  end subroutine set_error

  !> This routine calculates the analytic solution in point xyz
  function solution(xyz, t) result(sol)
    real(dp), intent(in) :: xyz(3), t
    real(dp)             :: sol, xyz_t(3)
    integer, parameter   :: sol_type = 1

    xyz_t = xyz - [vel_x, vel_y, vel_z] * t

    select case (sol_type)
    case (1)
       sol = sin(xyz_t(1))**4 * cos(xyz_t(3))**4 ! * cos(xyz_t(3))**4
    case (2)
       xyz_t = modulo(xyz_t, domain_len) / domain_len
       if (norm2(xyz_t - 0.5_dp) < 0.1_dp) then
          sol = 1.0_dp
       else
          sol = 0.0_dp
       end if
    end select
  end function solution

!not used!  !> This routine computes the drift diffusion upwind fluxes for all boxes
!not used!  subroutine fluxes_upwind1(boxes, id, flux_args)
!not used!    type(box3_t), intent(inout) :: boxes(:)
!not used!    integer, intent(in)         :: id
!not used!    real(dp), intent(in)        :: flux_args(:)
!not used!    real(dp)                    :: inv_dr
!not used!    integer                     :: nc
!not used!
!not used!    nc     = boxes(id)%n_cell
!not used!    inv_dr = 1/boxes(id)%dr
!not used!
!not used!    ! Diffusion
!not used!    boxes(id)%fx(:,:,:,i_phi) = boxes(id)%cc(0:nc, 1:nc, 1:nc, i_phi) &
!not used!         - boxes(id)%cc(1:nc+1, 1:nc, 1:nc, i_phi)
!not used!    boxes(id)%fx(:,:,:,i_phi) = boxes(id)%fx(:,:,:,i_phi) * flux_args(1) * inv_dr
!not used!
!not used!    boxes(id)%fy(:,:,:,i_phi) = boxes(id)%cc(1:nc, 0:nc, 1:nc, 1) &
!not used!         - boxes(id)%cc(1:nc, 1:nc+1, 1:nc, 1)
!not used!    boxes(id)%fy(:,:,:,i_phi) = boxes(id)%fy(:,:,:,i_phi) * flux_args(1) * inv_dr
!not used!
!not used!    boxes(id)%fz(:,:,:,i_phi) = boxes(id)%cc(1:nc, 1:nc, 0:nc, 1) &
!not used!         - boxes(id)%cc(1:nc, 1:nc, 1:nc+1, 1)
!not used!    boxes(id)%fz(:,:,:,i_phi) = boxes(id)%fz(:,:,:,i_phi) * flux_args(1) * inv_dr
!not used!
!not used!    ! Drift (1st order upwind, which is very diffusive!)
!not used!    if (flux_args(2) > 0) then
!not used!       boxes(id)%fx(:,:,:,i_phi) = boxes(id)%fx(:,:,:,i_phi) + &
!not used!            flux_args(2) * boxes(id)%cc(0:nc, 1:nc, 1:nc, 1)
!not used!    else
!not used!       boxes(id)%fx(:,:,:,i_phi) = boxes(id)%fx(:,:,:,i_phi) + &
!not used!            flux_args(2) * boxes(id)%cc(1:nc+1, 1:nc, 1:nc, 1)
!not used!    end if
!not used!    if (flux_args(3) > 0) then
!not used!       boxes(id)%fy(:,:,:,i_phi) = boxes(id)%fy(:,:,:,i_phi) + &
!not used!            flux_args(3) * boxes(id)%cc(1:nc, 0:nc, 1:nc, 1)
!not used!    else
!not used!       boxes(id)%fy(:,:,:,i_phi) = boxes(id)%fy(:,:,:,i_phi) + &
!not used!            flux_args(3) * boxes(id)%cc(1:nc, 1:nc+1, 1:nc, 1)
!not used!    end if
!not used!    if (flux_args(4) > 0) then
!not used!       boxes(id)%fz(:,:,:,i_phi) = boxes(id)%fz(:,:,:,i_phi) + &
!not used!            flux_args(4) * boxes(id)%cc(1:nc, 1:nc, 0:nc, 1)
!not used!    else
!not used!       boxes(id)%fz(:,:,:,i_phi) = boxes(id)%fz(:,:,:,i_phi) + &
!not used!            flux_args(4) * boxes(id)%cc(1:nc, 1:nc, 1:nc+1, 1)
!not used!    end if
!not used!  end subroutine fluxes_upwind1

  !> This routine computes the central differences in x-, y- and z-direction
  subroutine fluxes_centdif(boxes, id, flux_args)
    type(box3_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp), intent(in)        :: flux_args(:)
    real(dp)                    :: inv_dr
    integer                     :: nc

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    ! Diffusion
    boxes(id)%fx(:,:,:,i_phi) = boxes(id)%cc(0:nc, 1:nc, 1:nc, i_phi) &
         - boxes(id)%cc(1:nc+1, 1:nc, 1:nc, i_phi)
    boxes(id)%fx(:,:,:,i_phi) = boxes(id)%fx(:,:,:,i_phi) * flux_args(1) * inv_dr

    boxes(id)%fy(:,:,:,i_phi) = boxes(id)%cc(1:nc, 0:nc, 1:nc, 1) &
         - boxes(id)%cc(1:nc, 1:nc+1, 1:nc, 1)
    boxes(id)%fy(:,:,:,i_phi) = boxes(id)%fy(:,:,:,i_phi) * flux_args(1) * inv_dr

    boxes(id)%fz(:,:,:,i_phi) = boxes(id)%cc(1:nc, 1:nc, 0:nc, 1) &
         - boxes(id)%cc(1:nc, 1:nc, 1:nc+1, 1)
    boxes(id)%fz(:,:,:,i_phi) = boxes(id)%fz(:,:,:,i_phi) * flux_args(1) * inv_dr

    ! Drift (centered differences)
    boxes(id)%fx(:,:,:,i_phi) = boxes(id)%fx(:,:,:,i_phi) + &
         flux_args(2) * 0.5_dp * &
         (boxes(id)%cc(0:nc, 1:nc, 1:nc, 1) + boxes(id)%cc(1:nc+1, 1:nc, 1:nc, 1))
    boxes(id)%fy(:,:,:,i_phi) = boxes(id)%fy(:,:,:,i_phi) + &
         flux_args(3) * 0.5_dp * &
         (boxes(id)%cc(1:nc, 0:nc, 1:nc, 1) + boxes(id)%cc(1:nc, 1:nc+1, 1:nc, 1))
    boxes(id)%fz(:,:,:,i_phi) = boxes(id)%fz(:,:,:,i_phi) + &
         flux_args(4) * 0.5_dp * &
         (boxes(id)%cc(1:nc, 1:nc, 0:nc, 1) + boxes(id)%cc(1:nc, 1:nc, 1:nc+1, 1))
  end subroutine fluxes_centdif

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
  subroutine fluxes_koren(boxes, id, flux_args)
    use m_a3_prolong
    type(box3_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp), intent(in)        :: flux_args(:)
    real(dp)                    :: tmp, inv_dr
    real(dp)                    :: gradp, gradc, gradn
    real(dp)                    :: gc_data(boxes(id)%n_cell, &
         boxes(id)%n_cell, a3_num_neighbors)
    integer                     :: i, j, k, nc

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    ! Get a second layer of ghost cell data (the 'normal' routines give just one
    ! layer of ghost cells). Use a3_gc2_prolong_linear on refinement boundaries and
    ! a3_bc2_neumann_zero on physical boundaries.
    call a3_gc2_box(boxes, &                 ! List of all the boxes
                    id, &                    ! Id of box for which we set ghost cells
                    i_phi, &                 ! Variable for which ghost cells are set
                    a3_gc2_prolong_linear, & ! Procedure called at refinement boundaries
                    a3_bc2_neumann_zero, &   ! Procedure called at physical boundaries
                    gc_data, &               ! The requested ghost cells
                    nc)                      ! box%n_cell

    ! x-fluxes interior, advective part with flux limiter
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc+1
             gradc = boxes(id)%cc(i, j, k, i_phi) - boxes(id)%cc(i-1, j, k, i_phi)
             if (flux_args(2) < 0.0_dp) then
                if (i == nc+1) then
                   tmp = gc_data(j, k, a3_neighb_highx)
                else
                   tmp = boxes(id)%cc(i+1, j, k, i_phi)
                end if

                gradn = tmp - boxes(id)%cc(i, j, k, i_phi)
                boxes(id)%fx(i, j, k, i_phi) = flux_args(2) * &
                     (boxes(id)%cc(i, j, k, i_phi) - koren_mlim(gradc, gradn))
             else                  ! flux_args(2) > 0
                if (i == 1) then
                   tmp = gc_data(j, k, a3_neighb_lowx)
                else
                   tmp = boxes(id)%cc(i-2, j, k, i_phi)
                end if

                gradp = boxes(id)%cc(i-1, j, k, i_phi) - tmp
                boxes(id)%fx(i, j, k, i_phi) = flux_args(2) * &
                     (boxes(id)%cc(i-1, j, k, i_phi) + koren_mlim(gradc, gradp))
             end if

             ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
             boxes(id)%fx(i, j, k, i_phi) = boxes(id)%fx(i, j, k, i_phi) - &
                  flux_args(1) * gradc * inv_dr
          end do
       end do
    end do

    ! y-fluxes interior, advective part with flux limiter
    do k = 1, nc
       do j = 1, nc+1
          do i = 1, nc
             gradc = boxes(id)%cc(i, j, k, i_phi) - boxes(id)%cc(i, j-1, k, i_phi)
             if (flux_args(3) < 0.0_dp) then
                if (j == nc+1) then
                   tmp = gc_data(i, k, a3_neighb_highy)
                else
                   tmp = boxes(id)%cc(i, j+1, k, i_phi)
                end if

                gradn = tmp - boxes(id)%cc(i, j, k, i_phi)
                boxes(id)%fy(i, j, k, i_phi) = flux_args(3) * &
                     (boxes(id)%cc(i, j, k, i_phi) - koren_mlim(gradc, gradn))
             else                  ! flux_args(3) > 0

                if (j == 1) then
                   tmp = gc_data(i, k, a3_neighb_lowy)
                else
                   tmp = boxes(id)%cc(i, j-2, k, i_phi)
                end if

                gradp = boxes(id)%cc(i, j-1, k, i_phi) - tmp
                boxes(id)%fy(i, j, k, i_phi) = flux_args(3) * &
                     (boxes(id)%cc(i, j-1, k, i_phi) + koren_mlim(gradc, gradp))
             end if

             ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
             boxes(id)%fy(i, j, k, i_phi) = boxes(id)%fy(i, j, k, i_phi) - &
                  flux_args(1) * gradc * inv_dr
          end do
       end do
    end do

    ! z-fluxes interior, advective part with flux limiter
    do k = 1, nc+1
       do j = 1, nc
          do i = 1, nc
             gradc = boxes(id)%cc(i, j, k, i_phi) - boxes(id)%cc(i, j, k-1, i_phi)
             if (flux_args(4) < 0.0_dp) then
                if (k == nc+1) then
                   tmp = gc_data(i, j, a3_neighb_highz)
                else
                   tmp = boxes(id)%cc(i, j, k+1, i_phi)
                end if

                gradn = tmp - boxes(id)%cc(i, j, k, i_phi)
                boxes(id)%fz(i, j, k, i_phi) = flux_args(4) * &
                     (boxes(id)%cc(i, j, k, i_phi) - koren_mlim(gradc, gradn))
             else                  ! flux_args(4) > 0

                if (k == 1) then
                   tmp = gc_data(i, j, a3_neighb_lowz)
                else
                   tmp = boxes(id)%cc(i, j, k-2, i_phi)
                end if

                gradp = boxes(id)%cc(i, j, k-1, i_phi) - tmp
                boxes(id)%fz(i, j, k, i_phi) = flux_args(4) * &
                     (boxes(id)%cc(i, j, k-1, i_phi) + koren_mlim(gradc, gradp))
             end if

             ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
             boxes(id)%fz(i, j, k, i_phi) = boxes(id)%fz(i, j, k, i_phi) - &
                  flux_args(1) * gradc * inv_dr
          end do
       end do
    end do
  end subroutine fluxes_koren

  !> This routine computes the update of the solution per box
  subroutine update_solution(box, dt)
    type(box3_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr
    integer                     :: nc

    nc                    = box%n_cell
    inv_dr                = 1/box%dr
    box%cc(1:nc, 1:nc, 1:nc, i_phi) = box%cc(1:nc, 1:nc, 1:nc, i_phi) + dt(1) * ( &
         (box%fx(1:nc, :, :, i_phi) - box%fx(2:nc+1, :, :, i_phi)) * inv_dr + &
         (box%fy(:, 1:nc, :, i_phi) - box%fy(:, 2:nc+1, :, i_phi)) * inv_dr + &
         (box%fz(:, :, 1:nc, i_phi) - box%fz(:, :, 2:nc+1, i_phi)) * inv_dr)
  end subroutine update_solution

  !> This routine computes the update of the solution per box
  subroutine average_phi(box)
    type(box3_t), intent(inout) :: box
    box%cc(:, :, :, i_phi) = 0.5_dp * (box%cc(:, :, :, i_phi) + box%cc(:, :, :, i_phi_old))
  end subroutine average_phi

  ! Linear prolongation of i_phi values to new children
  subroutine prolong_to_new_children(tree, refine_info)
    use m_a3_prolong
    type(a3_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: refine_info
    integer                      :: lvl, i, id, p_id

    do lvl = 1, tree%highest_lvl
       do i = 1, size(refine_info%lvls(lvl)%add)
          id = refine_info%lvls(lvl)%add(i)
          p_id = tree%boxes(id)%parent

          ! Linear prolongation will not strictly conserve phi
          call a3_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_phi)
       end do

       do i = 1, size(refine_info%lvls(lvl)%add)
          id = refine_info%lvls(lvl)%add(i)
          ! Fill ghost cells for variables i_phi on the sides of a box, using
          ! a3_gc_interp_lim on refinement boundaries and a3_bc_neumann_zero on physical boundaries
          call a3_gc_box(tree%boxes, & ! List of all the boxes 
                         id, &         ! Id of box for which we set ghost cells
                         i_phi, &      ! Variable for which ghost cells are set
                         a3_gc_interp_lim, & ! Procedure called at refinement boundaries
                         a3_bc_neumann_zero) ! Procedure called at physical boundaries
       end do
    end do
  end subroutine prolong_to_new_children

end program drift_diffusion_3d
