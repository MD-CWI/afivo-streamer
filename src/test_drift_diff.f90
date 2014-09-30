program test_drift_diff
  use m_afivo

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  type(a2_t)         :: tree
  integer            :: i, id
  integer            :: ix_list(2, 1)
  integer            :: nb_list(4, 1)
  integer, parameter :: box_size    = 8
  integer            :: n_boxes_max = 10*1000
  real(dp)           :: dr, dt, time, end_time
  real(dp)           :: time_per_adapt, time_in_loop
  real(dp)           :: diff_coeff, vel_x, vel_y, dr_min(2)
  character(len=40)  :: fname
  logical            :: done_with_loop

  dr = 4.0_dp / box_size

  ! Initialize tree
  call a2_init(tree, n_boxes_max, box_size, n_var_cell=1, n_var_face=1, &
       dr = dr, r_min = [0.0_dp, 0.0_dp])

  id = 1
  ix_list(:, id) = [1,1]         ! Set index of box
  nb_list(:, id) = id            ! Box is periodic, so its own neighbor

  ! Set up the initial conditions
  call a2_set_base(tree, ix_list, nb_list)
  do i = 1, 3
     call a2_adjust_refinement(tree, ref_func_init)
  end do

  call a2_loop_box(tree, set_init_cond)
  call a2_gc_sides(tree, [1], a2_sides_prolong1, have_no_bc)
  call a2_gc_corners(tree, [1], a2_corners_prolong1, have_no_bc)

  do i = 1, 5
     call a2_adjust_refinement(tree, ref_func)
     call a2_loop_boxes(tree, prolong_to_new_children)
  end do

  call a2_loop_box(tree, set_init_cond)
  call a2_gc_sides(tree, [1], a2_sides_prolong1, have_no_bc)
  call a2_gc_corners(tree, [1], a2_corners_prolong1, have_no_bc)

  i              = 0
  time           = 0
  time_per_adapt = 0.02_dp
  end_time       = 1.0_dp
  diff_coeff     = 0.1_dp
  vel_x          = 2.0_dp
  vel_y          = 2.0_dp

  do while (time < end_time)
     i = i + 1
     print *, "i = ", i, "n_boxes", tree%n_boxes

     write(fname, "(A,I0,A)") "test_drift_diff_", i, ".vtu"
     call a2_write_tree(tree, trim(fname), (/"my_var"/), i, time)

     ! Advance time_per_adapt
     done_with_loop = .false.
     time_in_loop   = 0

     do while (.not. done_with_loop)
        ! Set diffusion and CFL limit for timestep
        dr_min = a2_min_dr(tree)
        dt = 0.9_dp / (2 * diff_coeff * sum(1/dr_min**2) + &
             sum( abs([vel_x, vel_y]) / dr_min ) + epsilon(1.0_dp))

        if (time_in_loop + dt > time_per_adapt) then
           dt = time_per_adapt - time_in_loop
           done_with_loop = .true.
        end if

        call a2_loop_box_arg(tree, calculate_fluxes, [diff_coeff, vel_x, vel_y])
        call a2_consistent_fluxes(tree, [1])
        call a2_loop_box_arg(tree, update_solution, [dt])
        call a2_gc_sides(tree, [1], a2_sides_prolong1, have_no_bc)
        call a2_gc_corners(tree, [1], a2_corners_prolong1, have_no_bc)
        time_in_loop = time_in_loop + dt
     end do

     time = time + time_per_adapt

     call a2_loop_boxes(tree, restrict_from_children)
     call a2_adjust_refinement(tree, ref_func)
     call a2_loop_boxes(tree, prolong_to_new_children)
     call a2_tidy_up(tree, 0.5_dp, 0.25_dp, 100*1000, .false.)
  end do

  call a2_destroy(tree)

contains

  integer function ref_func_init(box)
    type(box2_t), intent(in) :: box
    if (box%lvl < 5) then
       ref_func_init = a5_do_ref
    else
       ref_func_init = a5_rm_ref
    end if
  end function ref_func_init

  integer function ref_func(box)
    type(box2_t), intent(in) :: box
    real(dp)                 :: diff
    integer                  :: i, j, nc

    nc   = box%cfg%n_cell
    diff = 0
    do j = 1, nc
       do i = 1, nc
          diff = max(diff, (box%cc(i+1, j, 1) - box%cc(i, j, 1))**2 + &
               (box%cc(i, j+1, 1) - box%cc(i, j, 1))**2 + &
               (box%cc(i, j, 1) - box%cc(i-1, j, 1))**2 + &
               (box%cc(i, j, 1) - box%cc(i, j-1, 1))**2)
       end do
    end do
    diff = sqrt(0.5_dp * diff)

    if (diff > 0.08_dp) then
       ref_func = a5_do_ref
    else
       ref_func = a5_rm_ref
    end if
  end function ref_func

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%cfg%n_cell
    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          if (norm2(xy - 2) < 1) then
             box%cc(i, j, 1) = 1
          else if (norm2(xy - 2) < 1.1_dp) then
             box%cc(i, j, 1) = (1.1_dp - norm2(xy - 2)) * 10
          else
             box%cc(i, j, 1) = 0
          end if
       end do
    end do
  end subroutine set_init_cond

  subroutine calculate_fluxes(box, flux_args)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in) :: flux_args(:)
    real(dp) :: inv_dr(2)
    integer :: nc

    nc     = box%cfg%n_cell
    inv_dr = 1/a2_dr(box)

    ! Diffusion
    box%fx(:,:,1) = box%cc(0:nc, 1:nc, 1) - box%cc(1:nc+1, 1:nc, 1)
    box%fx(:,:,1) = box%fx(:,:,1) * flux_args(1) * inv_dr(1)
    box%fy(:,:,1) = box%cc(1:nc, 0:nc, 1) - box%cc(1:nc, 1:nc+1, 1)
    box%fy(:,:,1) = box%fy(:,:,1) * flux_args(1) * inv_dr(2)

    ! Drift (1st order upwind)
    if (flux_args(2) > 0) then
       box%fx(:,:,1) = box%fx(:,:,1) + flux_args(2) * box%cc(0:nc, 1:nc, 1)
    else
       box%fx(:,:,1) = box%fx(:,:,1) + flux_args(2) * box%cc(1:nc+1, 1:nc, 1)
    end if
    if (flux_args(3) > 0) then
       box%fy(:,:,1) = box%fy(:,:,1) + flux_args(3) * box%cc(1:nc, 0:nc, 1)
    else
       box%fy(:,:,1) = box%fy(:,:,1) + flux_args(3) * box%cc(1:nc, 1:nc+1, 1)
    end if
  end subroutine calculate_fluxes

  subroutine update_solution(box, dt)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr(2)
    integer                     :: nc

    nc = box%cfg%n_cell
    inv_dr = 1/a2_dr(box)

    nc = box%cfg%n_cell
    box%cc(1:nc, 1:nc, 1) = box%cc(1:nc, 1:nc, 1) + dt(1) * ( &
         (box%fx(1:nc, :, 1) - box%fx(2:nc+1, :, 1)) * inv_dr(1) + &
         (box%fy(:, 1:nc, 1) - box%fy(:, 2:nc+1, 1)) * inv_dr(2))
  end subroutine update_solution

  subroutine prolong_to_new_children(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in) :: id

    if (btest(boxes(id)%tag, a5_bit_new_children)) then
       call a2_prolong1_from(boxes, id, [1], .true.)
       boxes(id)%tag = ibclr(boxes(id)%tag, a5_bit_new_children)
    end if
  end subroutine prolong_to_new_children

  subroutine restrict_from_children(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id

    if (a2_has_children(boxes(id))) then
       call a2_restrict_to(boxes, id, [1])
    end if
  end subroutine restrict_from_children

  subroutine have_no_bc(boxes, id, i, ivs)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, i, ivs(:)
    stop "We have no boundary conditions in this example"
  end subroutine have_no_bc

end program test_drift_diff
