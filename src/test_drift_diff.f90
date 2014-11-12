program test_drift_diff
  use omp_lib
  use m_afivo_2d

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  type(a2_t)         :: tree
  integer            :: i, n, n_steps, id, lvl
  integer            :: ix_list(2, 1)
  integer            :: nb_list(4, 1)
  integer, parameter :: box_size    = 8
  integer, parameter :: i_phi       = 1
  integer, parameter :: i_phi_old = 2
  integer            :: n_boxes_max = 10*1000
  integer            :: n_lvls_max  = 20
  integer            :: n_changes
  real(dp)           :: dr, dt, time, end_time
  real(dp)           :: time_per_adapt, time_in_loop
  real(dp)           :: diff_coeff, vel_x, vel_y, dr_min(2)
  character(len=40)  :: fname
  logical            :: done_with_loop
  logical :: forward_euler = .false.

  dr = 4.0_dp / box_size

  ! Initialize tree
  call a2_init(tree, n_lvls_max, n_boxes_max, box_size, n_var_cell=2, &
       n_var_face=1, dr = dr, r_min = [0.0_dp, 0.0_dp])

  id             = 1
  ix_list(:, id) = [1,1] ! Set index of box
  nb_list(:, id) = id    ! Box is periodic, so its own neighbor

  ! Set up the initial conditions. First adapt mesh.
  call a2_set_base(tree, ix_list, nb_list)
  do i = 1, 20
     call a2_adjust_refinement(tree, ref_func_init, n_changes)
     if (n_changes == 0) exit
  end do

  call a2_loop_box(tree, set_init_cond)
  call a2_gc_sides(tree, [i_phi], a2_sides_prolong1, have_no_bc)
  call a2_gc_corners(tree, [i_phi], a2_corners_prolong1, have_no_bc)

  do i = 1, 20
     call a2_adjust_refinement(tree, ref_func, n_changes)
     call a2_loop_boxes(tree, prolong_to_new_children)

     if (n_changes == 0) exit
  end do

  ! Set final initial condition
  call a2_loop_box(tree, set_init_cond)
  call a2_restrict_tree(tree, [i_phi])
  call a2_gc_sides(tree, [i_phi], a2_sides_extrap, have_no_bc)
  call a2_gc_corners(tree, [i_phi], a2_corners_extrap, have_no_bc)

  i              = 0
  time           = 0
  time_per_adapt = 0.02_dp
  end_time       = 1.0_dp
  diff_coeff     = 0.0_dp
  vel_x          = 1.0_dp
  vel_y          = 0.0_dp

  !$omp parallel private(n)
  do
     !$omp single
     i = i + 1
     print *, "i = ", i, "n_boxes", tree%max_id, time
     write(fname, "(A,I0,A)") "test_drift_diff_", i, ".vtu"

     dr_min  = a2_min_dr(tree)
     dt      = 0.5_dp / (2 * diff_coeff * sum(1/dr_min**2) + &
          sum( abs([vel_x, vel_y]) / dr_min ) + epsilon(1.0_dp))
     n_steps = ceiling(time_per_adapt/dt)
     dt      = time_per_adapt / n_steps
     time    = time + time_per_adapt
     !$omp end single

     call a2_write_tree(tree, trim(fname), (/"phi", "old"/), i, time)

     !$omp barrier
     if (time > end_time) exit

     if (forward_euler) then
        do n = 1, n_steps
           call a2_loop_boxes_arg(tree, calculate_fluxes_limiter, [diff_coeff, vel_x, vel_y])
           call a2_loop_box_arg(tree, update_solution, [dt])
           call a2_restrict_tree(tree, [i_phi])
           call a2_gc_sides(tree, [i_phi], a2_sides_extrap, have_no_bc)
        end do
     else                    ! Modified midpoint
        do n = 1, n_steps
           ! Copy previous solution
           call a2_tree_copy_cc(tree, i_phi, i_phi_old)

           ! Take a half time step
           call a2_loop_boxes_arg(tree, calculate_fluxes_limiter, [diff_coeff, vel_x, vel_y])
           call a2_loop_box_arg(tree, update_solution, [0.5_dp * dt])
           call a2_restrict_tree(tree, [i_phi])
           call a2_gc_sides(tree, [i_phi], a2_sides_extrap, have_no_bc)

           ! Calculate fluxes
           call a2_loop_boxes_arg(tree, calculate_fluxes_limiter, [diff_coeff, vel_x, vel_y])

           ! Copy back old phi, and take full time step
           call a2_tree_copy_cc(tree, i_phi_old, i_phi)
           call a2_loop_box_arg(tree, update_solution, [dt])
           call a2_restrict_tree(tree, [i_phi])
           call a2_gc_sides(tree, [i_phi], a2_sides_extrap, have_no_bc)
        end do
     end if

     call a2_restrict_tree(tree, [i_phi])
     call a2_gc_sides(tree, [i_phi], a2_sides_extrap, have_no_bc)
     call a2_gc_corners(tree, [i_phi], a2_corners_extrap, have_no_bc)

     call a2_adjust_refinement(tree, ref_func)
     ! call a2_tidy_up(tree, 0.5_dp, 0.25_dp, 100*1000, .false.)
     call a2_loop_boxes(tree, prolong_to_new_children)
  end do
  !$omp end parallel

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

    nc   = box%n_cell
    diff = 0
    do j = 1, nc
       do i = 1, nc
          diff = max(diff, (box%cc(i+1, j, i_phi) - box%cc(i, j, i_phi))**2 + &
               (box%cc(i, j+1, i_phi) - box%cc(i, j, i_phi))**2 + &
               (box%cc(i, j, i_phi) - box%cc(i-1, j, i_phi))**2 + &
               (box%cc(i, j, i_phi) - box%cc(i, j-1, i_phi))**2)
       end do
    end do
    diff = sqrt(0.5_dp * diff)

    if (diff > 0.05_dp .and. box%lvl < 6) then
       ref_func = a5_do_ref
    else
       ref_func = a5_rm_ref
    end if
  end function ref_func

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    box%cc(i, j, i_phi) = 0

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          if (norm2(xy - 2) < 1) then
             box%cc(i, j, i_phi) = 1
          else if (norm2(xy - 2) < 1.1_dp) then
             box%cc(i, j, i_phi) = (1.1_dp - norm2(xy - 2)) * 10
          else
             box%cc(i, j, i_phi) = 0
          end if
       end do
    end do
  end subroutine set_init_cond

  subroutine calculate_fluxes(box, flux_args)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: flux_args(:)
    real(dp)                    :: inv_dr
    integer                     :: nc

    nc     = box%n_cell
    inv_dr = 1/box%dr

    ! Diffusion
    box%fx(:,:,i_phi) = box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi)
    box%fx(:,:,i_phi) = box%fx(:,:,i_phi) * flux_args(1) * inv_dr
    box%fy(:,:,i_phi) = box%cc(1:nc, 0:nc, 1) - box%cc(1:nc, 1:nc+1, 1)
    box%fy(:,:,i_phi) = box%fy(:,:,i_phi) * flux_args(1) * inv_dr

    ! Drift (1st order upwind, which is very diffusive!)
    if (flux_args(2) > 0) then
       box%fx(:,:,i_phi) = box%fx(:,:,i_phi) + &
            flux_args(2) * box%cc(0:nc, 1:nc, 1)
    else
       box%fx(:,:,i_phi) = box%fx(:,:,i_phi) + &
            flux_args(2) * box%cc(1:nc+1, 1:nc, 1)
    end if
    if (flux_args(3) > 0) then
       box%fy(:,:,i_phi) = box%fy(:,:,i_phi) + &
            flux_args(3) * box%cc(1:nc, 0:nc, 1)
    else
       box%fy(:,:,i_phi) = box%fy(:,:,i_phi) + &
            flux_args(3) * box%cc(1:nc, 1:nc+1, 1)
    end if
  end subroutine calculate_fluxes

  elemental function limiter_koren(theta)
    real(dp), intent(in) :: theta
    real(dp)             :: limiter_koren
    real(dp), parameter  :: one_sixth = 1.0_dp / 6.0_dp
    limiter_koren = max(0.0d0, min(1.0_dp, theta, (1.0_dp + 2.0_dp * theta) * one_sixth))
  end function limiter_koren

  elemental function limiter_grad_ratio(numerator, denominator)
    real(dp), intent(in) :: numerator, denominator
    real(dp)             :: limiter_grad_ratio
    real(dp), parameter  :: eps = epsilon(1.0d0)
    ! Avoid division by zero, and ensure that at zero gradients we have a ratio of 1
    limiter_grad_ratio = (sign(eps, numerator) + numerator) / &
         (denominator + sign(eps, denominator))
  end function limiter_grad_ratio

  subroutine calculate_fluxes_limiter(boxes, id, flux_args)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp), intent(in)        :: flux_args(:)
    real(dp)                    :: inv_dr, theta
    real(dp)                    :: gradp, gradc, gradn
    integer                     :: i, j, nc, nb_id

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    ! x-fluxes interior
    do j = 1, nc
       do i = 2, nc
          ! Advective part with flux limiter
          if (flux_args(2) < 0.0_dp) then
             gradc = boxes(id)%cc(i, j, i_phi) - boxes(id)%cc(i-1, j, i_phi)
             gradn = boxes(id)%cc(i+1, j, i_phi) - boxes(id)%cc(i, j, i_phi)
             theta     = limiter_grad_ratio(gradc, gradn)
             boxes(id)%fx(i, j, i_phi) = flux_args(2) * &
                  (boxes(id)%cc(i, j, i_phi) - limiter_koren(theta) * gradn)
          else
             gradp = boxes(id)%cc(i-1, j, i_phi) - boxes(id)%cc(i-2, j, i_phi)
             gradc = boxes(id)%cc(i, j, i_phi) - boxes(id)%cc(i-1, j, i_phi)
             theta     = limiter_grad_ratio(gradc, gradp)
             boxes(id)%fx(i, j, i_phi) = flux_args(2) * &
                  (boxes(id)%cc(i-1, j, i_phi) + limiter_koren(theta) * gradp)
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
          boxes(id)%fx(i, j, i_phi) = boxes(id)%fx(i, j, i_phi) - flux_args(1) * gradc
       end do
    end do

    ! x-lo side
    i = 1
    nb_id = boxes(id)%neighbors(a2_nb_lx)

    do j = 1, nc
       if (flux_args(2) < 0.0_dp) then ! flux towards boundary
          gradc = boxes(id)%cc(i, j, i_phi) - boxes(id)%cc(i-1, j, i_phi)
          gradn = boxes(id)%cc(i+1, j, i_phi) - boxes(id)%cc(i, j, i_phi)

          theta     = limiter_grad_ratio(gradc, gradn)
          boxes(id)%fx(i, j, i_phi) = flux_args(2) * &
               (boxes(id)%cc(i, j, i_phi) - limiter_koren(theta) * gradn)
       else if (nb_id > a5_no_box) then
          gradp = boxes(id)%cc(i-1, j, i_phi) - boxes(nb_id)%cc(nc-1, j, i_phi)
          gradc = boxes(id)%cc(i, j, i_phi) - boxes(id)%cc(i-1, j, i_phi)

          theta     = limiter_grad_ratio(gradc, gradp)
          boxes(id)%fx(i, j, i_phi) = flux_args(2) * &
               (boxes(id)%cc(i-1, j, i_phi) + limiter_koren(theta) * gradp)
       else
          ! First order upwind
          boxes(id)%fx(i, j, i_phi) = flux_args(2) * boxes(id)%cc(i-1, j, i_phi)
       end if

       boxes(id)%fx(i, j, i_phi) = boxes(id)%fx(i, j, i_phi) - flux_args(1) * gradc
    end do

    ! x-hi side
    i = nc+1
    nb_id = boxes(id)%neighbors(a2_nb_hx)

    do j = 1, nc
       if (flux_args(2) > 0.0_dp) then ! flux towards boundary
          gradp = boxes(id)%cc(i-1, j, i_phi) - boxes(id)%cc(i-2, j, i_phi)
          gradc = boxes(id)%cc(i, j, i_phi) - boxes(id)%cc(i-1, j, i_phi)

          theta     = limiter_grad_ratio(gradc, gradp)
          boxes(id)%fx(i, j, i_phi) = flux_args(2) * (boxes(id)%cc(i-1, j, i_phi) + limiter_koren(theta) * gradp)
       else if (nb_id > a5_no_box) then
          gradn = boxes(nb_id)%cc(2, j, i_phi) - boxes(id)%cc(i, j, i_phi)
          gradc = boxes(id)%cc(i, j, i_phi) - boxes(id)%cc(i-1, j, i_phi)

          theta     = limiter_grad_ratio(gradc, gradn)
          boxes(id)%fx(i, j, i_phi) = flux_args(2) * &
               (boxes(id)%cc(i, j, i_phi) - limiter_koren(theta) * gradn)
       else
          ! First order upwind
          boxes(id)%fx(i, j, i_phi) = flux_args(2) * boxes(id)%cc(i, j, i_phi)
       end if

       boxes(id)%fx(i, j, i_phi) = boxes(id)%fx(i, j, i_phi) - flux_args(1) * gradc
    end do

    boxes(id)%fy = 0

    ! ! y-fluxes
    ! do i = 1, nc
    !    do j = 2, nc

    !       gradp = box%cc(i, j-1, i_phi) - box%cc(i, j-2, i_phi)
    !       gradc = box%cc(i, j, i_phi) - box%cc(i, j-1, i_phi)
    !       gradn = box%cc(i, j+1, i_phi) - box%cc(i, j, i_phi)

    !       ! Advective part with flux limiter
    !       if (flux_args(3) < 0.0_dp) then
    !          theta     = limiter_grad_ratio(gradc, gradn)
    !          box%fy(i, j, i_phi) = flux_args(3) * (box%cc(i, j, i_phi) - limiter_koren(theta) * gradn)
    !       else
    !          theta     = limiter_grad_ratio(gradc, gradp)
    !          box%fy(i, j, i_phi) = flux_args(3) * (box%cc(i, j-1, i_phi) + limiter_koren(theta) * gradp)
    !       end if

    !       ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
    !       box%fy(i, j, i_phi) = box%fy(i, j, i_phi) - flux_args(1) * gradc
    !    end do
    ! end do

  end subroutine calculate_fluxes_limiter

  subroutine update_solution(box, dt)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr
    integer                     :: nc

    nc                    = box%n_cell
    inv_dr                = 1/box%dr
    box%cc(1:nc, 1:nc, 1) = box%cc(1:nc, 1:nc, 1) + dt(1) * ( &
         (box%fx(1:nc, :, 1) - box%fx(2:nc+1, :, 1)) * inv_dr + &
         (box%fy(:, 1:nc, 1) - box%fy(:, 2:nc+1, 1)) * inv_dr)
  end subroutine update_solution

  subroutine prolong_to_new_children(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in) :: id

    if (btest(boxes(id)%tag, a5_bit_new_children)) then
       call a2_prolong1_from(boxes, id, [1], .true.)
       boxes(id)%tag = ibclr(boxes(id)%tag, a5_bit_new_children)
    end if
  end subroutine prolong_to_new_children

  subroutine have_no_bc(boxes, id, i, ivs)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, i, ivs(:)
    stop "We have no boundary conditions in this example"
  end subroutine have_no_bc

end program test_drift_diff
