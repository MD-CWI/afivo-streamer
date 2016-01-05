!> \example test_drift_diff_2d.f90
!> A drift-diffusion example for m_a2_t
program test_drift_diff
  use m_a2_t
  use m_a2_core
  use m_a2_io
  use m_a2_utils
  use m_a2_gc
  use m_a2_restrict

  implicit none

  integer, parameter :: box_size    = 8
  integer, parameter :: i_phi       = 1
  integer, parameter :: i_phi_old   = 2
  integer, parameter :: i_err       = 3

  real(dp), parameter :: domain_len = 2 * acos(-1.0_dp)
  real(dp), parameter :: dr = domain_len / box_size

  type(a2_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: i, id
  integer            :: ix_list(2, 1)
  integer            :: nb_list(4, 1)
  integer            :: output_cnt
  real(dp)           :: dt, time, end_time, p_err, n_err
  real(dp)           :: dt_output
  real(dp)           :: diff_coeff, vel_x, vel_y, dr_min(2)
  character(len=40)  :: fname

  logical            :: write_out
  integer            :: time_step_method = 2

  print *, "Initialize tree"
  call a2_init(tree, box_size, n_var_cell=3, n_var_face=1, dr = dr)

  ! Set up geometry
  id             = 1
  ix_list(:, id) = [1,1] ! Set index of box
  nb_list(:, id) = id    ! Box is periodic, so its own neighbor

  print *, "Create the base mesh"
  call a2_set_base(tree, ix_list, nb_list)

  output_cnt = 0
  time       = 0
  dt_output  = 0.1_dp
  end_time   = 5.0_dp
  diff_coeff = 0.0_dp
  vel_x      = 1.0_dp
  vel_y      = 1.0_dp

  print *, "Set up the initial conditions"
  do
     ! We should only set the finest level, but this also works
     call a2_loop_box(tree, set_init_cond)
     call a2_gc_tree(tree, i_phi, a2_gc_interp_lim, have_no_bc)
     call a2_adjust_refinement(tree, set_ref_flags, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  print *, "Restrict the initial conditions"
  call a2_restrict_tree(tree, i_phi)
  call a2_gc_tree(tree, i_phi, a2_gc_interp_lim, have_no_bc)

  do
     dr_min  = a2_min_dr(tree)
     dt      = 0.5_dp / (2 * diff_coeff * sum(1/dr_min**2) + &
          sum( abs([vel_x, vel_y]) / dr_min ) + epsilon(1.0_dp))

     if (time > end_time) exit

     select case (time_step_method)
     case (1)
        ! Forward Euler
        call a2_loop_boxes_arg(tree, fluxes_koren, [diff_coeff, vel_x, vel_y])
        call a2_consistent_fluxes(tree, [1])
        call a2_loop_box_arg(tree, update_solution, [dt])
        call a2_restrict_tree(tree, i_phi)
        call a2_gc_tree(tree, i_phi, a2_gc_interp_lim, have_no_bc)
        time = time + dt
     case (2)
        ! Copy previous solution
        call a2_tree_copy_cc(tree, i_phi, i_phi_old)

        ! Two forward Euler steps over dt
        do i = 1, 2
           call a2_loop_boxes_arg(tree, fluxes_koren, [diff_coeff, vel_x, vel_y])
           call a2_consistent_fluxes(tree, [1])
           call a2_loop_box_arg(tree, update_solution, [dt])
           call a2_restrict_tree(tree, i_phi)
           call a2_gc_tree(tree, i_phi, a2_gc_interp_lim, have_no_bc)
        end do

        ! Take average of phi_old and phi
        call a2_loop_box(tree, average_phi)
        time = time + dt
     end select

     if (output_cnt * dt_output <= time) then
        write_out = .true.
        output_cnt = output_cnt + 1
        write(fname, "(A,I0)") "test_drift_diff_2d_", output_cnt
     else
        write_out = .false.
     end if

     if (write_out) then
        call a2_loop_box_arg(tree, set_error, [time])
        call a2_write_silo(tree, trim(fname), &
             (/"phi", "tmp", "err"/), output_cnt, time, &
             fc_names=["fx", "fy"], ixs_fc=[1])
        call a2_tree_max_cc(tree, i_err, p_err)
        call a2_tree_min_cc(tree, i_err, n_err)
        print *, "max error", max(p_err, abs(n_err))
        call a2_tree_sum_cc(tree, i_phi, p_err)
        print *, "sum phi", p_err
     end if

     call a2_adjust_refinement(tree, set_ref_flags, ref_info)
     call prolong_to_new_children(tree, ref_info)
     call a2_gc_tree(tree, i_phi, a2_gc_interp_lim, have_no_bc)
     call a2_tidy_up(tree, 0.8_dp, 0.5_dp, 10000, .false.)
  end do

  call a2_destroy(tree)

contains

  subroutine set_ref_flags(boxes, id, ref_flags)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flags(:)
    real(dp)                 :: diff
    integer                  :: nc

    nc   = boxes(id)%n_cell
    diff = max( &
         maxval(abs(boxes(id)%cc(1:nc+1, 1:nc, i_phi) - &
         boxes(id)%cc(0:nc, 1:nc, i_phi))), &
         maxval(abs(boxes(id)%cc(1:nc, 1:nc+1, i_phi) - &
         boxes(id)%cc(1:nc, 0:nc, i_phi))))

    if (boxes(id)%lvl < 3 .or. diff > 0.05_dp) then
       ref_flags(id) = a5_do_ref
    else if (diff > 0.2_dp * 0.05_dp) then
       ref_flags(id) = a5_kp_ref
    else if (boxes(id)%lvl > 4) then
       ref_flags(id) = a5_rm_ref
    end if

    ! if (boxes(id)%lvl < 6) ref_flags(id) = a5_do_ref
    ! if (boxes(id)%lvl < 6 .and. boxes(id)%r_min(1) < 4) ref_flags(id) = a5_do_ref
  end subroutine set_ref_flags

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 0, nc+1
       do i = 0, nc+1
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, i_phi) = solution(xy, 0.0_dp)
       end do
    end do
  end subroutine set_init_cond

  subroutine set_error(box, time)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: time(:)
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          xy = a2_r_cc(box, [i,j])
          box%cc(i, j, i_err) = &
               box%cc(i, j, i_phi) - solution(xy, time(1))
       end do
    end do
  end subroutine set_error

  function solution(xy, t) result(sol)
    real(dp), intent(in) :: xy(2), t
    real(dp) :: sol, xy_t(2)

    xy_t = xy - [vel_x, vel_y] * t
    sol = sin(xy_t(1))**4 * cos(xy_t(2))**4

    ! xy_t = modulo(xy_t, domain_len) / domain_len
    ! if (norm2(xy_t - [0.5_dp, 0.5_dp]) < 0.25_dp) then
    !    sol = 1
    ! else
    !    sol = 0
    ! end if
  end function solution

  subroutine fluxes_upwind1(boxes, id, flux_args)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp), intent(in)        :: flux_args(:)
    real(dp)                    :: inv_dr
    integer                     :: nc

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    ! Diffusion
    boxes(id)%fx(:,:,i_phi) = boxes(id)%cc(0:nc, 1:nc, i_phi) - boxes(id)%cc(1:nc+1, 1:nc, i_phi)
    boxes(id)%fx(:,:,i_phi) = boxes(id)%fx(:,:,i_phi) * flux_args(1) * inv_dr
    boxes(id)%fy(:,:,i_phi) = boxes(id)%cc(1:nc, 0:nc, 1) - boxes(id)%cc(1:nc, 1:nc+1, 1)
    boxes(id)%fy(:,:,i_phi) = boxes(id)%fy(:,:,i_phi) * flux_args(1) * inv_dr

    ! Drift (1st order upwind, which is very diffusive!)
    if (flux_args(2) > 0) then
       boxes(id)%fx(:,:,i_phi) = boxes(id)%fx(:,:,i_phi) + &
            flux_args(2) * boxes(id)%cc(0:nc, 1:nc, 1)
    else
       boxes(id)%fx(:,:,i_phi) = boxes(id)%fx(:,:,i_phi) + &
            flux_args(2) * boxes(id)%cc(1:nc+1, 1:nc, 1)
    end if
    if (flux_args(3) > 0) then
       boxes(id)%fy(:,:,i_phi) = boxes(id)%fy(:,:,i_phi) + &
            flux_args(3) * boxes(id)%cc(1:nc, 0:nc, 1)
    else
       boxes(id)%fy(:,:,i_phi) = boxes(id)%fy(:,:,i_phi) + &
            flux_args(3) * boxes(id)%cc(1:nc, 1:nc+1, 1)
    end if
  end subroutine fluxes_upwind1

  subroutine fluxes_centdif(boxes, id, flux_args)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp), intent(in)        :: flux_args(:)
    real(dp)                    :: inv_dr
    integer                     :: nc

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    ! Diffusion
    boxes(id)%fx(:,:,i_phi) = (boxes(id)%cc(0:nc, 1:nc, i_phi) - &
         boxes(id)%cc(1:nc+1, 1:nc, i_phi)) * flux_args(1) * inv_dr
    boxes(id)%fy(:,:,i_phi) = (boxes(id)%cc(1:nc, 0:nc, i_phi) - &
         boxes(id)%cc(1:nc, 1:nc+1, i_phi))  * flux_args(1) * inv_dr

    boxes(id)%fx(:,:,i_phi) = boxes(id)%fx(:,:,i_phi) + flux_args(2) * 0.5_dp * &
         (boxes(id)%cc(0:nc, 1:nc, 1) + boxes(id)%cc(1:nc+1, 1:nc, 1))
    boxes(id)%fy(:,:,i_phi) = boxes(id)%fy(:,:,i_phi) + flux_args(3) * 0.5_dp * &
         (boxes(id)%cc(1:nc, 0:nc, 1) + boxes(id)%cc(1:nc, 1:nc+1, 1))
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

  subroutine fluxes_koren(boxes, id, flux_args)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp), intent(in)        :: flux_args(:)
    real(dp)                    :: inv_dr
    real(dp)                    :: tmp, gradp, gradc, gradn
    real(dp)                    :: gc_data(boxes(id)%n_cell, a2_num_neighbors)
    integer                     :: i, j, nc

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    call a2_gc2_box(boxes, id, i_phi, a2_gc2_prolong1, &
         have_no_bc2, gc_data, nc)

    ! x-fluxes interior, advective part with flux limiter
    do j = 1, nc
       do i = 1, nc+1
          gradc = boxes(id)%cc(i, j, i_phi) - boxes(id)%cc(i-1, j, i_phi)
          if (flux_args(2) < 0.0_dp) then
             if (i == nc+1) then
                tmp = gc_data(j, a2_nb_hx)
             else
                tmp = boxes(id)%cc(i+1, j, i_phi)
             end if

             gradn = tmp - boxes(id)%cc(i, j, i_phi)
             boxes(id)%fx(i, j, i_phi) = flux_args(2) * &
                  (boxes(id)%cc(i, j, i_phi) - koren_mlim(gradc, gradn))
          else                  ! flux_args(2) > 0
             if (i == 1) then
                tmp = gc_data(j, a2_nb_lx)
             else
                tmp = boxes(id)%cc(i-2, j, i_phi)
             end if

             gradp = boxes(id)%cc(i-1, j, i_phi) - tmp
             boxes(id)%fx(i, j, i_phi) = flux_args(2) * &
                  (boxes(id)%cc(i-1, j, i_phi) + koren_mlim(gradc, gradp))
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
          boxes(id)%fx(i, j, i_phi) = boxes(id)%fx(i, j, i_phi) - &
               flux_args(1) * gradc * inv_dr
       end do
    end do

    ! y-fluxes interior, advective part with flux limiter
    do j = 1, nc+1
       do i = 1, nc
          gradc = boxes(id)%cc(i, j, i_phi) - boxes(id)%cc(i, j-1, i_phi)
          if (flux_args(3) < 0.0_dp) then
             if (j == nc+1) then
                tmp = gc_data(i, a2_nb_hy)
             else
                tmp = boxes(id)%cc(i, j+1, i_phi)
             end if

             gradn = tmp - boxes(id)%cc(i, j, i_phi)
             boxes(id)%fy(i, j, i_phi) = flux_args(3) * &
                  (boxes(id)%cc(i, j, i_phi) - koren_mlim(gradc, gradn))
          else                  ! flux_args(3) > 0
             if (j == 1) then
                tmp = gc_data(i, a2_nb_ly)
             else
                tmp = boxes(id)%cc(i, j-2, i_phi)
             end if

             gradp = boxes(id)%cc(i, j-1, i_phi) - tmp
             boxes(id)%fy(i, j, i_phi) = flux_args(3) * &
                  (boxes(id)%cc(i, j-1, i_phi) + koren_mlim(gradc, gradp))
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
          boxes(id)%fy(i, j, i_phi) = boxes(id)%fy(i, j, i_phi) - &
               flux_args(1) * gradc * inv_dr
       end do
    end do
  end subroutine fluxes_koren

  subroutine update_solution(box, dt)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr
    integer                     :: nc

    nc                    = box%n_cell
    inv_dr                = 1/box%dr
    box%cc(1:nc, 1:nc, i_phi) = box%cc(1:nc, 1:nc, i_phi) + dt(1) * ( &
         (box%fx(1:nc, :, i_phi) - box%fx(2:nc+1, :, i_phi)) * inv_dr + &
         (box%fy(:, 1:nc, i_phi) - box%fy(:, 2:nc+1, i_phi)) * inv_dr)
  end subroutine update_solution

  subroutine average_phi(box)
    type(box2_t), intent(inout) :: box
    box%cc(:, :, i_phi) = 0.5_dp * (box%cc(:, :, i_phi) + box%cc(:, :, i_phi_old))
  end subroutine average_phi

  subroutine prolong_to_new_children(tree, ref_info)
    use m_a2_prolong
    type(a2_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id

    do lvl = 1, tree%max_lvl
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          ! Linear prolongation will not strictly conserve phi
          call a2_prolong1_to(tree%boxes, id, i_phi)
       end do

       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a2_gc_box(tree%boxes, id, i_phi, &
               a2_gc_interp_lim, have_no_bc)
       end do
    end do
  end subroutine prolong_to_new_children

  subroutine have_no_bc(boxes, id, i, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, i, iv
    stop "We have no boundary conditions in this example"
    boxes(id)%cc(1, i, iv) = 0    ! Prevent warning unused
  end subroutine have_no_bc

  subroutine have_no_bc2(boxes, id, nb, iv, gc_data, nc)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv, nc
    real(dp), intent(out)       :: gc_data(nc)
    stop "We have no boundary conditions in this example"
    boxes(id)%cc(1, nb, iv) = 0    ! Prevent warning unused
    gc_data = nb                  ! idem
  end subroutine have_no_bc2

end program test_drift_diff
