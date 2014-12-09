program test_drift_diff
  use m_afivo_2d

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: box_size    = 16

  ! Cell-centered variables
  integer, parameter :: n_var_cell = 8
  integer, parameter :: i_elec     = 1
  integer, parameter :: i_pion     = 2
  integer, parameter :: i_elec_old = 3
  integer, parameter :: i_pion_old = 4
  integer, parameter :: i_pot      = 5
  integer, parameter :: i_tmp      = 6
  integer, parameter :: i_rhs      = 7
  integer, parameter :: i_res      = 8
  character(len=10)  :: cc_names(n_var_cell) = &
       [character(len=10) :: "elec", "pion", "elec_old", &
       "pion_old", "pot", "fld", "rhs", "res"]

  ! Face-centered variables
  integer, parameter :: n_var_face = 3
  integer, parameter :: i_flux     = 1
  integer, parameter :: i_fld      = 2
  integer, parameter :: i_fld_old  = 3

  type(a2_t)         :: tree
  integer            :: i, n, n_steps, id
  integer            :: ix_list(2, 1)
  integer            :: nb_list(4, 1)
  integer            :: n_boxes_init = 10*1000
  integer            :: n_lvls_max  = 20
  integer            :: n_changes, output_cnt
  real(dp)           :: dr, dt, time, end_time
  real(dp)           :: dt_adapt, dt_output
  character(len=40)  :: fname
  logical            :: write_out

  real(dp), parameter :: diff_coeff = 0.1_dp
  real(dp), parameter :: mobility = 0.03_dp

  dr = 4.0_dp / box_size

  ! Initialize tree
  call a2_init(tree, n_lvls_max, n_boxes_init, box_size, n_var_cell=n_var_cell, &
       n_var_face=n_var_face, dr = dr, r_min = [0.0_dp, 0.0_dp])

  ! Set up geometry
  id             = 1
  ix_list(:, id) = [1,1] ! Set index of box
  nb_list(:, id) = -1    ! We use -1 as physical boundary

  ! Create the base mesh
  call a2_set_base(tree, ix_list, nb_list)

  ! Set up the initial conditions
  do i = 1, 20
     ! We should only set the finest level, but this also works
     call a2_loop_box(tree, set_init_cond)
     call compute_fld(tree)
     call a2_adjust_refinement(tree, ref_func, n_changes)
     if (n_changes == 0) exit
  end do

  ! Restrict the initial conditions
  call a2_restrict_tree(tree, i_elec)
  call a2_restrict_tree(tree, i_pion)

  i          = 0
  output_cnt = 0
  time       = 0
  dt_adapt   = 0.01_dp
  dt_output  = 0.05_dp
  end_time   = 5.0_dp

  !$omp parallel private(n)
  do
     !$omp single
     i       = i + 1
     dt      = 0.5_dp * get_max_dt(tree)

     n_steps = ceiling(dt_adapt/dt)
     dt      = dt_adapt / n_steps
     time    = time + dt_adapt

     if (output_cnt * dt_output < time) then
        write_out = .true.
        output_cnt = output_cnt + 1
        write(fname, "(A,I0,A)") "test_drift_diff_", output_cnt, ".vtu"
     else
        write_out = .false.
     end if
     !$omp end single

     if (write_out) call a2_write_tree(tree, trim(fname), &
          cc_names, output_cnt, time)

     if (time > end_time) exit

     do n = 1, n_steps
        ! Copy previous solution
        call a2_tree_copy_cc(tree, i_elec, i_elec_old)
        call a2_tree_copy_cc(tree, i_pion, i_pion_old)

        ! Take a half time step
        call a2_loop_boxes(tree, fluxes_koren)
        call a2_loop_box_arg(tree, update_solution, [0.5_dp * dt])
        call a2_restrict_tree(tree, i_elec)
        call a2_restrict_tree(tree, i_pion)
        call a2_gc_sides(tree, i_elec, a2_sides_extrap, sides_bc)
        call a2_gc_sides(tree, i_pion, a2_sides_extrap, sides_bc)

        ! Calculate fluxes
        call a2_loop_boxes(tree, fluxes_koren)

        ! Copy back old elec/pion, and take full time step
        call a2_tree_copy_cc(tree, i_elec_old, i_elec)
        call a2_tree_copy_cc(tree, i_pion_old, i_pion)
        call a2_loop_box_arg(tree, update_solution, [dt])
        call a2_restrict_tree(tree, i_elec)
        call a2_restrict_tree(tree, i_pion)
        call a2_gc_sides(tree, i_elec, a2_sides_extrap, sides_bc)
        call a2_gc_sides(tree, i_pion, a2_sides_extrap, sides_bc)
     end do

     call a2_restrict_tree(tree, i_elec)
     call a2_restrict_tree(tree, i_pion)
     call a2_gc_sides(tree, i_elec, a2_sides_extrap, sides_bc)
     call a2_gc_sides(tree, i_pion, a2_sides_extrap, sides_bc)
     call a2_gc_corners(tree, i_elec, a2_corners_extrap, sides_bc)
     call a2_gc_corners(tree, i_pion, a2_corners_extrap, sides_bc)

     call a2_adjust_refinement(tree, ref_func)
     call a2_loop_boxes(tree, prolong_to_new_children)
  end do
  !$omp end parallel

  call a2_destroy(tree)

contains

  integer function ref_func(box)
    type(box2_t), intent(in) :: box
    if (box%lvl < 4) then
       ref_func = a5_do_ref
    else
       ref_func = a5_rm_ref
    end if
  end function ref_func

  real(dp) function get_max_dt(tree)
    type(a2_t), intent(in) :: tree
    get_max_dt = 1.0e-12_dp
  end function get_max_dt

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 0, nc+1
       do i = 0, nc+1
          xy = a2_r_cc(box, [i,j])
          if (norm2(xy - 2) < 1) then
             box%cc(i, j, i_elec) = 1e15_dp
             box%cc(i, j, i_pion) = 1e15_dp
          else
             box%cc(i, j, i_elec) = 0
             box%cc(i, j, i_pion) = 0
          end if
       end do
    end do
  end subroutine set_init_cond

  subroutine compute_fld(tree)
    type(a2_t), intent(in) :: tree

    real(dp), parameter :: UC_eps0 = 8.8541878176d-12 ! permitivity of vacuum (SI)
    real(dp), parameter :: UC_elem_charge = 1.6022d-19 ! the elementary charge in Coulombs
    real(dp), parameter :: fac = UC_eps0 / UC_elem_charge

    ! Set rhs


    ! Solve Poisson's eq. using multigrid

    !
  end subroutine compute_fld

  elemental function limiter_koren(theta)
    real(dp), intent(in) :: theta
    real(dp)             :: limiter_koren
    real(dp), parameter  :: one_sixth = 1.0_dp / 6.0_dp
    limiter_koren = max(0.0d0, min(1.0_dp, theta, (1.0_dp + 2.0_dp * theta) * one_sixth))
  end function limiter_koren

  elemental function ratio(numerator, denominator)
    real(dp), intent(in) :: numerator, denominator
    real(dp)             :: ratio
    real(dp), parameter  :: eps = epsilon(1.0d0)
    ! Avoid division by zero, and ensure that at zero gradients we have a ratio of 1
    ratio = (sign(eps, numerator) + numerator) / &
         (denominator + sign(eps, denominator))
  end function ratio

  subroutine fluxes_koren(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp)                    :: inv_dr, theta
    real(dp)                    :: gradp, gradc, gradn
    integer                     :: i, j, nc, nb_id

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    ! x-fluxes interior, advective part with flux limiter
    do j = 1, nc
       do i = 1, nc+1
          gradc = boxes(id)%cc(i, j, i_elec) - boxes(id)%cc(i-1, j, i_elec)
          if (boxes(id)%fx(i, j, i_fld) * mobility < 0.0_dp) then

             if (i == nc+1) then
                nb_id = boxes(id)%neighbors(a2_nb_hx)
                if (nb_id > a5_no_box) then
                   gradn = boxes(nb_id)%cc(2, j, i_elec) - boxes(id)%cc(i, j, i_elec)
                else
                   gradn = 0
                end if
             else
                gradn = boxes(id)%cc(i+1, j, i_elec) - boxes(id)%cc(i, j, i_elec)
             end if

             theta = ratio(gradc, gradn)
             boxes(id)%fx(i, j, i_elec) = boxes(id)%fx(i, j, i_fld) * mobility * &
                  (boxes(id)%cc(i, j, i_elec) - limiter_koren(theta) * gradn)
          else                  ! boxes(id)%fx(i, j, i_fld) * mobility > 0

             if (i == 1) then
                nb_id = boxes(id)%neighbors(a2_nb_lx)
                if (nb_id > a5_no_box) then
                   gradp = boxes(id)%cc(i-1, j, i_elec) - boxes(nb_id)%cc(nc-1, j, i_elec)
                else
                   gradp = 0
                end if
             else
                gradp = boxes(id)%cc(i-1, j, i_elec) - boxes(id)%cc(i-2, j, i_elec)
             end if

             theta = ratio(gradc, gradp)
             boxes(id)%fx(i, j, i_elec) = boxes(id)%fx(i, j, i_fld) * mobility * &
                  (boxes(id)%cc(i-1, j, i_elec) + limiter_koren(theta) * gradp)
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
          boxes(id)%fx(i, j, i_elec) = boxes(id)%fx(i, j, i_elec) - &
               diff_coeff * gradc * inv_dr
       end do
    end do

    ! y-fluxes interior, advective part with flux limiter
    do j = 1, nc+1
       do i = 1, nc
          gradc = boxes(id)%cc(i, j, i_elec) - boxes(id)%cc(i, j-1, i_elec)
          if (boxes(id)%fy(i, j, i_fld) * mobility < 0.0_dp) then

             if (j == nc+1) then
                nb_id = boxes(id)%neighbors(a2_nb_hy)
                if (nb_id > a5_no_box) then
                   gradp = boxes(nb_id)%cc(i, 2, i_elec) - boxes(id)%cc(i, j, i_elec)
                else
                   gradp = 0
                end if
             else
                gradn = boxes(id)%cc(i, j+1, i_elec) - boxes(id)%cc(i, j, i_elec)
             end if

             theta = ratio(gradc, gradn)
             boxes(id)%fy(i, j, i_elec) = boxes(id)%fy(i, j, i_fld) * mobility * &
                  (boxes(id)%cc(i, j, i_elec) - limiter_koren(theta) * gradn)
          else                  ! boxes(id)%fy(i, j, i_fld) * mobility > 0

             if (j == 1) then
                nb_id = boxes(id)%neighbors(a2_nb_ly)
                if (nb_id > a5_no_box) then
                   gradn = boxes(id)%cc(i, j-1, i_elec) - boxes(nb_id)%cc(i, nc-1, i_elec)
                else
                   gradn = 0
                end if
             else
                gradn = boxes(id)%cc(i, j-1, i_elec) - boxes(id)%cc(i, j-2, i_elec)
             end if

             theta = ratio(gradc, gradn)
             boxes(id)%fy(i, j, i_elec) = boxes(id)%fy(i, j, i_fld) * mobility * &
                  (boxes(id)%cc(i, j-1, i_elec) + limiter_koren(theta) * gradn)
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
          boxes(id)%fy(i, j, i_elec) = boxes(id)%fy(i, j, i_elec) - &
               diff_coeff * gradc * inv_dr
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
    box%cc(1:nc, 1:nc, 1) = box%cc(1:nc, 1:nc, 1) + dt(1) * ( &
         (box%fx(1:nc, :, 1) - box%fx(2:nc+1, :, 1)) * inv_dr + &
         (box%fy(:, 1:nc, 1) - box%fy(:, 2:nc+1, 1)) * inv_dr)
  end subroutine update_solution

  subroutine prolong_to_new_children(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in) :: id

    if (btest(boxes(id)%tag, a5_bit_new_children)) then
       call a2_prolong1_from(boxes, id, i_elec, .true.)
       call a2_prolong1_from(boxes, id, i_pion, .true.)
       boxes(id)%tag = ibclr(boxes(id)%tag, a5_bit_new_children)
    end if
  end subroutine prolong_to_new_children

  subroutine sides_bc(boxes, id, nb, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    integer                     :: nc

    nc = boxes(id)%n_cell

    if (boxes(id)%neighbors(nb) == -1) then
       select case (nb)
       case (a2_nb_lx)
          ! Neumann zero
          boxes(id)%cc(0, 1:nc, iv) = boxes(id)%cc(1, 1:nc, iv)
       case (a2_nb_hx)
          ! Neumann zero
          boxes(id)%cc(nc+1, 1:nc, iv) = boxes(id)%cc(nc, 1:nc, iv)
       case (a2_nb_ly)
          ! Dirichlet zero
          boxes(id)%cc(1:nc, 0, iv) = -boxes(id)%cc(1:nc, 1, iv)
       case (a2_nb_hy)
          ! Dirichlet one
          boxes(id)%cc(1:nc, nc+1, iv) = 2 - boxes(id)%cc(1:nc, nc, iv)
       end select
    end if
  end subroutine sides_bc

end program test_drift_diff
