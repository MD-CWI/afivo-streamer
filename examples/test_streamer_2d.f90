program test_drift_diff
  use m_afivo_2d
  use m_mg_2d
  use m_mg_diel

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: box_size    = 8

  ! Cell-centered variables
  integer, parameter :: n_var_cell = 9
  integer, parameter :: i_elec     = 1
  integer, parameter :: i_pion     = 2
  integer, parameter :: i_elec_old = 3
  integer, parameter :: i_pion_old = 4
  integer, parameter :: i_phi      = 5
  integer, parameter :: i_tmp      = 6
  integer, parameter :: i_rhs      = 7
  integer, parameter :: i_res      = 8
  integer, parameter :: i_eps      = 9
  character(len=10)  :: cc_names(n_var_cell) = &
       [character(len=10) :: "elec", "pion", "elec_old", &
       "pion_old", "phi", "fld", "rhs", "res", "eps"]

  ! Face-centered variables
  integer, parameter :: n_var_face = 2
  integer, parameter :: i_flux     = 1
  integer, parameter :: i_fld      = 2

  integer, parameter :: n_fmg_cycles = 1

  type(a2_t)         :: tree
  type(mg2_t)        :: mg
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
  real(dp), parameter :: mobility = -0.03_dp
  real(dp), parameter :: domain_len = 32.0e-3_dp

  dr = domain_len / box_size

  ! Initialize tree
  call a2_init(tree, n_lvls_max, n_boxes_init, box_size, n_var_cell=n_var_cell, &
       n_var_face=n_var_face, dr = dr, r_min = [0.0_dp, 0.0_dp], coarsen_to=4)

  ! Set the multigrid options
  mg%i_phi        = i_phi
  mg%i_tmp        = i_tmp
  mg%i_rhs        = i_rhs
  mg%i_res        = i_res
  mg%i_eps        = i_eps
  mg%n_cycle_base = 8

  mg%sides_bc     => sides_bc_pot
  mg%box_op       => lpl_box_diel
  mg%box_gsrb     => gsrb_lpl_box_diel
  mg%box_corr     => corr_lpl_box_diel

  call mg2_init_mg(mg)

  ! Set up geometry
  id             = 1
  ix_list(:, id) = [1,1] ! Set index of box
  nb_list(:, id) = -1    ! We use -1 as physical boundary

  ! Create the base mesh
  call a2_set_base(tree, ix_list, nb_list)

  ! Set up the initial conditions
  do i = 1, 10
     call a2_loop_box(tree, set_init_cond)
     call a2_restrict_tree(tree, i_rhs)
     call a2_restrict_tree(tree, i_phi)
     call compute_fld(tree, n_fmg_cycles)
     call a2_adjust_refinement(tree, ref_func, n_changes)
     if (n_changes == 0) exit
  end do

  call a2_loop_box(tree, set_init_cond)

  i                = 0
  output_cnt       = 0
  time             = 0
  dt_adapt         = 1.0e-11_dp
  dt_output        = 2.5e-10_dp
  end_time         = 25.0e-9_dp

  !$omp parallel private(n)
  do
     !$omp single
     i             = i + 1
     dt            = 0.5_dp * get_max_dt(tree)

     n_steps = ceiling(dt_adapt/dt)
     dt      = dt_adapt / n_steps
     time    = time + dt_adapt

     if (output_cnt * dt_output < time) then
        write_out = .true.
        output_cnt = output_cnt + 1
        write(fname, "(A,I0,A)") "test_str2d_", output_cnt, ".vtu"
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

        ! Two forward Euler steps over dt
        do i = 1, 2
           call a2_loop_boxes(tree, fluxes_koren)
           call a2_consistent_fluxes(tree, [i_flux])
           call compute_fld(tree, n_fmg_cycles)
           call a2_loop_box_arg(tree, update_solution, [dt])
           call a2_restrict_tree(tree, i_elec)
           call a2_restrict_tree(tree, i_pion)
           call a2_gc_sides(tree, i_elec, a2_sides_interp, sides_bc_dens)
           call a2_gc_sides(tree, i_pion, a2_sides_interp, sides_bc_dens)
        end do

        ! Take average of phi_old and phi
        call a2_loop_box(tree, average_dens)
     end do

     call a2_adjust_refinement(tree, ref_func)
     call a2_loop_boxes(tree, prolong_to_new_children)
     call compute_fld(tree, n_fmg_cycles)
  end do
  !$omp end parallel

  call a2_destroy(tree)

contains

  integer function ref_func(box)
    type(box2_t), intent(in) :: box
    integer :: nc
    real(dp) :: max_fld, max_dns, dr, alpha
    nc = box%n_cell
    max_fld = maxval(abs(box%cc(1:nc, 1:nc, i_tmp)))
    max_dns = maxval(box%cc(1:nc, 1:nc, i_elec))
    dr = box%dr
    alpha = get_alpha(max_fld)
    ! print *, alpha, 1/alpha, dr

    if (box%lvl < 4) then
       ref_func = a5_do_ref
    else if (box%lvl > 10) then
       ref_func = a5_rm_ref
    else if (max_dns > 1.0e15_dp .and. alpha * dr > 1) then
       ref_func = a5_do_ref
    else if (max_dns > 1.0e15_dp .and. dr > 1e-4_dp .and. max_fld > 1e6_dp) then
       ref_func = a5_do_ref
    else
       ref_func = a5_rm_ref
    end if
  end function ref_func

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2), xy0(2), xy1(2), sigma
    real(dp)                    :: bg_dens, dns

    nc = box%n_cell
    xy0 = [0.3_dp, 0.5_dp] * domain_len
    xy1 = xy0
    xy1(2) = xy1(2) + 2e-3_dp   ! 2 mm
    sigma = 2e-4_dp
    bg_dens = 5e13

    do j = 0, nc+1
       do i = 0, nc+1
          xy = a2_r_cc(box, [i,j])
          dns = 5e19_dp * rod_dens(xy, xy0, xy1, sigma, 3)

          if (xy(2) > 0.75_dp * domain_len) then
             box%cc(i, j, i_eps) = 5.0_dp
             box%cc(i, j, i_elec) = 0
             box%cc(i, j, i_pion) = 0
          else
             box%cc(i, j, i_eps) = 1.0_dp
             box%cc(i, j, i_elec) = bg_dens + dns
             box%cc(i, j, i_pion) = bg_dens + dns
          end if

       end do
    end do

    box%cc(:, :, i_phi) = 0
  end subroutine set_init_cond

  real(dp) function get_max_dt(tree)
    type(a2_t), intent(in) :: tree
    real(dp), parameter    :: UC_eps0        = 8.8541878176d-12
    real(dp), parameter    :: UC_elem_charge = 1.6022d-19
    real(dp)               :: max_fld, max_dns, dr_min(2)
    real(dp)               :: dt_cfl, dt_dif, dt_drt
    integer                :: lvl, i, id

    max_fld = -huge(1.0_dp)
    max_dns = -huge(1.0_dp)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          max_fld = max(max_fld, maxval(tree%boxes(id)%cc(:,:, i_tmp)))
          max_dns = max(max_dns, maxval(tree%boxes(id)%cc(:,:, i_elec)))
       end do
    end do

    dr_min = a2_min_dr(tree)
    ! CFL condition
    dt_cfl = dr_min(1) / (abs(mobility) * max_fld)
    ! Diffusion condition
    dt_dif = dr_min(1)**2 / (diff_coeff)
    ! Dielectric relaxation time
    dt_drt = UC_eps0 / (UC_elem_charge * abs(mobility) * max_dns)
    get_max_dt = 0.5_dp * min(dt_cfl, dt_dif, dt_drt)
  end function get_max_dt

  real(dp) function rod_dens(xy, xy0, xy1, sigma, falloff_type)
    real(dp), intent(in) :: xy(2), xy0(2), xy1(2), sigma
    integer, intent(in) :: falloff_type
    real(dp) :: line_len2, distance, temp
    real(dp) :: projection(2)

    line_len2 = sum((xy1 - xy0)**2)
    temp = sum((xy - xy0) * (xy1 - xy0)) / line_len2

    if (temp < 0.0_dp) then
       distance = sqrt(sum((xy-xy0)**2))
    else if (temp > 1.0_dp) then
       distance = sqrt(sum((xy-xy1)**2))
    else
       projection = xy0 + temp * (xy1 - xy0)
       distance = sqrt(sum((xy-projection)**2))
    end if

    select case (falloff_type)
    case (1)                    ! Sigmoid
       rod_dens    = 2 / (1 + exp(distance / sigma))
    case (2)                    ! Gaussian
       rod_dens    = exp(-(distance/sigma)**2)
    case (3)                    ! Smooth-step
       if (distance < sigma) then
          rod_dens = 1
       else if (distance < 2 * sigma) then
          temp = distance/sigma - 1
          rod_dens = (1- (3 * temp**2 - 2 * temp**3))
       else
          rod_dens = 0.0_dp
       end if
    case default
       rod_dens = 0.0_dp
    end select
  end function rod_dens

  real(dp) function gaussian_2d(x, x0, sigma)
    real(dp), intent(in) :: x(2), x0(2), sigma
    gaussian_2d = exp(-0.5_dp/sigma**2 * sum((x-x0)**2))
  end function gaussian_2d

  subroutine compute_fld(tree, n_fmg)
    type(a2_t), intent(inout) :: tree
    integer, intent(in) :: n_fmg

    real(dp), parameter :: UC_eps0 = 8.8541878176d-12
    real(dp), parameter :: UC_elem_charge = 1.6022d-19
    real(dp), parameter :: fac = UC_elem_charge / UC_eps0
    integer :: lvl, i, id, nc

    nc = tree%n_cell
    ! Set rhs
    do lvl = 1, tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          tree%boxes(id)%cc(1:nc, 1:nc, i_rhs) = fac * (&
               tree%boxes(id)%cc(1:nc, 1:nc, i_elec) - &
               tree%boxes(id)%cc(1:nc, 1:nc, i_pion))
       end do
       !$omp end do nowait
    end do

    ! Restrict the rhs
    call a2_restrict_tree(tree, i_rhs)

    ! Perform n_fmg full-multigrid cycles
    do i = 1, n_fmg
       call mg2_fas_fmg(tree, mg)
    end do

    ! Compute field from potential
    call a2_loop_box(tree, fld_from_pot)
  end subroutine compute_fld

  subroutine fld_from_pot(box)
    type(box2_t), intent(inout) :: box
    integer                     :: nc
    real(dp)                    :: inv_dr

    nc = box%n_cell
    inv_dr = 1 / box%dr
    box%fx(:, :, i_fld) = inv_dr * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi))
    box%fy(:, :, i_fld) = inv_dr * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi))

    ! Compute fields at the boundaries of the box, where eps can change
    box%fx(1, :, i_fld) = 2 * inv_dr * &
         (box%cc(0, 1:nc, i_phi) - box%cc(1, 1:nc, i_phi)) * &
         box%cc(0, 1:nc, i_eps) / &
         (box%cc(1, 1:nc, i_eps) + box%cc(0, 1:nc, i_eps))
    box%fx(nc+1, :, i_fld) = 2 * inv_dr * &
         (box%cc(nc, 1:nc, i_phi) - box%cc(nc+1, 1:nc, i_phi)) * &
         box%cc(nc+1, 1:nc, i_eps) / &
         (box%cc(nc+1, 1:nc, i_eps) + box%cc(nc, 1:nc, i_eps))
    box%fy(:, 1, i_fld) = 2 * inv_dr * &
         (box%cc(1:nc, 0, i_phi) - box%cc(1:nc, 1, i_phi)) * &
         box%cc(1:nc, 0, i_eps) / &
         (box%cc(1:nc, 1, i_eps) + box%cc(1:nc, 0, i_eps))
    box%fy(:, nc+1, i_fld) = 2 * inv_dr * &
         (box%cc(1:nc, nc, i_phi) - box%cc(1:nc, nc+1, i_phi)) * &
         box%cc(1:nc, nc+1, i_eps) / &
         (box%cc(1:nc, nc+1, i_eps) + box%cc(1:nc, nc, i_eps))

    box%cc(1:nc, 1:nc, i_tmp) = sqrt(&
         0.25_dp * (box%fx(1:nc, 1:nc, i_fld) + box%fx(2:nc+1, 1:nc, i_fld))**2 + &
         0.25_dp * (box%fy(1:nc, 1:nc, i_fld) + box%fy(1:nc, 2:nc+1, i_fld))**2)
  end subroutine fld_from_pot

  elemental function limiter_koren(theta)
    real(dp), intent(in) :: theta
    real(dp)             :: limiter_koren
    real(dp), parameter  :: one_sixth = 1.0_dp / 6.0_dp
    limiter_koren = max(0.0d0, min(1.0_dp, theta, &
         (1.0_dp + 2.0_dp * theta) * one_sixth))
  end function limiter_koren

  elemental function ratio(numerator, denominator)
    real(dp), intent(in) :: numerator, denominator
    real(dp)             :: ratio
    if (denominator /= 0.0_dp) then
       ratio = numerator / denominator
    else
       ratio = numerator * huge(1.0_dp)
    end if
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
                   gradn = boxes(nb_id)%cc(i, 2, i_elec) - boxes(id)%cc(i, j, i_elec)
                else
                   gradn = 0
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
                   gradp = boxes(id)%cc(i, j-1, i_elec) - boxes(nb_id)%cc(i, nc-1, i_elec)
                else
                   gradp = 0
                end if
             else
                gradp = boxes(id)%cc(i, j-1, i_elec) - boxes(id)%cc(i, j-2, i_elec)
             end if

             theta = ratio(gradc, gradp)
             boxes(id)%fy(i, j, i_elec) = boxes(id)%fy(i, j, i_fld) * mobility * &
                  (boxes(id)%cc(i, j-1, i_elec) + limiter_koren(theta) * gradp)
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
          boxes(id)%fy(i, j, i_elec) = boxes(id)%fy(i, j, i_elec) - &
               diff_coeff * gradc * inv_dr
       end do
    end do
  end subroutine fluxes_koren

  subroutine average_dens(box)
    type(box2_t), intent(inout) :: box
    box%cc(:, :, i_elec) = 0.5_dp * (box%cc(:, :, i_elec) + box%cc(:, :, i_elec_old))
    box%cc(:, :, i_pion) = 0.5_dp * (box%cc(:, :, i_pion) + box%cc(:, :, i_pion_old))
  end subroutine average_dens

  subroutine update_solution(box, dt)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr, src, fld
    integer                     :: i, j, nc

    nc                    = box%n_cell
    inv_dr                = 1/box%dr
    do j = 1, nc
       do i = 1, nc
          fld = box%cc(i,j, i_tmp)
          src = abs(mobility * fld) * get_alpha(fld) * &
               dt(1) * box%cc(i, j, i_elec)
          box%cc(i, j, i_elec) = box%cc(i, j, i_elec) + src
          box%cc(i, j, i_pion) = box%cc(i, j, i_pion) + src
       end do
    end do

    box%cc(1:nc, 1:nc, i_elec) = box%cc(1:nc, 1:nc, i_elec) + dt(1) * ( &
         (box%fx(1:nc, :, i_elec) - box%fx(2:nc+1, :, i_elec)) * inv_dr + &
         (box%fy(:, 1:nc, i_elec) - box%fy(:, 2:nc+1, i_elec)) * inv_dr)

  end subroutine update_solution

  real(dp) function get_alpha(fld)
    real(dp), intent(in) :: fld
    ! Breakdown fld of 3 MV/m
    get_alpha = max(0.0_dp, 1e5_dp * exp(1 - 1e7_dp/(abs(fld)+epsilon(1.0_dp))) - 9697.2_dp)
  end function get_alpha

  subroutine prolong_to_new_children(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in) :: id

    if (btest(boxes(id)%tag, a5_bit_new_children)) then
       call a2_prolong1_from(boxes, id, i_elec, .true.)
       call a2_prolong1_from(boxes, id, i_pion, .true.)
       call a2_prolong1_from(boxes, id, i_phi, .true.)
       call a2_prolong0_from(boxes, id, i_eps, .true.)
       boxes(id)%tag = ibclr(boxes(id)%tag, a5_bit_new_children)
    end if
  end subroutine prolong_to_new_children

  subroutine sides_bc_pot(boxes, id, nb, iv)
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
          ! Dirichlet
          boxes(id)%cc(1:nc, nc+1, iv) = 2 * 2.5e6_dp * domain_len &
               - boxes(id)%cc(1:nc, nc, iv)
       end select
    end if
  end subroutine sides_bc_pot

  subroutine sides_bc_dens(boxes, id, nb, iv)
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
          ! Neumann zero
          boxes(id)%cc(1:nc, 0, iv) = boxes(id)%cc(1:nc, 1, iv)
       case (a2_nb_hy)
          ! Neumann zero
          boxes(id)%cc(1:nc, nc+1, iv) = boxes(id)%cc(1:nc, nc, iv)
       end select
    end if
  end subroutine sides_bc_dens

end program
