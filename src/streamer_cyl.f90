program streamer_cyl

  use m_a2_t
  use m_a2_core
  use m_a2_gc
  use m_a2_utils
  use m_a2_restrict
  use m_a2_mg
  use m_a2_io
  use m_write_silo
  use m_streamer

  implicit none

  character(len=name_len) :: sim_name, output_dir
  character(len=name_len) :: cfg_name, tmp_name, prev_name

  type(a2_t)        :: tree               ! This contains the full grid information
  type(mg2_t)       :: mg                 ! Multigrid option struct
  type(ref_info_t)  :: ref_info

  real(dp)          :: fld_mod_t0, fld_sin_amplitude
  real(dp)          :: fld_sin_freq, fld_lin_deriv
  real(dp)          :: GLOBAL_fld_val

  character(len=200) :: fname_axis, fname_stats

  call create_cfg(sim_cfg)

  sim_name = ""
  prev_name = ""
  do n = 1, command_argument_count()
     call get_command_argument(n, cfg_name)
     call CFG_read_file(sim_cfg, trim(cfg_name))

     call CFG_get(sim_cfg, "sim_name", tmp_name)
     if (sim_name == "") then
        sim_name = tmp_name
     else if (tmp_name /= "" .and. tmp_name /= prev_name) then
        sim_name = trim(sim_name) // "_" // trim(tmp_name)
     end if
     prev_name = tmp_name
  end do

  call CFG_get(sim_cfg, "end_time", end_time)
  call CFG_get(sim_cfg, "box_size", box_size)
  call CFG_get(sim_cfg, "output_dir", output_dir)
  call CFG_get(sim_cfg, "domain_len", domain_len)
  call CFG_get(sim_cfg, "applied_fld", applied_fld)
  call CFG_get(sim_cfg, "fld_mod_t0", fld_mod_t0)
  call CFG_get(sim_cfg, "fld_sin_amplitude", fld_sin_amplitude)
  call CFG_get(sim_cfg, "fld_sin_freq", fld_sin_freq)
  call CFG_get(sim_cfg, "fld_lin_deriv", fld_lin_deriv)
  call CFG_get(sim_cfg, "dt_output", dt_output)
  call CFG_get(sim_cfg, "num_steps_amr", n_steps_amr)
  call CFG_get(sim_cfg, "dt_max", dt_max)
  call CFG_get(sim_cfg, "epsilon_diel", epsilon_diel)

  if (trim(output_dir) == "") stop "No output directory given"
  tmp_name = trim(output_dir) // "/" // trim(sim_name) // "_config.txt"
  call CFG_write(sim_cfg, trim(tmp_name))
  print *, "Settings written to ", trim(tmp_name)

  ! Initialize the transport coefficients
  call init_transport_coeff(sim_cfg)

  ! Set the initial conditions from the configuration
  applied_voltage = -domain_len * applied_fld
  call get_init_cond(sim_cfg, init_cond, 2)

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the multigrid options. First define the variables to use
  mg%i_phi        = i_phi
  mg%i_tmp        = i_fld
  mg%i_rhs        = i_rhs
  mg%i_eps        = i_eps

  ! The number of cycles at the lowest level
  mg%n_cycle_base = 4

  ! Routines to use for ...
  mg%sides_bc    => sides_bc_pot ! Filling ghost cell on physical boundaries
  mg%box_op      => mg2_auto_op
  mg%box_corr    => mg2_auto_corr
  mg%box_gsrb    => mg2_auto_gsrb

  ! This routine always needs to be called when using multigrid
  call mg2_init_mg(mg)

  output_cnt = 0          ! Number of output files written
  time       = 0          ! Simulation time (all times are in s)

  ! Set up the initial conditions
  do
     call a2_loop_box(tree, set_init_cond)
     call compute_fld(tree, n_fmg_cycles, .true.)
     call a2_adjust_refinement(tree, set_ref_flags, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  if (photoi_enabled) &
       call set_photoionization(tree, photoi_eta, photoi_num_photons)

  do
     ! Get a new time step, which is at most dt_amr
     dt = get_max_dt(tree)

     if (dt < 1e-14) then
        print *, "dt getting too small, instability?"
        time = end_time + 1.0_dp
     end if

     ! Every dt_output, write output
     if (output_cnt * dt_output <= time) then
        write_out = .true.
        output_cnt = output_cnt + 1
        write(fname, "(A,I6.6)") trim(sim_name) // "_", output_cnt
        fname_axis = trim(output_dir) // "/" // trim(fname) // "_axis.txt"
        fname_stats = trim(output_dir) // "/" // trim(sim_name) // ".txt"
     else
        write_out = .false.
     end if

     if (time > end_time) exit

     ! We perform n_steps between mesh-refinements
     do n = 1, n_steps_amr

        if (photoi_enabled) &
             call set_photoionization(tree, photoi_eta, photoi_num_photons, dt)

        ! Copy previous solution
        call a2_tree_copy_cc(tree, i_elec, i_elec_old)
        call a2_tree_copy_cc(tree, i_pion, i_pion_old)

        ! Two forward Euler steps over dt
        do i = 1, 2
           time = time + dt

           ! First calculate fluxes
           call a2_loop_boxes_arg(tree, fluxes_koren, [dt], .true.)
           call a2_consistent_fluxes(tree, [f_elec])

           ! Update the solution
           call a2_loop_box_arg(tree, update_solution, [dt], .true.)

           ! Restrict the electron and ion densities to lower levels
           call a2_restrict_tree(tree, i_elec)
           call a2_restrict_tree(tree, i_pion)

           ! Fill ghost cells
           call a2_gc_tree(tree, i_elec, a2_gc_interp_lim, a2_gc_neumann)
           call a2_gc_tree(tree, i_pion, a2_gc_interp_lim, a2_gc_neumann)

           ! Compute new field on first iteration
           if (i == 1) call compute_fld(tree, n_fmg_cycles, .false.)
        end do

        time = time - dt        ! Go back one time step

        ! Take average of phi_old and phi (explicit trapezoidal rule)
        call a2_loop_box(tree, average_dens)

        ! Compute field with new density
        call compute_fld(tree, n_fmg_cycles, .false.)
     end do

     if (write_out) then
        call a2_write_silo(tree, fname, &
             cc_names, output_cnt, time, dir=output_dir, &
             fc_names=["fld_r", "fld_z"], ixs_fc=[f_fld])
        call write_streamer_properties(tree, fname_stats, &
             fname_axis, output_cnt, time)
     end if

     call a2_adjust_refinement(tree, set_ref_flags, ref_info)

     if (ref_info%n_add > 0 .or. ref_info%n_rm > 0) then
        ! For boxes which just have been refined, set data on their children
        call prolong_to_new_boxes(tree, ref_info)

        ! Compute the field on the new mesh
        call compute_fld(tree, n_fmg_cycles, .false.)

        ! This will every now-and-then clean up the data in the tree
        call a2_tidy_up(tree, 0.9_dp, 0.25_dp, 5000, .false.)
     end if

  end do

  call a2_destroy(tree)

contains

  ! Initialize the AMR tree
  subroutine init_tree(tree)
    type(a2_t), intent(inout) :: tree

    ! Variables used below to initialize tree
    real(dp)                  :: dr
    integer                   :: id
    integer                   :: ix_list(2, 1) ! Spatial indices of initial boxes
    integer                   :: nb_list(4, 1) ! Neighbors of initial boxes
    integer                   :: n_boxes_init = 1000

    dr = domain_len / box_size

    ! Initialize tree
    call a2_init(tree, box_size, n_var_cell, n_var_face, dr, &
         coarsen_to=2, n_boxes=n_boxes_init, coord=a5_cyl)

    ! Set up geometry
    id             = 1          ! One box ...
    ix_list(:, id) = [1,1]      ! With index 1,1 ...
    nb_list(:, id) = -1         ! And neighbors -1 (physical boundary)

    ! Create the base mesh
    call a2_set_base(tree, ix_list, nb_list)

  end subroutine init_tree

  ! Refinement function
  subroutine set_ref_flags(boxes, id, ref_flags)
    use m_geom
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flags(:)
    integer                  :: nc, n
    real(dp)                 :: crv_phi, dr2, max_fld, max_dns
    real(dp)                 :: boxlen, dist, alpha, adx

    nc        = boxes(id)%n_cell
    dr2       = boxes(id)%dr**2
    crv_phi   = dr2 * maxval(abs(boxes(id)%cc(1:nc, 1:nc, i_rhs)))
    max_fld   = maxval(boxes(id)%cc(1:nc, 1:nc, i_fld))
    max_dns   = maxval(boxes(id)%cc(1:nc, 1:nc, i_elec))
    alpha     = LT_get_col(td_tbl, i_alpha, max_fld)
    adx       = boxes(id)%dr * alpha

    if (adx < 0.1_dp .and. boxes(id)%dr < 2.0e-5_dp) &
         ref_flags(id) = a5_rm_ref
    if (adx < 0.1_dp .and. crv_phi < 4.0_dp .and. boxes(id)%dr < 5.0e-5_dp) &
         ref_flags(id) = a5_rm_ref

    if (time < 5.0e-9_dp) then
       boxlen = boxes(id)%n_cell * boxes(id)%dr

       do n = 1, init_cond%n_cond
          dist = GM_dist_line(a2_r_center(boxes(id)), &
               init_cond%seed_r0(:, n), &
               init_cond%seed_r1(:, n), 2)
          if (dist - init_cond%seed_width(n) < boxlen &
               .and. boxes(id)%dr > 0.2_dp * init_cond%seed_width(n)) then
             ref_flags(id) = a5_do_ref
          end if
       end do
    end if

    if (adx > 1.0_dp .and. crv_phi > 0.0_dp) ref_flags(id) = a5_do_ref
  end subroutine set_ref_flags

  subroutine set_init_cond(box)
    use m_geom
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, n, nc
    real(dp)                    :: xy(2)
    real(dp)                    :: dens

    nc = box%n_cell
    box%cc(:, :, i_elec) = init_cond%bg_dens

    do j = 0, nc+1
       do i = 0, nc+1
          xy   = a2_r_cc(box, [i,j])

          do n = 1, init_cond%n_cond
             dens = init_cond%seed_dens(n) * &
                  GM_dens_line(xy, init_cond%seed_r0(:, n), &
                  init_cond%seed_r1(:, n), 2, &
                  init_cond%seed_width(n), &
                  init_cond%seed_falloff(n))
             box%cc(i, j, i_elec) = box%cc(i, j, i_elec) + dens
          end do
       end do
    end do

    box%cc(:, :, i_pion) = box%cc(:, :, i_elec)
    box%cc(:, :, i_phi) = 0     ! Inital potential set to zero

    call set_box_eps(box)
  end subroutine set_init_cond

  subroutine set_box_eps(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell

    do j = 0, nc+1
       do i = 0, nc+1
          xy = a2_r_cc(box, [i,j]) / domain_len

          if (xy(1) < 0.25_dp) then
             box%cc(i, j, i_eps) = epsilon_diel
          else
             box%cc(i, j, i_eps) = 1.0_dp
          end if
       end do
    end do
  end subroutine set_box_eps

  ! Get maximum time step based on e.g. CFL criteria
  real(dp) function get_max_dt(tree)
    type(a2_t), intent(in) :: tree
    real(dp), parameter    :: UC_eps0        = 8.8541878176d-12
    real(dp), parameter    :: UC_elem_charge = 1.6022d-19
    real(dp)               :: max_fld, min_fld, max_dns, dr_min
    real(dp)               :: mobility, diff_coeff, alpha, max_mobility
    real(dp)               :: dt_cfl, dt_dif, dt_drt, dt_alpha

    call a2_tree_max_cc(tree, i_fld, max_fld)
    call a2_tree_min_cc(tree, i_fld, min_fld)
    call a2_tree_max_cc(tree, i_elec, max_dns)

    dr_min       = a2_min_dr(tree)
    mobility     = LT_get_col(td_tbl, i_mobility, max_fld)
    max_mobility = LT_get_col(td_tbl, i_mobility, min_fld)
    diff_coeff   = LT_get_col(td_tbl, i_diffusion, max_fld)
    alpha        = LT_get_col(td_tbl, i_alpha, max_fld)

    ! CFL condition. Note there should be a factor sqrt(0.5), but instead we
    ! rely on the CFL number 0.5 at the bottom
    dt_cfl = dr_min / (mobility * max_fld)

    ! Diffusion condition
    dt_dif = 0.25_dp * dr_min**2 / diff_coeff

    ! Dielectric relaxation time
    dt_drt = UC_eps0 / (UC_elem_charge * max_mobility * &
         max(epsilon(1.0_dp), max_dns))

    ! Ionization limit
    dt_alpha =  1 / max(mobility * max_fld * alpha, epsilon(1.0_dp))

    get_max_dt = 0.5_dp * min(1/(1/dt_cfl + 1/dt_dif), dt_alpha, dt_max)
  end function get_max_dt

  ! Compute electric field on the tree. First perform multigrid to get electric
  ! potential, then take numerical gradient to geld field.
  subroutine compute_fld(tree, n_cycles, no_guess)
    use m_units_constants
    type(a2_t), intent(inout) :: tree
    integer, intent(in)       :: n_cycles
    logical, intent(in)       :: no_guess
    real(dp), parameter       :: fac = UC_elem_charge / UC_eps0
    integer                   :: lvl, i, id, nc

    nc = tree%n_cell

    ! Set the source term (rhs)
    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          tree%boxes(id)%cc(:, :, i_rhs) = fac * (&
               tree%boxes(id)%cc(:, :, i_elec) - &
               tree%boxes(id)%cc(:, :, i_pion))
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    applied_voltage = -domain_len * get_fld(time)

    ! Perform n_cycles fmg cycles (logicals: store residual, first call)
    do i = 1, n_cycles
       call mg2_fas_fmg(tree, mg, .true., no_guess .and. i == 1)
    end do

    ! Compute field from potential
    call a2_loop_box(tree, fld_from_pot)

    ! Set the field norm also in ghost cells
    call a2_gc_tree(tree, i_fld, a2_gc_interp, a2_gc_neumann)
  end subroutine compute_fld

  real(dp) function get_fld(time)
    use m_units_constants
    real(dp), intent(in) :: time

    if (time > fld_mod_t0) then
       get_fld = applied_fld + (time - fld_mod_t0) * fld_lin_deriv + &
            fld_sin_amplitude * &
            sin((time - fld_mod_t0) * 2 * UC_pi * fld_sin_freq)
    else
       get_fld = applied_fld
    end if
  end function get_fld

  ! Compute electric field from electrical potential
  subroutine fld_from_pot(box)
    type(box2_t), intent(inout) :: box
    integer                     :: nc
    real(dp)                    :: inv_dr

    nc     = box%n_cell
    inv_dr = 1 / box%dr

    box%fx(:, :, f_fld) = inv_dr * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi))
    box%fy(:, :, f_fld) = inv_dr * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi))

    ! Compute fields at the boundaries of the box, where eps can change (have to
    ! be careful that there is enough refinement)
    box%fx(1, :, f_fld) = 2 * inv_dr * &
         (box%cc(0, 1:nc, i_phi) - box%cc(1, 1:nc, i_phi)) * &
         box%cc(0, 1:nc, i_eps) / &
         (box%cc(1, 1:nc, i_eps) + box%cc(0, 1:nc, i_eps))
    box%fx(nc+1, :, f_fld) = 2 * inv_dr * &
         (box%cc(nc, 1:nc, i_phi) - box%cc(nc+1, 1:nc, i_phi)) * &
         box%cc(nc+1, 1:nc, i_eps) / &
         (box%cc(nc+1, 1:nc, i_eps) + box%cc(nc, 1:nc, i_eps))
    box%fy(:, 1, f_fld) = 2 * inv_dr * &
         (box%cc(1:nc, 0, i_phi) - box%cc(1:nc, 1, i_phi)) * &
         box%cc(1:nc, 0, i_eps) / &
         (box%cc(1:nc, 1, i_eps) + box%cc(1:nc, 0, i_eps))
    box%fy(:, nc+1, f_fld) = 2 * inv_dr * &
         (box%cc(1:nc, nc, i_phi) - box%cc(1:nc, nc+1, i_phi)) * &
         box%cc(1:nc, nc+1, i_eps) / &
         (box%cc(1:nc, nc+1, i_eps) + box%cc(1:nc, nc, i_eps))

    box%cc(1:nc, 1:nc, i_fld) = sqrt(&
         0.25_dp * (box%fx(1:nc, 1:nc, f_fld) + box%fx(2:nc+1, 1:nc, f_fld))**2 + &
         0.25_dp * (box%fy(1:nc, 1:nc, f_fld) + box%fy(1:nc, 2:nc+1, f_fld))**2)
  end subroutine fld_from_pot

  ! This fills ghost cells near physical boundaries for the potential
  subroutine sides_bc_pot(boxes, id, nb, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    integer                     :: nc

    nc = boxes(id)%n_cell

    select case (nb)
    case (a2_nb_lx)             ! Neumann
       boxes(id)%cc(0, 1:nc, iv) = boxes(id)%cc(1, 1:nc, iv)
    case (a2_nb_hx)             ! Neumann
       boxes(id)%cc(nc+1, 1:nc, iv) = boxes(id)%cc(nc, 1:nc, iv)
    case (a2_nb_ly)             ! Grounded
       boxes(id)%cc(1:nc, 0, iv) = -boxes(id)%cc(1:nc, 1, iv)
    case (a2_nb_hy)             ! Applied voltage
       boxes(id)%cc(:, nc+1, iv) = 2 * applied_voltage &
            - boxes(id)%cc(:, nc, iv)
    end select
  end subroutine sides_bc_pot

  ! Compute the electron fluxes due to drift and diffusion
  subroutine fluxes_koren(boxes, id, dt_vec)
    use m_units_constants
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp), intent(in)        :: dt_vec(:)
    real(dp)                    :: fac, inv_dr, tmp, gradp, gradc, gradn
    real(dp)                    :: mobility, diff_coeff, v_drift
    real(dp)                    :: fld, fld_avg
    real(dp)                    :: gc_data(boxes(id)%n_cell, a2_num_neighbors)
    integer                     :: i, j, nc
    type(LT_loc_t) :: loc

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr
    fac    = -0.8_dp * UC_eps0 / (UC_elem_charge * dt_vec(1))

    call a2_gc2_box(boxes, id, i_elec, a2_gc2_prolong1, &
         a2_gc2_neumann, gc_data, nc)

    ! x-fluxes interior, advective part with flux limiter
    do j = 1, nc
       do i = 1, nc+1
          fld_avg   = 0.5_dp * (boxes(id)%cc(i, j, i_fld) + &
               boxes(id)%cc(i-1, j, i_fld))
          loc        = LT_get_loc(td_tbl, fld_avg)
          mobility   = LT_get_col_at_loc(td_tbl, i_mobility, loc)
          diff_coeff = LT_get_col_at_loc(td_tbl, i_diffusion, loc)
          fld        = boxes(id)%fx(i, j, f_fld)
          v_drift    = -mobility * fld
          gradc      = boxes(id)%cc(i, j, i_elec) - boxes(id)%cc(i-1, j, i_elec)

          if (v_drift < 0.0_dp) then
             if (i == nc+1) then
                tmp = gc_data(j, a2_nb_hx)
             else
                tmp = boxes(id)%cc(i+1, j, i_elec)
             end if
             gradn = tmp - boxes(id)%cc(i, j, i_elec)
             boxes(id)%fx(i, j, f_elec) = v_drift * &
                  (boxes(id)%cc(i, j, i_elec) - koren_mlim(gradc, gradn))
             if (boxes(id)%fx(i, j, f_elec) < fac * fld) &
                  boxes(id)%fx(i, j, f_elec) = fac * fld
          else                  ! v_drift > 0
             if (i == 1) then
                tmp = gc_data(j, a2_nb_lx)
             else
                tmp = boxes(id)%cc(i-2, j, i_elec)
             end if
             gradp = boxes(id)%cc(i-1, j, i_elec) - tmp
             boxes(id)%fx(i, j, f_elec) = v_drift * &
                  (boxes(id)%cc(i-1, j, i_elec) + koren_mlim(gradc, gradp))
             if (boxes(id)%fx(i, j, f_elec) > fac * fld) &
                  boxes(id)%fx(i, j, f_elec) = fac * fld
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be
          ! scaled by 1/dx
          boxes(id)%fx(i, j, f_elec) = boxes(id)%fx(i, j, f_elec) - &
               diff_coeff * gradc * inv_dr
       end do
    end do

    ! y-fluxes interior, advective part with flux limiter
    do j = 1, nc+1
       do i = 1, nc
          fld_avg    = 0.5_dp * (boxes(id)%cc(i, j, i_fld) + &
               boxes(id)%cc(i, j-1, i_fld))
          loc        = LT_get_loc(td_tbl, fld_avg)
          mobility   = LT_get_col_at_loc(td_tbl, i_mobility, loc)
          diff_coeff = LT_get_col_at_loc(td_tbl, i_diffusion, loc)
          fld        = boxes(id)%fy(i, j, f_fld)
          v_drift    = -mobility * fld
          gradc      = boxes(id)%cc(i, j, i_elec) - boxes(id)%cc(i, j-1, i_elec)

          if (v_drift < 0.0_dp) then
             if (j == nc+1) then
                tmp = gc_data(i, a2_nb_hy)
             else
                tmp = boxes(id)%cc(i, j+1, i_elec)
             end if
             gradn = tmp - boxes(id)%cc(i, j, i_elec)
             boxes(id)%fy(i, j, f_elec) = v_drift * &
                  (boxes(id)%cc(i, j, i_elec) - koren_mlim(gradc, gradn))
             if (boxes(id)%fy(i, j, f_elec) < fac * fld) &
                  boxes(id)%fy(i, j, f_elec) = fac * fld
          else                  ! v_drift > 0
             if (j == 1) then
                tmp = gc_data(i, a2_nb_ly)
             else
                tmp = boxes(id)%cc(i, j-2, i_elec)
             end if
             gradp = boxes(id)%cc(i, j-1, i_elec) - tmp
             boxes(id)%fy(i, j, f_elec) = v_drift * &
                  (boxes(id)%cc(i, j-1, i_elec) + koren_mlim(gradc, gradp))
             if (boxes(id)%fy(i, j, f_elec) > fac * fld) &
                  boxes(id)%fy(i, j, f_elec) = fac * fld
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be
          ! scaled by 1/dx
          boxes(id)%fy(i, j, f_elec) = boxes(id)%fy(i, j, f_elec) - &
               diff_coeff * gradc * inv_dr
       end do
    end do

  end subroutine fluxes_koren

  ! Take average of new and old electron/ion density for explicit trapezoidal rule
  subroutine average_dens(box)
    type(box2_t), intent(inout) :: box
    box%cc(:, :, i_elec) = 0.5_dp * (box%cc(:, :, i_elec) + box%cc(:, :, i_elec_old))
    box%cc(:, :, i_pion) = 0.5_dp * (box%cc(:, :, i_pion) + box%cc(:, :, i_pion_old))
  end subroutine average_dens

  ! Advance solution over dt based on the fluxes / source term, using forward Euler
  subroutine update_solution(box, dt)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr, src, sflux, fld
    real(dp)                    :: alpha, eta, dflux(2), rfac(2)
    integer                     :: i, j, nc, ioff
    type(LT_loc_t)              :: loc

    nc     = box%n_cell
    inv_dr = 1/box%dr
    ioff   = (box%ix(1)-1) * nc

    do j = 1, nc
       do i = 1, nc
          ! Weighting of flux contribution for cylindrical coordinates
          rfac = [i+ioff-1, i+ioff] / (i+ioff-0.5_dp)

          fld      = box%cc(i,j, i_fld)
          loc      = LT_get_loc(td_tbl, fld)
          alpha    = LT_get_col_at_loc(td_tbl, i_alpha, loc)
          eta      = LT_get_col_at_loc(td_tbl, i_eta, loc)

          ! Set source term equal to ||flux|| * (alpha - eta)
          dflux(1) = box%fx(i, j, f_elec) + box%fx(i+1, j, f_elec)
          dflux(2) = box%fy(i, j, f_elec) + box%fy(i, j+1, f_elec)
          src = 0.5_dp * norm2(dflux) * (alpha - eta)

          if (photoi_enabled) &
               src = src + box%cc(i,j, i_pho)

          ! Contribution of flux
          sflux = (box%fy(i, j, f_elec) - box%fy(i, j+1, f_elec) + &
               rfac(1) * box%fx(i, j, f_elec) - &
               rfac(2) * box%fx(i+1, j, f_elec)) * inv_dr

          box%cc(i, j, i_elec) = box%cc(i, j, i_elec) + (src + sflux) * dt(1)
          box%cc(i, j, i_pion) = box%cc(i, j, i_pion) + src * dt(1)
       end do
    end do
  end subroutine update_solution

  subroutine set_photoionization(tree, eta, num_photons, dt)
    use m_units_constants

    type(a2_t), intent(inout) :: tree
    real(dp), intent(in)      :: eta
    real(dp), intent(in), optional :: dt
    integer, intent(in)       :: num_photons
    real(dp), parameter       :: p_quench = 30.0D0 * UC_torr_to_bar
    real(dp)                  :: quench_fac

    ! Compute quench factor, because some excited species will be quenched by
    ! collisions, preventing the emission of a UV photon
    quench_fac = p_quench / (gas_pressure + p_quench)

    ! Set photon production rate per cell, which is proportional to the
    ! ionization rate.
    call a2_loop_box_arg(tree, set_photoi_rate, [eta * quench_fac], .true.)

    call PH_set_src_2d(tree, photoi_tbl, sim_rng, num_photons, &
         i_pho, i_pho, 0.25e-3_dp, .true., .true., 1e-9_dp, dt)

  end subroutine set_photoionization

  subroutine set_photoi_rate(box, coeff)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: coeff(:)
    integer                     :: i, j, nc
    real(dp)                    :: fld, alpha, mobility, dr, tmp
    type(LT_loc_t)              :: loc

    nc = box%n_cell

    do j = 1, nc
       do i = 1, nc
          dr       = box%dr
          fld      = box%cc(i, j, i_fld)
          loc      = LT_get_loc(td_tbl, fld)
          alpha    = LT_get_col_at_loc(td_tbl, i_alpha, loc)
          mobility = LT_get_col_at_loc(td_tbl, i_mobility, loc)

          tmp = fld * mobility * alpha * box%cc(i, j, i_elec) * coeff(1)
          if (tmp < 0) tmp = 0
          box%cc(i, j, i_pho) = tmp
       end do
    end do
  end subroutine set_photoi_rate

  ! For each box that gets refined, set data on its children using this routine
  subroutine prolong_to_new_boxes(tree, ref_info)
    use m_a2_prolong
    type(a2_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id

    do lvl = 1, tree%max_lvl
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a2_prolong1_to(tree%boxes, id, i_elec)
          call a2_prolong1_to(tree%boxes, id, i_pion)
          call a2_prolong1_to(tree%boxes, id, i_phi)
          call set_box_eps(tree%boxes(id))
       end do

       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a2_gc_box(tree%boxes, id, i_elec, &
               a2_gc_interp_lim, a2_gc_neumann)
          call a2_gc_box(tree%boxes, id, i_pion, &
               a2_gc_interp_lim, a2_gc_neumann)
          call a2_gc_box(tree%boxes, id, i_phi, &
               mg2_sides_rb, sides_bc_pot)
       end do
    end do
  end subroutine prolong_to_new_boxes

  subroutine write_streamer_properties(tree, fname_stats, fname_axis, output_cnt, time)
    type(a2_t), intent(in)       :: tree
    character(len=*), intent(in) :: fname_axis, fname_stats
    integer, intent(in) :: output_cnt
    real(dp), intent(in) :: time

    real(dp), allocatable        :: props(:), axis_data(:,:)
    character(len=15), allocatable :: prop_names(:)
    integer                      :: n
    integer, parameter           :: unit_1 = 777, unit_2 = 778

    call get_streamer_properties(tree, props, prop_names)
    call get_cc_axis(tree, [i_elec, i_pion], [f_fld, f_elec], axis_data)

    if (output_cnt == 1) then
       open(unit_1, file=trim(fname_stats), action="write")
       write(unit_1, *) "ix_out      time     ", prop_names
       close(unit_1)
    else
       open(unit_1, file=trim(fname_stats), action="write", &
            position="append")
       write(unit_1, *) output_cnt, time, props
       close(unit_1)
    end if

    open(unit_2, file=trim(fname_axis), action="write")
    write(unit_1, *) "z, n_e, n_i, fld_z, flux_z"
    do n = 1, size(axis_data, 2)
       write(unit_2, *) axis_data(:, n)
    end do
    close(unit_2)
    deallocate(axis_data)

    print *, "Written ", trim(fname_axis)
  end subroutine write_streamer_properties

  subroutine get_streamer_properties(tree, props, prop_names)
    type(a2_t), intent(in)                        :: tree
    real(dp), intent(inout), allocatable          :: props(:)
    character(len=15), intent(inout), allocatable :: prop_names(:)

    integer :: ip
    integer, parameter :: n_props = 38
    real(dp)       :: rz(2), phi_head, z_head
    real(dp)       :: alpha, mu, Er_max_norm, Ez_max
    type(a2_loc_t) :: loc_ez, loc_er, loc_dens, loc
    integer        :: id, ix(2)

    allocate(props(n_props))
    allocate(prop_names(n_props))

    ip = 1
    prop_names(ip) = "E_bg"
    props(ip)      = get_fld(time)

    ip             = ip + 1
    prop_names(ip) = "Ez_max"
    call a2_reduction_loc(tree, box_maxfld_z, reduce_max, &
         -1.0e99_dp, props(ip), loc_ez)
    Ez_max = props(ip)

    ip             = ip + 1
    prop_names(ip) = "Er_max"
    call a2_reduction_loc(tree, box_maxfld_r, reduce_max, &
         -1.0e99_dp, props(ip), loc_er)

    ! Radius of streamer is defined as location of maximum r-field
    rz       = a2_r_loc(tree, loc_er)

    ip             = ip + 1
    prop_names(ip) = "Er_max(r)"
    props(ip)      = rz(1)

    ip             = ip + 1
    prop_names(ip) = "Er_max(z)"
    props(ip)      = rz(2)
    z_head         = rz(2)

    ip             = ip + 1
    prop_names(ip) = "Er_max(E)"
    props(ip)      = tree%boxes(loc_er%id)%cc(loc_er%ix(1), &
         loc_er%ix(2), i_fld)
    Er_max_norm    = props(ip)

    alpha = LT_get_col(td_tbl, i_alpha, Ez_max)
    mu = LT_get_col(td_tbl, i_alpha, Ez_max)
    ip             = ip + 1
    prop_names(ip) = "Ez_max(alpha)"
    props(ip)      = alpha
    ip             = ip + 1
    prop_names(ip) = "Ez_max(S)"
    props(ip)      = alpha * mu * Ez_max

    alpha = LT_get_col(td_tbl, i_alpha, Er_max_norm)
    mu = LT_get_col(td_tbl, i_alpha, Er_max_norm)
    ip             = ip + 1
    prop_names(ip) = "Er_max(alpha)"
    props(ip)      = alpha
    ip             = ip + 1
    prop_names(ip) = "Er_max(S)"
    props(ip)      = alpha * mu * Er_max_norm

    if (time > 1.0e-9_dp .and. rz(2) > 0.9_dp * domain_len .or. &
         rz(2) < 0.1_dp * domain_len) &
         stop "Simulation has reached boundary"

    ! Get electron density and potential at location of radius
    loc_dens       = a2_get_loc(tree, [0.0_dp, rz(2)])
    id             = loc_dens%id
    ix             = loc_dens%ix
    ip             = ip + 1
    prop_names(ip) = "n_e(head)"
    props(ip)      = tree%boxes(id)%cc(ix(1), ix(2), i_elec)

    ! Set phi to potential difference
    phi_head       = tree%boxes(id)%cc(ix(1), ix(2), i_phi)
    ip             = ip + 1
    prop_names(ip) = "dphi"
    props(ip)      = phi_head - (rz(2)/domain_len) * applied_voltage

    ! Height of streamer
    rz             = a2_r_cc(tree%boxes(loc_ez%id), loc_ez%ix)
    ip             = ip + 1
    prop_names(ip) = "Ez_max(z)"
    props(ip)      = rz(2)

    do i = 4, 14, 2
       GLOBAL_fld_val = i * 1e6_dp
       ip = ip + 1
       write(prop_names(ip), "(A,I0,A)") "r_max(", i, "e6)"
       call a2_reduction_loc(tree, box_fld_maxr, reduce_max, &
            0.0_dp, props(ip), loc)

       ip = ip + 1
       write(prop_names(ip), "(A,I0,A)") "dphi_r(", i, "e6)"
       if (loc%id > a5_no_box) then
          props(ip) = phi_head - &
               tree%boxes(loc%id)%cc(loc%ix(1), loc%ix(2), i_phi)
       else
          props(ip) = 0
       end if

       ip = ip + 1
       write(prop_names(ip), "(A,I0,A)") "z_max(", i, "e6)"
       call a2_reduction_loc(tree, box_fld_maxz, reduce_max, &
            0.0_dp, props(ip), loc)
       props(ip) = props(ip) - z_head

       ip = ip + 1
       write(prop_names(ip), "(A,I0,A)") "dphi_z(", i, "e6)"
       if (loc%id > a5_no_box) then
          props(ip) = phi_head - &
               tree%boxes(loc%id)%cc(loc%ix(1), loc%ix(2), i_phi)

       else
          props(ip) = 0
       end if
    end do

    print *, "n_props", ip
  end subroutine get_streamer_properties

  real(dp) function reduce_max(a, b)
    real(dp), intent(in) :: a, b
    reduce_max = max(a,b)
  end function reduce_max

  subroutine box_maxfld_r(box, val, ix)
    type(box2_t), intent(in) :: box
    real(dp), intent(out)    :: val
    integer, intent(out)     :: ix(2)
    integer                  :: nc

    nc = box%n_cell
    ix = maxloc(abs(box%fx(1:nc, :, f_fld) + box%fx(2:nc+1, :, f_fld)))
    val = 0.5_dp * abs(box%fx(ix(1), ix(2), f_fld) + &
         box%fx(ix(1)+1, ix(2), f_fld))
  end subroutine box_maxfld_r

  subroutine box_fld_maxr(box, val, ix)
    type(box2_t), intent(in) :: box
    real(dp), intent(out)    :: val
    integer, intent(out)     :: ix(2)
    integer                  :: i, j, nc
    real(dp)                 :: rz(2)

    val = 0
    nc = box%n_cell

    do j = 1, nc
       do i = 1, nc
          if (box%cc(i, j, i_fld) > GLOBAL_fld_val) then
             rz = a2_r_cc(box, [i,j])
             if (rz(1) > val) then
                val = rz(1)
                ix = [i, j]
             end if
          end if
       end do
    end do
  end subroutine box_fld_maxr

  subroutine box_fld_maxz(box, val, ix)
    type(box2_t), intent(in) :: box
    real(dp), intent(out)    :: val
    integer, intent(out)     :: ix(2)
    integer                  :: i, j, nc
    real(dp)                 :: rz(2)

    val = 0
    nc = box%n_cell

    do j = 1, nc
       do i = 1, nc
          if (box%cc(i, j, i_fld) > GLOBAL_fld_val) then
             rz = a2_r_cc(box, [i,j])
             if (rz(2) > val) then
                val = rz(2)
                ix = [i, j]
             end if
          end if
       end do
    end do
  end subroutine box_fld_maxz

  subroutine box_maxfld_z(box, val, ix)
    type(box2_t), intent(in) :: box
    real(dp), intent(out)    :: val
    integer, intent(out)     :: ix(2)
    integer                  :: nc

    nc = box%n_cell
    ix = maxloc(abs(box%fy(:, 1:nc, f_fld) + box%fy(:, 2:nc+1, f_fld)))
    val = 0.5_dp * abs(box%fy(ix(1), ix(2), f_fld) + &
         box%fy(ix(1), ix(2)+1, f_fld))
  end subroutine box_maxfld_z

  subroutine get_cc_axis(tree, ixs_cc, ixs_fc, axis_data)
    type(a2_t), intent(in)               :: tree
    integer, intent(in)                  :: ixs_cc(:), ixs_fc(:)
    real(dp), allocatable, intent(inout) :: axis_data(:, :)

    type(a2_loc_t)                       :: loc
    real(dp)                             :: z, box_z, box_dz
    integer                              :: i, id, nc, cnt, n_cc, n_fc

    n_cc = size(ixs_cc)
    n_fc = size(ixs_fc)
    nc   = tree%n_cell
    z    = 0
    cnt  = 0

    ! Determine how many boxes lie on the axis
    do
       loc = a2_get_loc(tree, [0.0_dp, z])
       if (loc%id == -1) exit

       cnt    = cnt + nc
       id     = loc%id
       box_z  = tree%boxes(id)%r_min(2)
       box_dz = tree%boxes(id)%dr
       z      = box_z + (nc+1) * box_dz
    end do

    ! Now store the actual axis data
    allocate(axis_data(n_cc+n_fc+1, cnt))
    cnt = 0
    z   = 0

    do
       loc = a2_get_loc(tree, [0.0_dp, z])
       if (loc%id == -1) exit

       id     = loc%id
       box_z  = tree%boxes(id)%r_min(2)
       box_dz = tree%boxes(id)%dr

       axis_data(1, cnt+1:cnt+nc) = &
            box_z + [((i-0.5_dp) * box_dz, i = 1, nc)]
       axis_data(2:n_cc+1, cnt+1:cnt+nc) = &
            transpose(tree%boxes(id)%cc(1, 1:nc, ixs_cc))
       axis_data(n_cc+2:, cnt+1:cnt+nc) = transpose( &
            0.5_dp * (tree%boxes(id)%fy(1, 1:nc, ixs_fc) + &
            tree%boxes(id)%fy(1, 2:nc+1, ixs_fc)))

       cnt    = cnt + nc
       z      = box_z + (nc+1) * box_dz
    end do
  end subroutine get_cc_axis

end program
