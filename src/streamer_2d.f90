program streamer_2d

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

  integer                :: i, n
  character(len=ST_slen) :: fname
  logical                :: write_out

  type(a2_t)             :: tree ! This contains the full grid information
  type(mg2_t)            :: mg   ! Multigrid option struct
  type(ref_info_t)       :: ref_info

  call ST_create_cfg(ST_cfg)
  call ST_read_cfg_files(ST_cfg)
  call ST_load_cfg(ST_cfg)

  ! Initialize the transport coefficients
  call ST_load_transport_data(ST_cfg, ST_td_tbl, ST_photoi_tbl)

  ! Set the initial conditions from the configuration
  call ST_get_init_cond(ST_cfg, ST_init_cond, 2)

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the multigrid options. First define the variables to use
  mg%i_phi        = i_phi
  mg%i_tmp        = i_fld
  mg%i_rhs        = i_rhs
  mg%i_eps        = i_eps

  ! The number of cycles at the lowest level
  mg%n_cycle_base = 8

  ! Routines to use for ...
  mg%sides_bc    => sides_bc_pot ! Filling ghost cell on physical boundaries
  mg%box_op      => mg2_auto_op
  mg%box_corr    => mg2_auto_corr
  mg%box_gsrb    => mg2_auto_gsrb

  ! This routine always needs to be called when using multigrid
  call mg2_init_mg(mg)

  ST_out_cnt = 0 ! Number of output files written
  ST_time    = 0 ! Simulation time (all times are in s)

  ! Set up the initial conditions
  do
     call a2_loop_box(tree, set_init_cond)
     call compute_fld(tree, n_fmg_cycles, .true.)
     call a2_adjust_refinement(tree, set_ref_flags, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  if (ST_photoi_enabled) &
       call set_photoionization(tree, ST_photoi_eta, ST_photoi_num_photons)

  do
     ! Get a new time step, which is at most dt_amr
     ST_dt = get_max_dt(tree)

     if (ST_dt < 1e-14) then
        print *, "ST_dt getting too small, instability?"
        ST_time = ST_end_time + 1.0_dp
     end if

     ! Every ST_dt_out, write output
     if (ST_out_cnt * ST_dt_out <= ST_time) then
        write_out = .true.
        ST_out_cnt = ST_out_cnt + 1
        write(fname, "(A,I6.6)") trim(ST_sim_name) // "_", ST_out_cnt
     else
        write_out = .false.
     end if

     if (write_out) call a2_write_silo(tree, fname, &
          cc_names, ST_out_cnt, ST_time, dir=ST_output_dir)

     if (ST_time > ST_end_time) exit

     ! We perform n_steps between mesh-refinements
     do n = 1, ST_steps_amr
        ST_time = ST_time + ST_dt

        if (ST_photoi_enabled) &
             call set_photoionization(tree, ST_photoi_eta, ST_photoi_num_photons, ST_dt)

        ! Copy previous solution
        call a2_tree_copy_cc(tree, i_elec, i_elec_old)
        call a2_tree_copy_cc(tree, i_pion, i_pion_old)

        ! Two forward Euler steps over ST_dt
        do i = 1, 2
           ! First calculate fluxes
           call a2_loop_boxes_arg(tree, fluxes_koren, [ST_dt], .true.)
           call a2_consistent_fluxes(tree, [f_elec])

           ! Update the solution
           call a2_loop_box_arg(tree, update_solution, [ST_dt], .true.)

           ! Restrict the electron and ion densities to lower levels
           call a2_restrict_tree(tree, i_elec)
           call a2_restrict_tree(tree, i_pion)

           ! Fill ghost cells
           call a2_gc_tree(tree, i_elec, a2_gc_interp_lim, a2_gc_neumann)
           call a2_gc_tree(tree, i_pion, a2_gc_interp_lim, a2_gc_neumann)

           ! Compute new field on first iteration
           if (i == 1) call compute_fld(tree, n_fmg_cycles, .false.)
        end do

        ! Take average of phi_old and phi (explicit trapezoidal rule)
        call a2_loop_box(tree, average_dens)

        ! Compute field with new density
        call compute_fld(tree, n_fmg_cycles, .false.)
     end do

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

    dr = ST_domain_len / ST_box_size

    ! Initialize tree
    call a2_init(tree, ST_box_size, n_var_cell, n_var_face, dr, &
         coarsen_to=2, n_boxes = n_boxes_init)

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
    alpha     = LT_get_col(ST_td_tbl, i_alpha, max_fld)
    adx       = boxes(id)%dr * alpha

    if (adx < 0.1_dp .and. boxes(id)%dr < 2.0e-5_dp) &
         ref_flags(id) = a5_rm_ref
    if (adx < 0.1_dp .and. crv_phi < 4.0_dp .and. boxes(id)%dr < 5.0e-5_dp) &
         ref_flags(id) = a5_rm_ref

    ! Refine around initial conditions
    if (ST_time < 2.0e-9_dp) then
       boxlen = boxes(id)%n_cell * boxes(id)%dr

       do n = 1, ST_init_cond%n_cond
          dist = GM_dist_line(a2_r_center(boxes(id)), &
               ST_init_cond%seed_r0(:, n), &
               ST_init_cond%seed_r1(:, n), 2)
          if (dist - ST_init_cond%seed_width(n) < boxlen &
               .and. boxes(id)%dr > 0.2_dp * ST_init_cond%seed_width(n)) then
             ref_flags(id) = a5_do_ref
          end if
       end do
    end if

    if (adx > 1.0_dp .and. crv_phi > 1) ref_flags(id) = a5_do_ref

  end subroutine set_ref_flags

  subroutine set_init_cond(box)
    use m_geom
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, n, nc
    real(dp)                    :: xy(2)
    real(dp)                    :: dens

    nc = box%n_cell
    box%cc(:, :, i_elec) = ST_init_cond%bg_dens


    do j = 0, nc+1
       do i = 0, nc+1
          xy   = a2_r_cc(box, [i,j])

          do n = 1, ST_init_cond%n_cond
             dens = ST_init_cond%seed_dens(n) * &
                  GM_dens_line(xy, ST_init_cond%seed_r0(:, n), &
                  ST_init_cond%seed_r1(:, n), 2, &
                  ST_init_cond%seed_width(n), &
                  ST_init_cond%seed_falloff(n))
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
          xy = a2_r_cc(box, [i,j]) / ST_domain_len

          if (xy(2) < 0.25_dp) then
             box%cc(i, j, i_eps) = ST_epsilon_diel
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
    mobility     = LT_get_col(ST_td_tbl, i_mobility, max_fld)
    max_mobility = LT_get_col(ST_td_tbl, i_mobility, min_fld)
    diff_coeff   = LT_get_col(ST_td_tbl, i_diffusion, max_fld)
    alpha        = LT_get_col(ST_td_tbl, i_alpha, max_fld)

    ! CFL condition
    dt_cfl = dr_min / (mobility * max_fld) ! Factor ~ sqrt(0.5)

    ! Diffusion condition
    dt_dif = 0.25_dp * dr_min**2 / diff_coeff

    ! Dielectric relaxation time
    dt_drt = UC_eps0 / (UC_elem_charge * max_mobility * max_dns)

    ! Ionization limit
    dt_alpha =  1 / max(mobility * max_fld * alpha, epsilon(1.0_dp))

    get_max_dt = 0.5_dp * min(1/(1/dt_cfl + 1/dt_dif), &
         dt_alpha, ST_dt_max)
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

    ST_applied_voltage = -ST_domain_len * ST_get_fld(ST_time)

    ! Perform n_cycles fmg cycles (logicals: store residual, first call)
    do i = 1, n_cycles
       call mg2_fas_fmg(tree, mg, .false., no_guess .and. i == 1)
    end do

    ! Compute field from potential
    call a2_loop_box(tree, fld_from_pot)

    ! Set the field norm also in ghost cells
    call a2_gc_tree(tree, i_fld, a2_gc_interp, a2_gc_neumann)
  end subroutine compute_fld

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
          loc        = LT_get_loc(ST_td_tbl, fld_avg)
          mobility   = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)
          diff_coeff = LT_get_col_at_loc(ST_td_tbl, i_diffusion, loc)
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
          loc        = LT_get_loc(ST_td_tbl, fld_avg)
          mobility   = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)
          diff_coeff = LT_get_col_at_loc(ST_td_tbl, i_diffusion, loc)
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
    real(dp)                    :: alpha, eta, dflux(2)
    integer                     :: i, j, nc
    type(LT_loc_t) :: loc

    nc                    = box%n_cell
    inv_dr                = 1/box%dr
    do j = 1, nc
       do i = 1, nc
          fld      = box%cc(i,j, i_fld)
          loc      = LT_get_loc(ST_td_tbl, fld)
          alpha    = LT_get_col_at_loc(ST_td_tbl, i_alpha, loc)
          eta      = LT_get_col_at_loc(ST_td_tbl, i_eta, loc)
          ! mobility = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)
          ! src = abs(mobility * fld) * (alpha-eta) * &
          !      box%cc(i, j, i_elec)

          dflux(1) = box%fx(i, j, f_elec) + box%fx(i+1, j, f_elec)
          dflux(2) = box%fy(i, j, f_elec) + box%fy(i, j+1, f_elec)
          src = 0.5_dp * norm2(dflux) * (alpha - eta)

          if (ST_photoi_enabled) &
               src = src + box%cc(i,j, i_pho)

          sflux = (box%fx(i, j, f_elec) - box%fx(i+1, j, f_elec) + &
               box%fy(i, j, f_elec) - box%fy(i, j+1, f_elec)) * inv_dr

          box%cc(i, j, i_elec) = box%cc(i, j, i_elec) + (src + sflux) * dt(1)
          box%cc(i, j, i_pion) = box%cc(i, j, i_pion) + src * dt(1)
       end do
    end do

    do j = 1, nc
       do i = 1, nc
          if (box%cc(i, j, i_eps) > 1.0_dp) then
             box%cc(i, j, i_pion) = box%cc(i, j, i_pion) - box%cc(i, j, i_elec)
             box%cc(i, j, i_elec) = 0
          end if
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
    quench_fac = p_quench / (ST_gas_pressure + p_quench)

    ! Set photon production rate per cell, which is proportional to the
    ! ionization rate.
    call a2_loop_box_arg(tree, set_photoi_rate, [eta * quench_fac], .true.)

    call PH_set_src_2d(tree, ST_photoi_tbl, ST_rng, num_photons, &
         i_pho, i_pho, 0.25e-3_dp, .false., .false., 1e-9_dp, dt)

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
          loc      = LT_get_loc(ST_td_tbl, fld)
          alpha    = LT_get_col_at_loc(ST_td_tbl, i_alpha, loc)
          mobility = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)

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

  ! This fills ghost cells near physical boundaries for the potential
  subroutine sides_bc_pot(boxes, id, nb, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    integer                     :: nc

    nc = boxes(id)%n_cell

    select case (nb)
    case (a2_nb_lx)
       boxes(id)%cc(0, 1:nc, iv) = boxes(id)%cc(1, 1:nc, iv)
    case (a2_nb_hx)
       boxes(id)%cc(nc+1, 1:nc, iv) = boxes(id)%cc(nc, 1:nc, iv)
    case (a2_nb_ly)             ! Grounded
       boxes(id)%cc(1:nc, 0, iv) = -boxes(id)%cc(1:nc, 1, iv)
    case (a2_nb_hy)
       boxes(id)%cc(1:nc, nc+1, iv) = 2 * ST_applied_voltage &
            - boxes(id)%cc(1:nc, nc, iv)
    end select
  end subroutine sides_bc_pot

end program streamer_2d
