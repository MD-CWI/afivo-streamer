!> Program to perform 3d streamer simulations in Cartesian and cylindrical coordinates
program streamer_3d

  use m_a3_all
  use m_streamer
  use m_field_3d
  use m_init_cond_3d
  use m_photoi_3d
  use m_flux_schemes_3d
  use m_units_constants

  implicit none
  
  integer, parameter     :: int8 = selected_int_kind(18)
  integer(int8)          :: t_start, t_current, count_rate
  real(dp)               :: wc_time, inv_count_rate, time_last_print
  integer                :: i, it
  character(len=ST_slen) :: fname
  logical                :: write_out
  real(dp)               :: dt_prev, med
  real(dp), parameter    :: fac = UC_elem_charge / UC_eps0
  
  

  type(CFG_t)            :: cfg ! The configuration for the simulation
  type(a3_t)            :: tree      ! This contains the full grid information
  type(mg3_t)           :: mg        ! Multigrid option struct
  type(ref_info_t)       :: ref_info

  integer :: output_cnt = 0 ! Number of output files written

  call CFG_update_from_arguments(cfg)
  call ST_initialize(cfg, 3)
  call photoi_initialize(cfg)

  call ST_load_transport_data(cfg)
  call field_initialize(cfg, mg)
  call init_cond_initialize(cfg, 3)

  fname = trim(ST_output_dir) // "/" // trim(ST_simulation_name) // "_out.cfg"
  call CFG_write(cfg, trim(fname))

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the multigrid options. First define the variables to use
  mg%i_phi     = i_phi
  mg%i_tmp     = i_electric_fld
  mg%i_rhs     = i_rhs
  mg%i_eps     = i_eps
  mg%sigma_rhs = sigma_rhs
  

  ! This automatically handles cylindrical symmetry
  mg%box_op => mg3_auto_op
  mg%box_gsrb => mg3_auto_gsrb
  mg%box_corr => mg3_auto_corr

  ! This routine always needs to be called when using multigrid
  call mg3_init_mg(mg)

  output_cnt      = 0         ! Number of output files written
  ST_time         = 0         ! Simulation time (all times are in s)
  med = a3_harm_w(1.0_dp, ST_epsilon_die, 0.5_dp)
  

  ! Set up the initial conditions
  do
     call a3_loop_box(tree, init_cond_set_box)
     call field_compute(tree, mg, .false.)
     call a3_adjust_refinement(tree, refine_routine, ref_info, &
          ST_refine_buffer_width, .true.)
     if (ref_info%n_add == 0) exit
  end do

  call init_cond_stochastic_density(tree)

  print *, "Number of threads", af_get_max_threads()
  call a3_print_info(tree)

  ! Start from small time step
  ST_dt   = ST_dt_min
  dt_prev = ST_dt
  
  ! Initial wall clock time
  call system_clock(t_start, count_rate)
  inv_count_rate = 1.0_dp / count_rate
  time_last_print = -1e10_dp

  

  do it = 1, huge(1)-1
     if (ST_time >= ST_end_time) exit
     
     ! Update wall clock time
     call system_clock(t_current)
     wc_time = (t_current - t_start) * inv_count_rate

     ! Every ST_print_status_interval, print some info about progress
     if (wc_time - time_last_print > ST_print_status_sec) then
        call print_status()
        time_last_print = wc_time
     end if

     ! Every ST_dt_output, write output
     if (output_cnt * ST_dt_output <= ST_time + ST_dt) then
        write_out  = .true.
        ST_dt      = output_cnt * ST_dt_output - ST_time
        output_cnt = output_cnt + 1
     else
        write_out = .false.
     end if

     if (photoi_enabled .and. mod(it, photoi_per_steps) == 0) then
        call photoi_set_src(tree, ST_dt)
     end if

     ! Copy previous solution
     call a3_tree_copy_cc(tree, i_electron, i_electron_old)
     call a3_tree_copy_cc(tree, i_pos_ion, i_pos_ion_old)
     

     ST_dt_matrix = ST_dt_max      ! Maximum time step
     ! Two forward Euler steps over ST_dt
     do i = 1, 2
        ST_time = ST_time + ST_dt
        

        ! First calculate fluxes
        call a3_loop_boxes(tree, fluxes_koren, .true.)
        call a3_consistent_fluxes(tree, [flux_elec])
        
        ! Update the solution
        call a3_loop_box_arg(tree, update_solution, [ST_dt, real(i, dp)], .true.)
        
        if (i == 1) then
           ! Compute new field on first iteration
           call field_compute(tree, mg, .true.)

           ! Coarse densities might be required for ghost cells
           call a3_restrict_tree(tree, i_electron, i_eps = i_eps, med = med)
           call a3_restrict_tree(tree, sigma_rhs, i_eps = i_eps, med = med, s_flag = .true.)
        end if
     end do
     

     ST_time = ST_time - ST_dt        ! Go back one time step
     

     ! Take average of n and sigma (explicit trapezoidal rule)
     call a3_loop_box(tree, average_density)

     ! Coarse densities might be required for ghost cells
     call a3_restrict_tree(tree, i_electron, i_eps = i_eps, med = med)
     call a3_restrict_tree(tree, sigma_rhs, i_eps = i_eps, med = med, s_flag = .true.)

     ! Compute field with new density
     call field_compute(tree, mg, .true.)

     ! Determine next time step
     ST_dt   = min(2 * dt_prev, ST_dt_safety_factor * minval(ST_dt_matrix))
     dt_prev = ST_dt

     if (ST_dt < ST_dt_min) then
        print *, "ST_dt getting too small, instability?", ST_dt
        error stop
     end if
     

     if (write_out) then
        ! Fill ghost cells before writing output
        call a3_gc_tree(tree, i_eps, i_electron, a3_gc_interp, a3_bc_dirichlet_mirror)
        call a3_gc_tree(tree, i_eps, i_pos_ion, a3_gc_interp, a3_bc_dirichlet_mirror)

        if (ST_output_src_term) then
           call a3_restrict_tree(tree, i_src, i_eps = i_eps, med = med)
           call a3_gc_tree(tree, i_eps, i_src, a3_gc_interp, a3_bc_neumann_zero)
        end if

        if (photoi_enabled) then
           call photoi_set_src(tree, ST_dt) ! Because the mesh could have changed
           call a3_restrict_tree(tree, i_photo, i_eps = i_eps, med = med)
           call a3_gc_tree(tree, i_eps, i_photo, a3_gc_interp, a3_bc_neumann_zero)
        end if

        write(fname, "(A,I6.6)") trim(ST_simulation_name) // "_", output_cnt
        call a3_write_silo(tree, fname, output_cnt, ST_time, &
             vars_for_output, dir=ST_output_dir)

        if (ST_datfile_write) then
           call a3_write_tree(tree, fname, ST_output_dir)
        end if

        write(fname, "(A,I6.6)") trim(ST_simulation_name) // "_log.txt"
        call write_log_file(tree, fname, output_cnt, ST_output_dir)

        if (ST_plane_write) then
           write(fname, "(A,I6.6)") trim(ST_simulation_name) // &
                "_plane_", output_cnt
           call a3_write_plane(tree, fname, [ST_plane_ivar], &
                ST_plane_rmin * ST_domain_len, &
                ST_plane_rmax * ST_domain_len, &
                ST_plane_npixels, ST_output_dir)
        end if

        if (ST_lineout_write) then
           write(fname, "(A,I6.6)") trim(ST_simulation_name) // &
                "_line_", output_cnt
           call a3_write_line(tree, trim(fname), &
                [i_electron, i_pos_ion, i_phi, i_electric_fld, i_eps], &
                r_min = ST_lineout_rmin(1:3) * ST_domain_len, &
                r_max = ST_lineout_rmax(1:3) * ST_domain_len, &
                n_points=ST_lineout_npoints, dir=ST_output_dir)
        end if
     end if
     

     if (mod(it, ST_refine_per_steps) == 0) then
        ! Restrict electron and ion densities before refinement
        call a3_restrict_tree(tree, i_electron, i_eps = i_eps, med = med)
        call a3_restrict_tree(tree, i_pos_ion, i_eps = i_eps, med = med)
        call a3_restrict_tree(tree, sigma_rhs, i_eps = i_eps, med = med, s_flag = .true.)

        if (ST_output_src_term .and. ST_output_src_decay_rate > 0) then
           ! Have to set src term on coarse grids as well, and fill ghost cells
           ! before prolongation
           call a3_restrict_tree(tree, i_src, i_eps = i_eps, med = med)
           call a3_gc_tree(tree, i_eps, i_src, a3_gc_interp, a3_bc_neumann_zero)
        end if

        ! Fill ghost cells before refinement
        call a3_gc_tree(tree, i_eps, i_electron, a3_gc_interp, a3_bc_dirichlet_mirror)
        call a3_gc_tree(tree, i_eps, i_pos_ion, a3_gc_interp, a3_bc_dirichlet_mirror)

        call a3_adjust_refinement(tree, refine_routine, ref_info, &
             ST_refine_buffer_width, .true.)

        if (ref_info%n_add > 0 .or. ref_info%n_rm > 0) then
           ! For boxes which just have been refined, set data on their children
           call prolong_to_new_boxes(tree, ref_info)

           ! Compute the field on the new mesh
           call field_compute(tree, mg, .true.)
        end if

     end if
  end do

  call print_status()
  call a3_destroy(tree)

contains

  !> Initialize the AMR tree
  subroutine init_tree(tree)
    type(a3_t), intent(inout) :: tree

    ! Variables used below to initialize tree
    real(dp)                  :: dr
    integer                   :: id
    integer                   :: ix_list(3, 1) ! Spatial indices of initial boxes
    integer                   :: n_boxes_init = 1000

    dr = ST_domain_len / ST_box_size

    ! Initialize tree
    if (ST_cylindrical) then
       call a3_init(tree, ST_box_size, n_var_cell, n_var_face, dr, &
            coarsen_to=2, n_boxes=n_boxes_init, coord=af_cyl, &
            cc_names=ST_cc_names)
    else
       call a3_init(tree, ST_box_size, n_var_cell, n_var_face, dr, &
            coarsen_to=2, n_boxes=n_boxes_init, cc_names=ST_cc_names, fc_names=ST_fc_names)
    end if

    ! Set up geometry
    id             = 1          ! One box ...
    ix_list(:, id) = 1          ! With index 1,1 ...

    ! Create the base mesh
    call a3_set_base(tree, 1, ix_list)

  end subroutine init_tree

  ! This routine sets the cell refinement flags for box
  subroutine refine_routine(box, cell_flags)
    use m_geometry
    use m_init_cond_3d
    type(box3_t), intent(in) :: box
    ! Refinement flags for the cells of the box
    integer, intent(out)     :: cell_flags(box%n_cell, box%n_cell, box%n_cell)
    integer                  :: i, j, k, n, nc
    real(dp)                 :: cphi, dx, dx2
    real(dp)                 :: alpha, adx, fld, eps_max, eps_min
    real(dp)                 :: dist, rmin(3), rmax(3), grad_n(3), Efld(3)

    nc      = box%n_cell
    dx      = box%dr
    dx2     = box%dr**2
    eps_max = maxval(box%cc(:, :, :, i_eps))
    eps_min = minval(box%cc(:, :, :, i_eps))

    do k = 1, nc; do j = 1, nc; do i = 1, nc
       fld   = box%cc(i, j, k, i_electric_fld)
       alpha = LT_get_col(ST_td_tbl, i_alpha, ST_refine_adx_fac * fld) / ST_refine_adx_fac
       adx   = box%dr * alpha / ST_refine_adx
       cphi = dx2 * abs(box%cc(i, j, k, i_rhs)) / ST_refine_cphi
       
       

         grad_n(1) = 0.5_dp*(gradient_3d(box, 1, [i, j, k], 1, i_electron) + &
                     gradient_3d(box, 1, [i, j, k], 2, i_electron))
         grad_n(2) = 0.5_dp*(gradient_3d(box, 1, [i, j, k], 3, i_electron) + &
                     gradient_3d(box, 1, [i, j, k], 4, i_electron))
         grad_n(3) = 0.5_dp*(gradient_3d(box, 1, [i, j, k], 5, i_electron) + &
                     gradient_3d(box, 1, [i, j, k], 6, i_electron))

         Efld(1)   = 0.5_dp*(gradient_3d(box, 2, [i, j, k], 1, i_phi) + &
                     gradient_3d(box, 2, [i, j, k], 2, i_phi))
         Efld(2)   = 0.5_dp*(gradient_3d(box, 2, [i, j, k], 3, i_phi) + &
                     gradient_3d(box, 2, [i, j, k], 4, i_phi))
         Efld(2)   = 0.5_dp*(gradient_3d(box, 2, [i, j, k], 5, i_phi) + &
                     gradient_3d(box, 2, [i, j, k], 6, i_phi))

       
       if (adx + cphi > 1.0_dp .or. (eps_max > eps_min .and. box%lvl < 7)) then
            cell_flags(i, j, k) = af_do_ref
       else if (adx < 0.125_dp .and. cphi < 1.0_dp .and. dx < ST_derefine_dx  &
          .and. (eps_max == eps_min .or. box%lvl >= 7)) then
            cell_flags(i, j, k) = af_rm_ref
       else
          cell_flags(i, j, k) = af_keep_ref
       end if
       

       ! Refine around the initial conditions
       if (ST_time < ST_refine_init_time) then
          do n = 1, init_conds%n_cond
             dist = GM_dist_line(a3_r_cc(box, [i, j, k]), &
                  init_conds%seed_r0(:, n), &
                  init_conds%seed_r1(:, n), 3)
             if (dist - init_conds%seed_width(n) < 2 * dx &
                  .and. box%dr > ST_refine_init_fac * &
                  init_conds%seed_width(n)) then
                cell_flags(i, j, k) = af_do_ref
             end if
          end do
       end if

    end do; end do; end do

    ! Check fixed refinements
    rmin = box%r_min
    rmax = box%r_min + box%dr * box%n_cell

    do n = 1, size(ST_refine_regions_dr)
       if (ST_time <= ST_refine_regions_tstop(n) .and. &
            dx > ST_refine_regions_dr(n) .and. all(&
            rmax >= ST_refine_regions_rmin(:, n) .and. &
            rmin <= ST_refine_regions_rmax(:, n))) then
          ! Mark just the center cell to prevent refining neighbors
          cell_flags(nc/2, nc/2, nc/2) = af_do_ref
       end if
    end do
    

    ! Make sure we don't have or get a too fine or too coarse a grid
    if (dx > ST_refine_max_dx) then
       cell_flags = af_do_ref
    else if (dx < 2 * ST_refine_min_dx) then
       where(cell_flags == af_do_ref) cell_flags = af_keep_ref
    end if
    

  end subroutine refine_routine

  !> Compute the electron fluxes due to drift and diffusion
  subroutine fluxes_koren(boxes, id)
    type(box3_t), intent(inout) :: boxes(:)
    integer, intent(in)          :: id
    real(dp)                     :: inv_dr, fld
    ! Velocities at cell faces
    real(dp), allocatable        :: v(:, :, :, :)
    ! Diffusion coefficients at cell faces
    real(dp), allocatable        :: dc(:, :, :, :)
    ! Cell-centered densities
    real(dp), allocatable        :: cc(:, :, :)
    real(dp), allocatable        :: eps(:, :, :)
    real(dp), allocatable        :: sigma(:, :, :, :)
    integer                      :: nc, i, j, k

    nc     = boxes(id)%n_cell
    inv_dr = 1.0_dp/boxes(id)%dr
    
    allocate(v(1:nc+1, 1:nc+1, 1:nc+1, 3))
    allocate(dc(1:nc+1, 1:nc+1, 1:nc+1, 3))
    allocate(cc(-1:nc+2, -1:nc+2, -1:nc+2))
    allocate(eps(-1:nc+2, -1:nc+2, -1:nc+2))
    allocate(sigma(0:nc+2, 0:nc+2, 0:nc+2, 3))

    ! Fill ghost cells on the sides of boxes (no corners)
    call a3_gc_box(boxes, id, i_eps, i_electron, a3_gc_interp, &
         a3_bc_dirichlet_mirror, .false.)
         
    call a3_gc_box_fc(boxes, id, i_eps, sigma_rhs, a3_gc_prolong_copy_fc, &
         a3_bc_dirichlet_zero_fc, med)
         
    call a3_gc2_box(boxes, id, i_eps, i_electron, a3_gc2_prolong_linear, &
         a3_bc2_dirichlet_mirror, cc, nc, med)

    

    ! We use the average field to compute the mobility and diffusion coefficient
    ! at the interface
    
    do k = 1, nc+1; do j = 1, nc+1; do i = 1, nc+1  

       fld = abs(boxes(id)%fc(i, j, k, 1, electric_fld))
       v(i, j, k, 1)  = -LT_get_col(ST_td_tbl, i_mobility, fld) * &
            boxes(id)%fc(i, j, k, 1, electric_fld)
       if (boxes(id)%cc(i, j, k, i_eps) > med .and. boxes(id)%cc(i-1, j, k, i_eps) <= med) then
         v(i, j, k, 1) = -LT_get_col(ST_td_tbl, i_mobility, fld) * &
              (boxes(id)%fc(i, j, k, 1, electric_fld)*boxes(id)%cc(i, j, k, i_eps) - &
              boxes(id)%fc(i, j, k, 1, sigma_rhs))/boxes(id)%cc(i-1, j, k, i_eps)
       end if
       dc(i, j, k, 1) = LT_get_col(ST_td_tbl, i_diffusion, fld)

       fld = abs(boxes(id)%fc(i, j, k, 2, electric_fld))
       v(i, j, k, 2)  = -LT_get_col(ST_td_tbl, i_mobility, fld) * &
            boxes(id)%fc(i, j, k, 2, electric_fld)
       if (boxes(id)%cc(i, j, k, i_eps) > med .and. boxes(id)%cc(i, j-1, k, i_eps) <= med) then
         v(i, j, k, 2) = -LT_get_col(ST_td_tbl, i_mobility, fld) * &
              (boxes(id)%fc(i, j, k, 2, electric_fld)*boxes(id)%cc(i, j, k, i_eps) - &
              boxes(id)%fc(i, j, k, 2, sigma_rhs))/boxes(id)%cc(i, j-1, k, i_eps)
        end if
       dc(i, j, k, 2) = LT_get_col(ST_td_tbl, i_diffusion, fld)

       fld = abs(boxes(id)%fc(i, j, k, 3, electric_fld))
       v(i, j, k, 3)  = -LT_get_col(ST_td_tbl, i_mobility, fld) * &
            boxes(id)%fc(i, j, k, 3, electric_fld)
       if (boxes(id)%cc(i, j, k, i_eps) > med .and. boxes(id)%cc(i-1, j, k, i_eps) <= med) then
         v(i, j, k, 1) = -LT_get_col(ST_td_tbl, i_mobility, fld) * &
              (boxes(id)%fc(i, j, k, 1, electric_fld)*boxes(id)%cc(i, j, k, i_eps) - &
              boxes(id)%fc(i, j, k, 1, sigma_rhs))/boxes(id)%cc(i-1, j, k, i_eps)
       end if
       dc(i, j, k, 3) = LT_get_col(ST_td_tbl, i_diffusion, fld)

    end do; end do; end do

    call flux_koren_3d(cc, v, nc, inv_dr)
    call flux_diff_3d(boxes(id), dc, i_electron)

    boxes(id)%fc(:, :, :, :, flux_elec) = v + dc
    

    if (ST_update_ions) then
       ! Use a constant diffusion coefficient for ions
       dc = ST_ion_diffusion

       ! Use a constant mobility for ions
       v(:, :, :, 1:3) = ST_ion_mobility * &
            boxes(id)%fc(:, :, :, 1:3, electric_fld)

       ! Fill ghost cells on the sides of boxes (no corners)
       call a3_gc_box(boxes, id, i_eps, i_pos_ion, a3_gc_interp, &
            a3_bc_dirichlet_mirror, .false.)

       call flux_koren_3d(cc, v, nc, inv_dr)
       call flux_diff_3d(boxes(id), dc, i_pos_ion)

       boxes(id)%fc(:, :, :, :, flux_ion) = v + dc
    end if
  end subroutine fluxes_koren

  !> Take average of new and old electron/ion density for explicit trapezoidal rule
  subroutine average_density(box)
    type(box3_t), intent(inout) :: box
    box%cc(:, :, :, i_electron) = 0.5_dp *  &
         (box%cc(:, :, :, i_electron) + box%cc(:, :, :, i_electron_old))
    box%cc(:, :, :, i_pos_ion) = 0.5_dp * &
         (box%cc(:, :, :, i_pos_ion)  + box%cc(:, :, :, i_pos_ion_old))
  end subroutine average_density

  !> Advance solution over dt based on the fluxes / source term, using forward Euler
  subroutine update_solution(box, args)
    use omp_lib
    use m_units_constants
    type(box3_t), intent(inout) :: box
    real(dp), intent(in)         :: args(:)
    real(dp)                     :: dt, inv_dr, src, fld, fld_vec(3), s_flow(3)
    real(dp)                     :: alpha, eta, f_elec, f_ion, mu, diff
    real(dp)                     :: dt_cfl, dt_dif, dt_drt, s_fac

    integer                      :: i, j, k, nc, tid, i_step
    type(LT_loc_t)               :: loc

    tid      = omp_get_thread_num() + 1
    dt       = args(1)
    i_step   = nint(args(2))

    nc     = box%n_cell
    inv_dr = 1/box%dr
    f_ion  = 0.0_dp
    

    do k = 1, nc; do j = 1, nc; do i = 1, nc
       fld   = box%cc(i, j, k, i_electric_fld)
       loc   = LT_get_loc(ST_td_tbl, fld)
       alpha = LT_get_col_at_loc(ST_td_tbl, i_alpha, loc)*a3_heaviside(box%cc(i, j, k, i_eps), med)
       eta   = LT_get_col_at_loc(ST_td_tbl, i_eta, loc)*a3_heaviside(box%cc(i, j, k, i_eps), med)
       mu    = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)
       f_elec = 0.0_dp

       ! Contribution of flux

       if (box%cc(i, j, k, i_eps) <= med) then
         f_elec = (sum(box%fc(i, j, k, 1:3, flux_elec)) - &
              box%fc(i+1, j, k, 1, flux_elec) - &
              box%fc(i, j+1, k, 2, flux_elec) - &
              box%fc(i, j, k+1, 3, flux_elec)) * inv_dr * dt
       else
         f_elec = 0.0_dp
       end if

       if (ST_update_ions) then

          f_ion = (sum(box%fc(i, j, k, 1:3, flux_ion)) - &
               box%fc(i+1, j, k, 1, flux_ion) - &
               box%fc(i, j+1, k, 2, flux_ion) - &
               box%fc(i, j, k+1, 3, flux_ion)) * inv_dr * dt

       end if

       ! Source term
       src = fld * mu * box%cc(i, j, k, i_electron) * (alpha - eta)
       if (photoi_enabled .and. box%cc(i, j, k, i_eps) <= med) src = src + box%cc(i, j, k, i_photo)

       if (i_step == 1 .and. ST_output_src_term) then
          if (ST_output_src_decay_rate < 0) then
             ! No time-averaging
             box%cc(i, j, k, i_src) = src
          else
             ! Approximate exp(-t*rate) as 1-t*rate
             box%cc(i, j, k, i_src) = (1 - dt * ST_output_src_decay_rate) * &
                  box%cc(i, j, k, i_src) + src * dt
          end if
       end if

       ! Convert to density
       src  = src * dt

       ! Add flux and source term
       if (box%cc(i, j, k, i_eps) <= med) then
         box%cc(i, j, k, i_electron) = box%cc(i, j, k, i_electron) + f_elec + src
       
         ! Add flux and source term
         box%cc(i, j, k, i_pos_ion)  = box%cc(i, j, k, i_pos_ion) + f_ion + src
       end if

       ! Possible determine new time step restriction
       if (i_step == 2) then

          fld_vec(1) = 0.5_dp * (box%fc(i, j, k, 1, electric_fld) + &
               box%fc(i+1, j, k, 1, electric_fld))
          fld_vec(2) = 0.5_dp * (box%fc(i, j, k, 2, electric_fld) + &
               box%fc(i, j+1, k, 2, electric_fld))
          fld_vec(3) = 0.5_dp * (box%fc(i, j, k, 3, electric_fld) + &
               box%fc(i, j, k+1, 3, electric_fld))

          diff = LT_get_col(ST_td_tbl, i_diffusion, fld)

          ! CFL condition
          dt_cfl = 1.0_dp/sum(abs(fld_vec * max(mu, abs(ST_ion_mobility), epsilon(1.0_dp))) * inv_dr)
          ! Diffusion condition
          dt_dif = box%dr**2 / max(2 * 3 * max(diff, ST_ion_diffusion), epsilon(1.0_dp))
          ! Dielectric relaxation time
          dt_drt = UC_eps0 / (UC_elem_charge * (mu + abs(ST_ion_mobility)) * &
               max(box%cc(i, j, k, i_electron), epsilon(1.0_dp)))

          ! Take the minimum of the CFL condition with Courant number 0.5 and
          ! the combined CFL-diffusion condition with Courant number 1.0. The
          ! 0.5 is emperical, to have good accuracy (and TVD/positivity) in
          ! combination with the explicit trapezoidal rule
          dt_cfl = min(0.5_dp * dt_cfl, 1/(1/dt_cfl + 1/dt_dif))

          ST_dt_matrix(ST_ix_cfl, tid) = min(ST_dt_matrix(ST_ix_cfl, tid), &
               dt_cfl)
          ST_dt_matrix(ST_ix_drt, tid) = min(ST_dt_matrix(ST_ix_drt, tid), &
               dt_drt)
          ST_dt_matrix(ST_ix_diff, tid) = min(ST_dt_matrix(ST_ix_diff, tid), &
               dt_dif)
       end if

    end do; end do; end do
  end subroutine update_solution

  !> For each box that gets refined, set data on its children using this routine
  subroutine prolong_to_new_boxes(tree, ref_info)
    use m_a3_prolong
    type(a3_t), intent(inout)    :: tree
    type(ref_info_t), intent(in)  :: ref_info
    integer                       :: lvl, i, id, p_id

    !$omp parallel private(lvl, i, id, p_id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          p_id = tree%boxes(id)%parent
          
          call define_DI(tree%boxes(id)) ! For dielectric re-shaping
          call a3_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_electron, i_eps = i_eps, med = med)
          call a3_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_pos_ion, i_eps = i_eps, med = med)
          call a3_prolong_sparse(tree%boxes(p_id), tree%boxes(id), i_phi, i_eps = i_eps)
          if (photoi_enabled) then
             call a3_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_photo, i_eps = i_eps, med = med)
          end if
          if (ST_output_src_term) then
             call a3_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_src, i_eps = i_eps, med = med)
          end if
          call define_DI(tree%boxes(id))
       end do
       !$omp end do

       !$omp do
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a3_gc_box(tree%boxes, id, i_eps, i_electron, &
               a3_gc_interp, a3_bc_dirichlet_mirror)
          call a3_gc_box(tree%boxes, id, i_eps, i_pos_ion, &
               a3_gc_interp, a3_bc_dirichlet_mirror)
          call a3_gc_box(tree%boxes, id, i_eps, i_phi, &
               mg%sides_rb, mg%sides_bc)
          if (photoi_enabled) then
             call a3_gc_box(tree%boxes, id, i_eps, i_photo, &
               a3_gc_interp, photoi_helmh_bc)
          end if
          if (ST_output_src_term) then
             call a3_gc_box(tree%boxes, id, i_eps, i_src, &
               a3_gc_interp, a3_bc_neumann_zero)
          end if
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine prolong_to_new_boxes

  subroutine write_log_file(tree, filename, out_cnt, dir)
    type(a3_t), intent(in)      :: tree
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: out_cnt
    character(len=*), intent(in) :: dir
    character(len=ST_slen)       :: fname
    character(len=20), save      :: fmt
    integer, parameter           :: my_unit = 123
    real(dp)                     :: velocity, dt
    real(dp), save               :: prev_pos(3) = 0
    real(dp)                     :: sum_elec, sum_pos_ion
    real(dp)                     :: max_elec, max_field, max_Er
    type(a3_loc_t)              :: loc_elec, loc_field, loc_Er

    call a3_prepend_directory(filename, dir, fname)

    call a3_tree_sum_cc(tree, i_electron, sum_elec)
    call a3_tree_sum_cc(tree, i_pos_ion, sum_pos_ion)
    call a3_tree_max_cc(tree, i_electron, max_elec, loc_elec)
    call a3_tree_max_cc(tree, i_electric_fld, max_field, loc_field)
    call a3_tree_max_fc(tree, 1, electric_fld, max_Er, loc_Er)

    dt = ST_dt_safety_factor * minval(ST_dt_matrix)

    if (out_cnt == 1) then
       open(my_unit, file=trim(fname), action="write")

       write(my_unit, *) "# it time dt v sum(n_e) sum(n_i) ", &
            "max(E) x y z max(n_e) x y z wc_time n_cells"
       fmt = "(I6,14E16.8,I12)"

       close(my_unit)

       ! Start with velocity zero
       prev_pos = a3_r_loc(tree, loc_field)
    end if

    velocity = norm2(a3_r_loc(tree, loc_field) - prev_pos) / ST_dt_output
    prev_pos = a3_r_loc(tree, loc_field)

    open(my_unit, file=trim(fname), action="write", &
         position="append")

    write(my_unit, fmt) out_cnt, ST_time, dt, velocity, sum_elec, &
         sum_pos_ion, max_field, a3_r_loc(tree, loc_field), max_elec, &
         a3_r_loc(tree, loc_elec), wc_time, a3_num_cells_used(tree)

    close(my_unit)

  end subroutine write_log_file
  
  subroutine print_status()
    write(*, "(F7.3,A,I0,A,E10.3,A,E10.3,A,E10.3)") &
             100 * ST_time / ST_end_time, "% it: ", it, &
             " t:", ST_time, " dt:", ST_dt, " wc:", wc_time
  end subroutine print_status
  

end program streamer_3d

