#include "../afivo/src/cpp_macros_$Dd.h"
!> Program to perform $Dd streamer simulations in Cartesian and cylindrical coordinates
program streamer_$Dd

  use m_a$D_all
  use m_streamer
  use m_field_$Dd
  use m_init_cond_$Dd
  use m_photoi_$Dd
  use m_flux_schemes_$Dd
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
  type(a$D_t)            :: tree      ! This contains the full grid information
  type(mg$D_t)           :: mg        ! Multigrid option struct

  type(ref_info_t)       :: ref_info

  integer :: output_cnt = 0 ! Number of output files written

  call CFG_update_from_arguments(cfg)
  call ST_initialize(cfg, $D)
  call photoi_initialize(cfg)

  call ST_load_transport_data(cfg)
  call field_initialize(cfg, mg)
  call init_cond_initialize(cfg, $D)

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
  mg%box_op => mg$D_auto_op
  mg%box_gsrb => mg$D_auto_gsrb
  mg%box_corr => mg$D_auto_corr

  ! This routine always needs to be called when using multigrid
  call mg$D_init_mg(mg)

  output_cnt      = 0         ! Number of output files written
  ST_time         = 0         ! Simulation time (all times are in s)
  med = a$D_harm_w(1.0_dp, ST_epsilon_die, 0.5_dp)
  

  ! Set up the initial conditions
  do
     call a$D_loop_box(tree, init_cond_set_box)
     call field_compute(tree, mg, .false.)
     call a$D_adjust_refinement(tree, refine_routine, ref_info, &
          ST_refine_buffer_width, .true.)
     if (ref_info%n_add == 0) exit
  end do

  call init_cond_stochastic_density(tree)

  print *, "Number of threads", af_get_max_threads()
  call a$D_print_info(tree)


  ! Start from small time step
  ST_dt   = ST_dt_min
  dt_prev = ST_dt
  
  ! Initial wall clock time
  call system_clock(t_start, count_rate)
  inv_count_rate = 1.0_dp / count_rate
  time_last_print = -1e10_dp

  

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
     call a$D_tree_copy_cc(tree, i_electron, i_electron_old)
     call a$D_tree_copy_cc(tree, i_pos_ion, i_pos_ion_old)
     

     ST_dt_matrix = ST_dt_max      ! Maximum time step
     ! Two forward Euler steps over ST_dt
     do i = 1, 2
        ST_time = ST_time + ST_dt
        

        ! First calculate fluxes
        call a$D_loop_boxes(tree, fluxes_koren, .true.)
        call a$D_consistent_fluxes(tree, [flux_elec])
        
        ! Update the solution
        call a$D_loop_box_arg(tree, update_solution, [ST_dt, real(i, dp)], .true.)
        
        if (i == 1) then
           ! Compute new field on first iteration
           call field_compute(tree, mg, .true.)

           ! Coarse densities might be required for ghost cells
           call a$D_restrict_tree(tree, i_electron, i_eps = i_eps, med = med)
           call a$D_restrict_tree(tree, sigma_rhs, i_eps = i_eps, med = med, s_flag = .true.)
        end if
     end do
     

     ST_time = ST_time - ST_dt        ! Go back one time step
     

     ! Take average of n and sigma (explicit trapezoidal rule)
     call a$D_loop_box(tree, average_density)

     ! Coarse densities might be required for ghost cells
     call a$D_restrict_tree(tree, i_electron, i_eps = i_eps, med = med)
     call a$D_restrict_tree(tree, sigma_rhs, i_eps = i_eps, med = med, s_flag = .true.)

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
        call a$D_gc_tree(tree, i_eps, i_electron, a$D_gc_interp_lim, a$D_bc_neumann_zero, med = med)
        call a$D_gc_tree(tree, i_eps, i_pos_ion, a$D_gc_interp_lim, a$D_bc_neumann_zero, med = med)

        if (ST_output_src_term) then
           call a$D_restrict_tree(tree, i_src, i_eps = i_eps, med = med)
           call a$D_gc_tree(tree, i_eps, i_src, a$D_gc_interp_lim, a$D_bc_neumann_zero, med = med)
        end if

        if (photoi_enabled) then
           call photoi_set_src(tree, ST_dt) ! Because the mesh could have changed
           call a$D_restrict_tree(tree, i_photo, i_eps = i_eps, med = med)
           call a$D_gc_tree(tree, i_eps, i_photo, a$D_gc_interp_lim, a$D_bc_neumann_zero, med = med)
        end if

        write(fname, "(A,I6.6)") trim(ST_simulation_name) // "_", output_cnt
        call a$D_write_silo(tree, fname, output_cnt, ST_time, &
             vars_for_output, dir=ST_output_dir)

        if (ST_datfile_write) then
           call a$D_write_tree(tree, fname, ST_output_dir)
        end if

        write(fname, "(A,I6.6)") trim(ST_simulation_name) // "_log.txt"
        call write_log_file(tree, fname, output_cnt, ST_output_dir)

        if (ST_plane_write) then
           write(fname, "(A,I6.6)") trim(ST_simulation_name) // &
                "_plane_", output_cnt
           call a$D_write_plane(tree, fname, [ST_plane_ivar], &
                ST_plane_rmin * ST_domain_len, &
                ST_plane_rmax * ST_domain_len, &
                ST_plane_npixels, ST_output_dir)
        end if

        if (ST_lineout_write) then
           write(fname, "(A,I6.6)") trim(ST_simulation_name) // &
                "_line_", output_cnt
           call a$D_write_line(tree, trim(fname), &
                [i_electron, i_pos_ion, i_phi, i_electric_fld, i_eps], &
                r_min = ST_lineout_rmin(1:$D) * ST_domain_len, &
                r_max = ST_lineout_rmax(1:$D) * ST_domain_len, &
                n_points=ST_lineout_npoints, dir=ST_output_dir)
        end if
     end if
     


     if (mod(it, ST_refine_per_steps) == 0) then
        ! Restrict electron and ion densities before refinement
        call a$D_restrict_tree(tree, i_electron, i_eps = i_eps, med = med)
        call a$D_restrict_tree(tree, i_pos_ion, i_eps = i_eps, med = med)
        call a$D_restrict_tree(tree, sigma_rhs, i_eps = i_eps, med = med, s_flag = .true.)

        if (ST_output_src_term .and. ST_output_src_decay_rate > 0) then
           ! Have to set src term on coarse grids as well, and fill ghost cells
           ! before prolongation
           call a$D_restrict_tree(tree, i_src, i_eps = i_eps, med = med)
           call a$D_gc_tree(tree, i_eps, i_src, a$D_gc_interp_lim, a$D_bc_neumann_zero, med = med)
        end if

        ! Fill ghost cells before refinement
        call a$D_gc_tree(tree, i_eps, i_electron, a$D_gc_interp_lim, a$D_bc_neumann_zero, med = med)
        call a$D_gc_tree(tree, i_eps, i_pos_ion, a$D_gc_interp_lim, a$D_bc_neumann_zero, med = med)

        call a$D_adjust_refinement(tree, refine_routine, ref_info, &
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
  call a$D_destroy(tree)

contains

  !> Initialize the AMR tree
  subroutine init_tree(tree)
    type(a$D_t), intent(inout) :: tree

    ! Variables used below to initialize tree
    real(dp)                  :: dr
    integer                   :: id
    integer                   :: ix_list($D, 2) ! Spatial indices of initial boxes
    integer                   :: n_boxes_init = 1000

    dr = ST_domain_len / ST_box_size

    ! Initialize tree
    if (ST_cylindrical) then
       call a$D_init(tree, ST_box_size, n_var_cell, n_var_face, dr, &
            coarsen_to=2, n_boxes=n_boxes_init, coord=af_cyl, &
            cc_names=ST_cc_names)
    else
       call a$D_init(tree, ST_box_size, n_var_cell, n_var_face, dr, &
            coarsen_to=2, n_boxes=n_boxes_init, cc_names=ST_cc_names, fc_names=ST_fc_names)
    end if

    ! Set up geometry
    id             = 2          
#if $D == 2
    ix_list(:, 1) = [1, 1]     
    ix_list(:, 2) = [1, 2]
#elif $D == 3
    ix_list(:, 1) = [1, 1, 1]     
    ix_list(:, 2) = [1, 1, 2]
#endif

    ! Create the base mesh
    call a$D_set_base(tree, 2, ix_list)

  end subroutine init_tree

  ! This routine sets the cell refinement flags for box
  subroutine refine_routine(box, cell_flags)
    use m_geometry
    use m_init_cond_$Dd
    type(box$D_t), intent(in) :: box
    ! Refinement flags for the cells of the box
    integer, intent(out)     :: cell_flags(DTIMES(box%n_cell))
    integer                  :: IJK, n, nc
    real(dp)                 :: cphi, dx, dx2
    real(dp)                 :: alpha, adx, fld, eps_max, eps_min
    real(dp)                 :: dist, rmin($D), rmax($D), grad_n($D), Efld($D)

    nc      = box%n_cell
    dx      = box%dr
    dx2     = box%dr**2
    eps_max = maxval(box%cc(DTIMES(:), i_eps))
    eps_min = minval(box%cc(DTIMES(:), i_eps))

    do KJI_DO(1,nc)
       fld   = box%cc(IJK, i_electric_fld)
       alpha = LT_get_col(ST_td_tbl, i_alpha, ST_refine_adx_fac * fld) / ST_refine_adx_fac
       adx   = box%dr * alpha / ST_refine_adx
       cphi = dx2 * abs(box%cc(IJK, i_rhs)) / ST_refine_cphi
       
              
       
       if (adx + cphi > 1.0_dp .or. (eps_max > eps_min .and. box%lvl < 6)) then
            cell_flags(IJK) = af_do_ref
       else if (adx < 0.125_dp .and. cphi < 1.0_dp .and. dx < ST_derefine_dx  &
          .and. (eps_max == eps_min .or. box%lvl >= 6)) then
            cell_flags(IJK) = af_rm_ref
       else
          cell_flags(IJK) = af_keep_ref
       end if
       

       ! Refine around the initial conditions
       if (ST_time < ST_refine_init_time) then
          do n = 1, init_conds%n_cond
             dist = GM_dist_line(a$D_r_cc(box, [IJK]), &
                  init_conds%seed_r0(:, n), &
                  init_conds%seed_r1(:, n), $D)
             if (dist - init_conds%seed_width(n) < 2 * dx &
                  .and. box%dr > ST_refine_init_fac * &
                  init_conds%seed_width(n)) then
                cell_flags(IJK) = af_do_ref
             end if
          end do
       end if

    end do; CLOSE_DO

    ! Check fixed refinements
    rmin = box%r_min
    rmax = box%r_min + box%dr * box%n_cell

    do n = 1, size(ST_refine_regions_dr)
       if (ST_time <= ST_refine_regions_tstop(n) .and. &
            dx > ST_refine_regions_dr(n) .and. all(&
            rmax >= ST_refine_regions_rmin(:, n) .and. &
            rmin <= ST_refine_regions_rmax(:, n))) then
          ! Mark just the center cell to prevent refining neighbors
          cell_flags(DTIMES(nc/2)) = af_do_ref
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
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)          :: id
    real(dp)                     :: inv_dr, fld
    ! Velocities at cell faces
    real(dp), allocatable        :: v(DTIMES(:), :)
    ! Diffusion coefficients at cell faces
    real(dp), allocatable        :: dc(DTIMES(:), :)
    ! Cell-centered densities
    real(dp), allocatable        :: cc(DTIMES(:))
    real(dp), allocatable        :: eps(DTIMES(:))
    real(dp), allocatable        :: sigma(DTIMES(:), :)
    integer                      :: nc, IJK

    nc     = boxes(id)%n_cell
    inv_dr = 1.0_dp/boxes(id)%dr
    
    allocate(v(DTIMES(1:nc+1), $D))
    allocate(dc(DTIMES(1:nc+1), $D))
    allocate(cc(DTIMES(-1:nc+2)))
    allocate(eps(DTIMES(-1:nc+2)))
    allocate(sigma(DTIMES(0:nc+2), $D))

    ! Fill ghost cells on the sides of boxes (no corners)
    call a$D_gc_box(boxes, id, i_eps, i_electron, a$D_gc_interp_lim, &
         a$D_bc_neumann_zero, .false., med = med)
         
    call a$D_gc_box_fc(boxes, id, i_eps, sigma_rhs, a$D_gc_prolong_copy_fc, &
         a$D_bc_dirichlet_zero_fc, med)
         
    call a$D_gc2_box(boxes, id, i_eps, i_electron, a$D_gc2_prolong_linear, &
         a$D_bc2_neumann_zero, cc, nc, med)

    

    ! We use the average field to compute the mobility and diffusion coefficient
    ! at the interface
    
    do KJI_DO(1,nc+1)  
#if $D == 2      
      fld = abs(boxes(id)%fc(IJK, 1, electric_fld))
      v(IJK, 1) = -LT_get_col(ST_td_tbl, i_mobility, fld) * &
           boxes(id)%fc(IJK, 1, electric_fld)
      dc(IJK, 1) = LT_get_col(ST_td_tbl, i_diffusion, fld)
      
      fld = abs(boxes(id)%fc(IJK, 2, electric_fld))
      v(IJK, 2) = -LT_get_col(ST_td_tbl, i_mobility, fld) * &
           boxes(id)%fc(IJK, 2, electric_fld)
      dc(IJK, 2) = LT_get_col(ST_td_tbl, i_diffusion, fld)

#elif $D == 3
       fld = abs(boxes(id)%fc(IJK, 1, electric_fld))
       v(IJK, 1)  = -LT_get_col(ST_td_tbl, i_mobility, fld) * &
            boxes(id)%fc(IJK, 1, electric_fld)
       dc(IJK, 1) = LT_get_col(ST_td_tbl, i_diffusion, fld)

       fld = abs(boxes(id)%fc(IJK, 2, electric_fld))
       v(IJK, 2)  = -LT_get_col(ST_td_tbl, i_mobility, fld) * &
            boxes(id)%fc(IJK, 2, electric_fld)
       dc(IJK, 2) = LT_get_col(ST_td_tbl, i_diffusion, fld)

       fld = abs(boxes(id)%fc(IJK, 3, electric_fld))
       v(IJK, 3)  = -LT_get_col(ST_td_tbl, i_mobility, fld) * &
            boxes(id)%fc(IJK, 3, electric_fld)
       dc(IJK, 3) = LT_get_col(ST_td_tbl, i_diffusion, fld)
#endif
    end do; CLOSE_DO

    call flux_koren_$Dd(cc, v, nc, inv_dr)
    call flux_diff_$Dd(boxes(id), dc, i_electron)

    boxes(id)%fc(DTIMES(:), :, flux_elec) = v + dc
    

    if (ST_update_ions) then
       ! Use a constant diffusion coefficient for ions
       dc = ST_ion_diffusion

       ! Use a constant mobility for ions
       v(DTIMES(:), 1:$D) = ST_ion_mobility * &
            boxes(id)%fc(DTIMES(:), 1:$D, electric_fld)

       ! Fill ghost cells on the sides of boxes (no corners)
       call a$D_gc_box(boxes, id, i_eps, i_pos_ion, a$D_gc_interp_lim, &
            a$D_bc_neumann_zero, .false., med = med)

       call flux_koren_$Dd(cc, v, nc, inv_dr)
       call flux_diff_$Dd(boxes(id), dc, i_pos_ion)

       boxes(id)%fc(DTIMES(:), :, flux_ion) = v + dc
    end if
  end subroutine fluxes_koren

  !> Take average of new and old electron/ion density for explicit trapezoidal rule
  subroutine average_density(box)
    type(box$D_t), intent(inout) :: box
    box%cc(DTIMES(:), i_electron) = 0.5_dp *  &
         (box%cc(DTIMES(:), i_electron) + box%cc(DTIMES(:), i_electron_old))
    box%cc(DTIMES(:), i_pos_ion) = 0.5_dp * &
         (box%cc(DTIMES(:), i_pos_ion)  + box%cc(DTIMES(:), i_pos_ion_old))
  end subroutine average_density

  !> Advance solution over dt based on the fluxes / source term, using forward Euler
  subroutine update_solution(box, args)
    use omp_lib
    use m_units_constants
    type(box$D_t), intent(inout) :: box
    real(dp), intent(in)         :: args(:)
    real(dp)                     :: dt, inv_dr, src, fld, fld_vec($D), s_flow($D)
    real(dp)                     :: alpha, eta, f_elec, f_ion, mu, diff, rand
    real(dp)                     :: dt_cfl, dt_dif, dt_drt, s_fac
#if $D == 2
    real(dp)                     :: rfac(2)
    integer                      :: ioff
#endif
    integer                      :: IJK, nc, tid, i_step, n, flag, ix($D), count_rate, time, dir
    type(LT_loc_t)               :: loc

    tid      = omp_get_thread_num() + 1
    dt       = args(1)
    i_step   = nint(args(2))

    nc     = box%n_cell
    inv_dr = 1/box%dr
    f_ion  = 0.0_dp
    
#if $D == 2
    ioff   = (box%ix(1)-1) * nc
#endif

    do KJI_DO(1,nc)
       fld   = box%cc(IJK, i_electric_fld)
       loc   = LT_get_loc(ST_td_tbl, fld)
       alpha = LT_get_col_at_loc(ST_td_tbl, i_alpha, loc)
       eta   = LT_get_col_at_loc(ST_td_tbl, i_eta, loc)
       mu    = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)
       f_elec = 0.0_dp
       flag = 0

       ! Contribution of flux
#if $D == 2
       if (ST_cylindrical) then
          ! Weighting of flux contribution for cylindrical coordinates
          rfac(:) = [i+ioff-1, i+ioff] / (i+ioff-0.5_dp)
       else
          rfac(:) = 1.0_dp
       end if
       
       if (box%cc(IJK, i_eps) < ST_epsilon_die) then
         
         if (box%cc(i-1, j, i_eps) == ST_epsilon_die) then
           s_flow(1) = - dt  * box%fc(IJK, 1, flux_elec) 
           box%fc(IJK, 1, sigma_rhs) = box%fc(IJK, 1, sigma_rhs) + fac * max(s_flow(1), 0.0_dp)
           f_elec = f_elec - (s_flow(1) + box%fc(i+1, j, 1, flux_elec) * dt ) * inv_dr
         else
           f_elec = f_elec + (box%fc(IJK, 1, flux_elec) - box%fc(i+1, j, 1, flux_elec)) * inv_dr * dt          
         end if
         if (box%cc(i, j-1, i_eps) == ST_epsilon_die) then
           s_flow(2) = - dt  * box%fc(IJK, 2, flux_elec) 
           box%fc(IJK, 2, sigma_rhs) = box%fc(IJK, 2, sigma_rhs) + fac * max(s_flow(2), 0.0_dp)
           f_elec = f_elec - (s_flow(2) + box%fc(i, j+1, 2, flux_elec) * dt) * inv_dr
         else
           f_elec = f_elec + (box%fc(IJK, 2, flux_elec) - box%fc(i, j+1, 2, flux_elec)) * inv_dr * dt          
         end if
         
       else 
         
         if (box%cc(i-1, j, i_eps) < ST_epsilon_die) then
           s_flow(1) = dt  * box%fc(IJK, 1, flux_elec) 
           box%fc(IJK, 1, sigma_rhs) = box%fc(IJK, 1, sigma_rhs) + fac * max(s_flow(1), 0.0_dp)         
         end if
         if (box%cc(i, j-1, i_eps) < ST_epsilon_die) then
           s_flow(2) = dt  * box%fc(IJK, 2, flux_elec) 
           box%fc(IJK, 2, sigma_rhs) = box%fc(IJK, 2, sigma_rhs) + fac * max(s_flow(2), 0.0_dp)
         end if
         
       end if
              
#elif $D == 3

     if (box%cc(IJK, i_eps) <= med) then
       if (box%cc(i-1, j, k, i_eps) > med) then
         if (i > 1) then
           s_flow(1) = dt * (box%fc(i-1, j, k, 1, flux_elec) - &
           box%fc(IJK, 1, flux_elec)) 
         else
           s_flow(1) = 0.0_dp
         end if
         box%fc(IJK, 1, sigma_rhs) = box%fc(IJK, 1, sigma_rhs) + fac * max(s_flow(1), 0.0_dp) 
  end if
  if (box%cc(i, j-1, k, i_eps) > med) then
    if (j > 1) then
      s_flow(2) = dt * (box%fc(i, j-1, k, 2, flux_elec) - &
          box%fc(IJK, 2, flux_elec))
    else
      s_flow(2) = 0.0_dp
    end if 
    box%fc(IJK, 2, sigma_rhs) = box%fc(IJK, 2, sigma_rhs) + fac * max(s_flow(2), 0.0_dp) 
  end if  
  if (box%cc(i, j, k-1, i_eps) > med) then
    if (k > 1) then
      s_flow(3) = dt * (box%fc(i, j, k-1, 3, flux_elec) - &
          box%fc(IJK, 3, flux_elec))
    else
      s_flow(3) = 0.0_dp
    end if 
    box%fc(IJK, 3, sigma_rhs) = box%fc(IJK, 3, sigma_rhs) + fac * max(s_flow(3), 0.0_dp) 
  end if  
  
  f_elec    = (sum(box%fc(IJK, :, flux_elec)) - box%fc(i, j+1, k, 2, flux_elec) - &
        box%fc(i+1, j, k, 1, flux_elec) - box%fc(i, j, k+1, 3, flux_elec) ) * inv_dr * dt 
else 
  f_elec    = 0.0_dp
  if (box%cc(i-1, j, k, i_eps) <= med) then
    s_flow(1) = dt * (box%fc(IJK, 1, flux_elec) - &
        box%fc(i+1, j, k, 1, flux_elec)) 
    box%fc(IJK, 1, sigma_rhs) = box%fc(IJK, 1, sigma_rhs) + fac * max(s_flow(1), 0.0_dp) 
  end if
  if (box%cc(i, j-1, k, i_eps) <= med) then
    s_flow(2) = dt * (box%fc(IJK, 2, flux_elec) - &
        box%fc(i, j+1, k, 2, flux_elec)) 
    box%fc(IJK, 2, sigma_rhs) = box%fc(IJK, 2, sigma_rhs) + fac * max(s_flow(2), 0.0_dp) 
  end if  
  if (box%cc(i, j, k-1, i_eps) <= med) then
    s_flow(3) = dt * (box%fc(IJK, 3, flux_elec) - &
        box%fc(i, j, k+1, 3, flux_elec)) 
    box%fc(IJK, 3, sigma_rhs) = box%fc(IJK, 3, sigma_rhs) + fac * max(s_flow(3), 0.0_dp) 
  end if 
end if

 

#endif

       if (ST_update_ions) then
#if $D == 2
          f_ion = (box%fc(i, j, 2, flux_ion) - box%fc(i, j+1, 2, flux_ion) + &
               rfac(1) * box%fc(i, j, 1, flux_ion) - &
               rfac(2) * box%fc(i+1, j, 1, flux_ion)) * inv_dr * dt
#elif $D == 3
          f_ion = (sum(box%fc(i, j, k, 1:3, flux_ion)) - &
               box%fc(i+1, j, k, 1, flux_ion) - &
               box%fc(i, j+1, k, 2, flux_ion) - &
               box%fc(i, j, k+1, 3, flux_ion)) * inv_dr * dt
#endif
       end if

       ! Source term
       src = fld * mu * box%cc(IJK, i_electron) * (alpha - eta)
       if (photoi_enabled .and. box%cc(IJK, i_eps) < ST_epsilon_die) then
         src = src + box%cc(IJK, i_photo)
       else if (photoi_enabled .and. box%cc(IJK, i_eps) == ST_epsilon_die) then
         do n = 1, 4    
           call system_clock(time, count_rate)         
           call random_seed(time)
           call random_number(rand)
           dir = int(1+rand*(a$D_num_neighbors-1))
           ix = a$D_neighb_dix(:, dir)
#if $D == 2
           if (box%cc(i+ix(1), j+ix(2), i_eps) < ST_epsilon_die) then
             flag = 1
             exit
           end if    
#elif $D == 3
           if (box%cc(i+ix(1), j+ix(2), k+ix(3), i_eps) < ST_epsilon_die) then
            flag = 1 
             exit        
           end if  
#endif
         end do
       else if (box%cc(IJK, i_eps) > med) then
           src = 0.0_dp
       end if


       if (i_step == 1 .and. ST_output_src_term) then
          if (ST_output_src_decay_rate < 0) then
             ! No time-averaging
             box%cc(IJK, i_src) = src * a$D_heaviside(box%cc(IJK, i_eps), ST_epsilon_die)
          else
             ! Approximate exp(-t*rate) as 1-t*rate
             box%cc(IJK, i_src) = ((1 - dt * ST_output_src_decay_rate) * &
                  box%cc(IJK, i_src) + src * dt )* a$D_heaviside(box%cc(IJK, i_eps), ST_epsilon_die)
          end if
       end if

       ! Convert to density
       src  = src * dt


       if (flag == 1) then ! Add photoemission effect
#if $D == 2
         box%cc(i+ix(1), j+ix(2), i_electron) = box%cc(i+ix(1), j+ix(2), i_electron) + &
                                              box%cc(IJK, i_photo) * ST_phe_yield * dt
         box%fc(max(i, i+ix(1)), max(j, j+ix(2)), 1+(dir-1)/2, sigma_rhs) = &
               box%fc(max(i, i+ix(1)), max(j, j+ix(2)), 1+(dir-1)/2, sigma_rhs) - box%cc(IJK, i_photo) * ST_phe_yield * box%dr * dt
#elif $D == 3
         box%cc(i+ix(1), j+ix(2), k+ix(3), i_electron) = box%cc(i+ix(1), j+ix(2), k+ix(3), i_electron) + &
                                                       box%cc(IJK, i_photo) * ST_phe_yield * dt
#endif
       else
         if (box%cc(IJK, i_eps) < ST_epsilon_die) then
           ! Add flux and source term
           box%cc(IJK, i_electron) = box%cc(IJK, i_electron) + f_elec + src
           box%cc(IJK, i_pos_ion)  = box%cc(IJK, i_pos_ion) + f_ion + src
         end if
       end if

       ! Possible determine new time step restriction
       if (i_step == 2) then
#if $D == 2
          fld_vec(1) = 0.5_dp * (box%fc(IJK, 1, electric_fld) + &
               box%fc(i+1, j, 1, electric_fld))
          fld_vec(2) = 0.5_dp * (box%fc(IJK, 2, electric_fld) + &
               box%fc(i, j+1, 2, electric_fld))
#elif $D == 3
          fld_vec(1) = 0.5_dp * (box%fc(IJK, 1, electric_fld) + &
               box%fc(i+1, j, k, 1, electric_fld))
          fld_vec(2) = 0.5_dp * (box%fc(IJK, 2, electric_fld) + &
               box%fc(i, j+1, k, 2, electric_fld))
          fld_vec(3) = 0.5_dp * (box%fc(IJK, 3, electric_fld) + &
               box%fc(i, j, k+1, 3, electric_fld))
#endif

          diff = LT_get_col(ST_td_tbl, i_diffusion, fld)

          ! CFL condition
          dt_cfl = 1.0_dp/sum(abs(fld_vec * max(mu, abs(ST_ion_mobility), epsilon(1.0_dp))) * inv_dr)
          ! Diffusion condition
          dt_dif = box%dr**2 / max(2 * $D * max(diff, ST_ion_diffusion), epsilon(1.0_dp))
          ! Dielectric relaxation time
          dt_drt = UC_eps0 / (UC_elem_charge * (mu + abs(ST_ion_mobility)) * &
               max(box%cc(IJK, i_electron), epsilon(1.0_dp)))

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

    end do; CLOSE_DO
  end subroutine update_solution

  !> For each box that gets refined, set data on its children using this routine
  subroutine prolong_to_new_boxes(tree, ref_info)
    use m_a$D_prolong
    type(a$D_t), intent(inout)    :: tree
    type(ref_info_t), intent(in)  :: ref_info
    integer                       :: lvl, i, id, p_id

    !$omp parallel private(lvl, i, id, p_id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          p_id = tree%boxes(id)%parent
          
          call define_DI(tree%boxes(id)) ! For dielectric re-shaping
          call a$D_prolong_copy_fc(tree%boxes(p_id), tree%boxes(id), sigma_rhs, i_eps = i_eps, med = med)
          call a$D_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_electron, i_eps = i_eps, med = med)
          call a$D_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_pos_ion, i_eps = i_eps, med = med)
          call a$D_prolong_sparse(tree%boxes(p_id), tree%boxes(id), i_phi, i_eps = i_eps)
          if (photoi_enabled) then
             call a$D_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_photo, i_eps = i_eps, med = med)
          end if
          if (ST_output_src_term) then
             call a$D_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_src, i_eps = i_eps, med = med)
          end if
          call define_DI(tree%boxes(id))
       end do
       !$omp end do

       !$omp do
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a$D_gc_box(tree%boxes, id, i_eps, i_electron, &
               a$D_gc_interp_lim, a$D_bc_neumann_zero, med = med)
          call a$D_gc_box(tree%boxes, id, i_eps, i_pos_ion, &
               a$D_gc_interp_lim, a$D_bc_neumann_zero, med = med)
          call a$D_gc_box(tree%boxes, id, i_eps, i_phi, &
               mg%sides_rb, mg%sides_bc, med = med)
          if (photoi_enabled) then
             call a$D_gc_box(tree%boxes, id, i_eps, i_photo, &
               a$D_gc_interp_lim, a$D_bc_neumann_zero, med = med)
          end if
          if (ST_output_src_term) then
             call a$D_gc_box(tree%boxes, id, i_eps, i_src, &
               a$D_gc_interp_lim, a$D_bc_neumann_zero, med = med)
          end if
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine prolong_to_new_boxes

  subroutine write_log_file(tree, filename, out_cnt, dir)
    type(a$D_t), intent(in)      :: tree
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: out_cnt
    character(len=*), intent(in) :: dir
    character(len=ST_slen)       :: fname
    character(len=30), save      :: fmt
    integer, parameter           :: my_unit = 123
    real(dp)                     :: velocity, dt
    real(dp), save               :: prev_pos($D) = 0
    real(dp)                     :: sum_elec, sum_pos_ion
    real(dp)                     :: max_elec, max_field, max_Er
    type(a$D_loc_t)              :: loc_elec, loc_field, loc_Er

    call a$D_prepend_directory(filename, dir, fname)

    call a$D_tree_sum_cc(tree, i_electron, sum_elec)
    call a$D_tree_sum_cc(tree, i_pos_ion, sum_pos_ion)
    call a$D_tree_max_cc(tree, i_electron, max_elec, loc_elec)
    call a$D_tree_max_cc(tree, i_electric_fld, max_field, loc_field)
    call a$D_tree_max_fc(tree, 1, electric_fld, max_Er, loc_Er)

    dt = ST_dt_safety_factor * minval(ST_dt_matrix)

    if (out_cnt == 1) then
       open(my_unit, file=trim(fname), action="write")
#if $D == 2
       write(my_unit, *) "# it time dt v sum(n_e) sum(n_i) ", &
            "max(E) x y max(n_e) x y max(E_r) x y wc_time n_cells min(dx) highest(lvl)"
       fmt = "(I6,15E16.8,I12,1E16.8,I3)"
#elif $D == 3
       write(my_unit, *) "# it time dt v sum(n_e) sum(n_i) ", &
            "max(E) x y z max(n_e) x y z wc_time n_cells min(dx) highest(lvl)"
       fmt = "(I6,14E16.8,I12,1E16.8,I3)"
#endif
       close(my_unit)

       ! Start with velocity zero
       prev_pos = a$D_r_loc(tree, loc_field)
    end if

    velocity = norm2(a$D_r_loc(tree, loc_field) - prev_pos) / ST_dt_output
    prev_pos = a$D_r_loc(tree, loc_field)

    open(my_unit, file=trim(fname), action="write", &
         position="append")
#if $D == 2
    write(my_unit, fmt) out_cnt, ST_time, dt, velocity, sum_elec, &
         sum_pos_ion, max_field, a$D_r_loc(tree, loc_field), max_elec, &
         a$D_r_loc(tree, loc_elec), max_Er, a$D_r_loc(tree, loc_Er), &
         wc_time, a$D_num_cells_used(tree), a$D_min_dr(tree),tree%highest_lvl
#elif $D == 3
    write(my_unit, fmt) out_cnt, ST_time, dt, velocity, sum_elec, &
         sum_pos_ion, max_field, a$D_r_loc(tree, loc_field), max_elec, &
         a$D_r_loc(tree, loc_elec), wc_time, a$D_num_cells_used(tree), a$D_min_dr(tree),tree%highest_lvl
#endif
    close(my_unit)

  end subroutine write_log_file
  

  subroutine print_status()
    write(*, "(F7.2,A,I0,A,E10.3,A,E10.3,A,E10.3,A,E10.3)") &
             100 * ST_time / ST_end_time, "% it: ", it, &
             " t:", ST_time, " dt:", ST_dt, " wc:", wc_time, &
             " ncell:", real(a$D_num_cells_used(tree), dp)
  end subroutine print_status


end program streamer_$Dd
