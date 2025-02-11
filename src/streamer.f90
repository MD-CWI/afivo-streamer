!> Program to perform streamer simulations with AMR
program streamer
#include "../afivo/src/cpp_macros.h"
  use iso_fortran_env, only: error_unit
  use m_config
  use m_af_all
  use m_streamer
  use m_field
  use m_init_cond
  use m_refine
  use m_photoi
  use m_chemistry
  use m_gas
  use m_coupling
  use m_fluid
  use m_dt
  use m_types
  use m_user_methods
  use m_output
  use m_dielectric
  use m_units_constants
  use m_model
  use m_analysis
  use omp_lib

  implicit none

  integer, parameter        :: int8       = selected_int_kind(18)
  integer, parameter        :: max_attemps_per_time_step = 10
  integer, parameter        :: datfile_version = 30
  integer(int8)             :: t_start, t_current, count_rate
  real(dp)                  :: wc_time = 0.0_dp, inv_count_rate
  real(dp)                  :: time_last_print, time_last_output
  integer                   :: i, it, n, coord_type, box_bytes
  integer                   :: n_steps_rejected
  integer, allocatable      :: ref_links(:, :)
  logical                   :: write_out
  real(dp)                  :: time, dt, dt_lim, photoi_prev_time
  real(dp)                  :: dt_gas_lim, dt_lim_step
  real(dp)                  :: fraction_steps_rejected
  real(dp)                  :: memory_limit_GB
  type(af_t)                :: tree_copy      ! Used when reading a tree from a file
  type(ref_info_t)          :: ref_info       ! Contains info about refinement changes
  integer                   :: output_cnt = 0 ! Number of output files written
  character(len=string_len) :: restart_from_file = undefined_str
  real(dp)                  :: max_field
  type(af_loc_t)            :: loc_field, loc_field_t0
  real(dp)                  :: pos_Emax(NDIM), pos_Emax_t0(NDIM)
  real(dp)                  :: breakdown_field_Td, current_output_dt
  real(dp)                  :: time_until_next_pulse
  real(dp)                  :: field_energy, field_energy_prev
  real(dp)                  :: tmp, field_energy_prev_time
  logical                   :: step_accepted, start_of_new_pulse

  ! To keep track of the computational cost of different parts
  real(dp) :: t1, t2, t3

  !> The configuration for the simulation
  type(CFG_t) :: cfg
  !> This contains the full grid information
  type(af_t)  :: tree

  call print_program_name()

  ! Parse command line configuration files and options
  call CFG_update_from_arguments(cfg)

  call CFG_add_get(cfg, "restart_from_file", restart_from_file, &
       "If set, restart simulation from a previous .dat file")

  memory_limit_GB = 4.0_dp**(NDIM-1) ! 1, 4, 16 GB for 1D, 2D, 3D
  call CFG_add_get(cfg, "memory_limit_GB", memory_limit_GB, &
       "Memory limit (GB)")

  call initialize_modules(cfg, tree, mg, restart_from_file /= undefined_str)

  call CFG_write(cfg, trim(output_name) // "_out.cfg", custom_first=.true.)

  call chemistry_get_breakdown_field(breakdown_field_Td, 1.0e3_dp)
  write(*, '(A,E12.4)') " Estimated breakdown field (Td): ", breakdown_field_Td

  ! Specify default methods for all the variables
  do i = 1, size(all_densities)
     call af_set_cc_methods(tree, all_densities(i), &
          bc_species, af_gc_interp_lim, ST_prolongation_method)
  end do

  if (.not. gas_constant_density) then
     if (gas_dynamics) then
        ! Let the gas density evolve in time
        call af_set_cc_methods(tree, i_gas_dens, &
             af_bc_neumann_zero, af_gc_interp, ST_prolongation_method)
     else
        ! The gas density is specified by a function
        call af_set_cc_methods(tree, i_gas_dens, &
             funcval=set_gas_density_from_user_function)
     end if
  end if

  do i = 1, tree%n_var_cell
     if (tree%cc_write_output(i) .and. .not. &
          (tree%has_cc_method(i) .or. i == i_phi)) then
        call af_set_cc_methods(tree, i, af_bc_neumann_zero, &
             af_gc_interp, ST_prolongation_method)
     end if
  end do

  ! These values are overwritten when restarting
  it               = 0
  time             = 0.0_dp ! Simulation time (all times are in s)
  global_time      = time
  photoi_prev_time = time   ! Time of last photoionization computation
  dt               = global_dt
  n_steps_rejected = 0
  fraction_steps_rejected = 0.0_dp
  pos_Emax_t0 = 0.0_dp ! Initial streamer position

  ! Initialize the tree (which contains all the mesh information)
  if (restart_from_file /= undefined_str) then
     tree_copy = tree           ! Store the settings
     call af_read_tree(tree, restart_from_file, read_sim_data)

     ! Restore some of the settings
     tree%cc_write_output = tree_copy%cc_write_output
     tree%cc_write_binary = tree_copy%cc_write_binary
     tree%fc_write_binary = tree_copy%fc_write_binary

     box_bytes = af_box_bytes(tree%n_cell, tree%n_var_cell, tree%n_var_face)
     tree%box_limit = nint(memory_limit_GB * 2.0_dp**30 / box_bytes)

     if (tree%n_cell /= ST_box_size) &
          error stop "restart_from_file: incompatible box size"

     if (tree%n_var_cell /= tree_copy%n_var_cell) then
        write(error_unit, *) "n_var_cell here:", tree_copy%n_var_cell
        write(error_unit, *) "n_var_cell file:", tree%n_var_cell
        error stop "restart_from_file: incompatible variable list"
     end if

     if (ST_use_dielectric) error stop "Restarting not support with dielectric"

     ! @todo more consistency checks

     ! This routine always needs to be called when using multigrid
     call mg_init(tree, mg)
  else
     if (ST_cylindrical) then
        coord_type = af_cyl
     else
        coord_type = af_xyz
     end if

     call af_init(tree, ST_box_size, ST_domain_origin + ST_domain_len, &
          ST_coarse_grid_size, periodic=ST_periodic, coord=coord_type, &
          r_min=ST_domain_origin, mem_limit_gb=memory_limit_GB)

     ! Set up the initial conditions
     call set_initial_conditions(tree, mg)

     ! Write initial output
     output_cnt = 0 ! Number of output files written
     call output_write(tree, output_cnt, 0.0_dp, write_sim_data)
  end if

  print *, "Simulation output: ", trim(output_name)
  print *, "Number of threads: ", af_get_max_threads()
  call af_print_info(tree)
  call af_stencil_print_info(tree)

  ! Initial wall clock time
  call system_clock(t_start, count_rate)
  inv_count_rate   = 1.0_dp / count_rate
  time_last_print  = -1e10_dp
  time_last_output = time

  call field_compute_energy(tree, field_energy_prev)
  field_energy_prev_time = time

  do
     it = it + 1
     if (time >= ST_end_time) exit

     if (associated(user_generic_method)) then
        call user_generic_method(tree, time)
     end if

     ! Initialize starting position of streamer
     if (ST_use_end_streamer_length .and. it == ST_initial_streamer_pos_steps_wait) then
        call af_tree_max_cc(tree, i_electric_fld, max_field, loc_field_t0)
        pos_Emax_t0 = af_r_loc(tree, loc_field_t0)
     end if

     ! Check if streamer length exceeds the defined maximal streamer length
     if (ST_use_end_streamer_length .and. it > ST_initial_streamer_pos_steps_wait) then
        call af_tree_max_cc(tree, i_electric_fld, max_field, loc_field)
        pos_Emax = af_r_loc(tree, loc_field)

        if (norm2(pos_Emax_t0 - pos_Emax) >= ST_end_streamer_length) then
           print *, "Streamer reached its desired length"
           exit
        end if
     end if

     ! Update wall clock time
     call system_clock(t_current)
     wc_time = (t_current - t_start) * inv_count_rate

     ! Every ST_print_status_interval, print some info about progress
     if (wc_time - time_last_print > output_status_delay) then
        call output_status(tree, time, wc_time, it, dt)
        time_last_print = wc_time
     end if

     time_until_next_pulse = field_pulse_period - modulo(time, field_pulse_period)

     if (abs(current_voltage) > 0.0_dp .or. &
          time_until_next_pulse < refine_prepulse_time) then
        current_output_dt = output_dt
        current_electrode_dx = refine_electrode_dx
     else
        current_output_dt = output_dt * output_dt_factor_pulse_off
        current_electrode_dx = electrode_derefine_factor*refine_electrode_dx
     end if

     ! Every output_dt, write output
     write_out = (time + dt >= time_last_output + current_output_dt)
     if (write_out) then
        ! Ensure that dt is non-negative, even when current_output_dt changes
        dt = max(0.0_dp, time_last_output + current_output_dt - time)
     end if

     ! Make sure to capture the start of the next pulse
     start_of_new_pulse = (dt >= time_until_next_pulse)
     if (start_of_new_pulse) then
        dt = max(time_until_next_pulse, dt_min)
     end if

     if (photoi_enabled .and. mod(it, photoi_per_steps) == 0) then
        t1 = omp_get_wtime()
        call photoi_set_src(tree, time - photoi_prev_time)
        t2 = omp_get_wtime()
        wc_time_photoi = wc_time_photoi + t2 - t1
        photoi_prev_time = time
     end if

     if (ST_use_electrode) then
        call set_electrode_densities(tree)
     end if

     ! Advance over dt, but make a copy first so that we can try again if dt was too large
     dt_lim = huge_real
     step_accepted = .false.
     do n = 1, max_attemps_per_time_step
        t1 = omp_get_wtime()
        call copy_current_state()
        t2 = omp_get_wtime()
        wc_time_copy_state = wc_time_copy_state + t2 - t1

        call af_advance(tree, dt, dt_lim_step, time, all_densities, &
             time_integrator, forward_euler)

        ! dt_lim is the minimum over all steps taken (including rejected ones)
        dt_lim = min(dt_lim, dt_lim_step)

        ! Check if dt was small enough for the new state
        step_accepted = (dt <= dt_lim_step)

        if (step_accepted) then
           exit
        else
           n_steps_rejected = n_steps_rejected + 1
           write (*, "(I0,A,I0,A,2E12.4,A,F6.2)") it, " Step rejected (#", &
                n_steps_rejected, ") (dt, dt_lim) = ", dt, dt_lim, &
                " frac", fraction_steps_rejected
           call output_status(tree, time, wc_time, it, dt)

           ! Go back to previous state and try with a smaller dt
           dt = dt_safety_factor * dt_lim_step
           time = global_time
           write_out = .false. ! Since we advance less far in time
           call restore_previous_state()
        end if
     end do

     ! The approximate fraction of rejected steps, over the last ~100
     fraction_steps_rejected = 0.99_dp * fraction_steps_rejected
     if (n > 1) fraction_steps_rejected = fraction_steps_rejected + 0.01_dp

     if (n == max_attemps_per_time_step+1) &
          error stop "All time steps were rejected"

     ! Update global variable based on current space-integrated data
     ST_global_rates = ST_global_rates + &
          sum(ST_current_rates(1:n_reactions, :), dim=2) * dt
     ST_global_JdotE = ST_global_JdotE + &
          sum(ST_current_JdotE(1, :)) * dt

     ! Estimate electric current according to Sato's equation V*I = sum(J.E),
     ! where J includes both the conduction current and the displacement
     ! current, see 10.1088/0022-3727/32/5/005.
     ! The latter is computed through the field energy, which contains
     ! some noise, so the current is only updated every N iterations.
     if (mod(it, current_update_per_steps) == 0) then
        call field_compute_energy(tree, field_energy)

        ! Time derivative of field energy
        tmp = (field_energy - field_energy_prev)/(time - field_energy_prev_time)
        field_energy_prev      = field_energy
        field_energy_prev_time = time

        ! Add J.E term
        if (abs(current_voltage) > 0.0_dp) then
           ST_global_JdotE_current = sum(ST_current_JdotE(1, :)) / current_voltage
           ST_global_displ_current = tmp/current_voltage
        else
           ST_global_JdotE_current = 0.0_dp
           ST_global_displ_current = 0.0_dp
        end if
     end if

     ! Make sure field is available for latest time state
     t1 = omp_get_wtime()
     call field_compute(tree, mg, 0, time, .true.)
     t2 = omp_get_wtime()
     wc_time_field = wc_time_field + t2 - t1

     if (gas_dynamics) then
        call coupling_add_fluid_source(tree, dt)

        ! Go back to time at beginning of step
        time = global_time

        call af_advance(tree, dt, dt_gas_lim, time, &
             gas_vars, time_integrator, gas_forward_euler)
        call coupling_update_gas_density(tree)
     else
        dt_gas_lim = dt_max
     end if

     ! Do not increase time step when many steps are rejected. Here we use a
     ! pretty arbitrary threshold of 0.1
     tmp = dt_max_growth_factor
     if (fraction_steps_rejected > 0.1_dp) tmp = 1.0_dp

     dt = min(tmp * global_dt, dt_safety_factor * min(dt_lim, dt_gas_lim))

     if (start_of_new_pulse) then
        ! Start a new pulse with a small time step
        dt = dt_min
        if (associated(user_new_pulse_conditions)) then
           call af_loop_box(tree, user_new_pulse_conditions)
        end if
     end if

     ! dt is modified when writing output, global_dt not
     global_dt   = dt
     global_time = time

     if (global_dt < dt_min) then
        write(error_unit, "(A,E12.4,A)") " Time step (dt =", global_dt, &
             ") getting too small"
        write(error_unit, *) "See the documentation on time integration"
        call output_status(tree, time, wc_time, it, dt)
        write_out = .true.
     end if

     t1 = omp_get_wtime()
     if (write_out) then
        output_cnt       = output_cnt + 1
        time_last_output = global_time
        call output_write(tree, output_cnt, wc_time, write_sim_data)
        if (ST_use_dielectric .and. surface_output) then
           call surface_write_output(tree, diel, [i_photon_flux, i_surf_dens], &
                ["photon_flux", "surf_dens  "], output_name, output_cnt)
        end if
     end if
     t2 = omp_get_wtime()
     wc_time_output = wc_time_output + t2 - t1

     if (global_dt < dt_min) error stop "dt too small"

     if (ST_cylindrical .and. ST_abort_axisymmetric_if_branching) then
        if (axisymmetric_is_branching(tree)) then
           error stop  "Branching detected, abort"
        end if
     end if

     if (mod(it, refine_per_steps) == 0) then
        ! Restrict species, for the ghost cells near refinement boundaries
        call af_restrict_tree(tree, all_densities)
        call af_gc_tree(tree, all_densities)

        if (gas_dynamics) then
           call af_restrict_tree(tree, gas_vars)
           call af_gc_tree(tree, gas_vars)
        end if

        if (ST_use_dielectric) then
           ! Make sure there are no refinement jumps across the dielectric
           call surface_get_refinement_links(diel, ref_links)
           call af_adjust_refinement(tree, refine_routine, ref_info, &
                refine_buffer_width, ref_links)
           call surface_update_after_refinement(tree, diel, ref_info)
        else
           call af_adjust_refinement(tree, refine_routine, ref_info, &
                refine_buffer_width)
        end if

        if (ref_info%n_add > 0 .or. ref_info%n_rm > 0) then
           ! Compute the field on the new mesh
           call field_compute(tree, mg, 0, time, .true.)

           ! Compute photoionization on new mesh
           if (photoi_enabled) then
              call photoi_set_src(tree, time - photoi_prev_time)
              photoi_prev_time = time
           end if
        end if
     end if

     t3 = omp_get_wtime()
     wc_time_refine = wc_time_refine + t3 - t2
  end do

  call output_status(tree, time, wc_time, it, dt)

  write(*, "(A)") "Computational cost breakdown (%)"
  write(*, "(7(A10))") "flux", "source", "copy", "field", "output", &
       "refine", "photoi"
  write(*, "(7(F10.2))") 1e2*wc_time_flux/wc_time, 1e2*wc_time_source/wc_time, &
       1e2*wc_time_copy_state/wc_time, 1e2*wc_time_field/wc_time, &
       1e2*wc_time_output/wc_time, 1e2*wc_time_refine/wc_time, &
       1e2*wc_time_photoi/wc_time

contains

  subroutine initialize_modules(cfg, tree, mg, restart)
    use m_user
    use m_table_data
    use m_transport_data

    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout)  :: tree
    type(mg_t), intent(inout)  :: mg
    logical, intent(in)        :: restart

    call model_initialize(cfg)
    call user_initialize(cfg, tree)
    call dt_initialize(cfg)
    global_dt = dt_min

    call table_data_initialize(cfg)
    call gas_initialize(tree, cfg)
    call transport_data_initialize(cfg)
    call chemistry_initialize(tree, cfg)
    call ST_initialize(tree, cfg, NDIM)
    call photoi_initialize(tree, cfg)
    call refine_initialize(cfg)
    call field_initialize(tree, cfg, mg)
    call init_cond_initialize(tree, cfg)
    call output_initialize(tree, cfg)
    call dielectric_initialize(tree, cfg)

    call output_initial_summary(tree)

  end subroutine initialize_modules

  subroutine set_initial_conditions(tree, mg)
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(inout) :: mg
    integer                   :: n, lvl, i, id, n_surface_variables

    ! Refine up to maximal grid spacing
    do lvl = 1, af_max_lvl-1
       if (all(af_lvl_dr(tree, lvl) <= refine_max_dx)) exit
    end do
    call af_refine_up_to_lvl(tree, lvl)

    ! Set initial conditions
    call af_loop_box(tree, init_cond_set_box)
    if (associated(user_initial_conditions)) then
       call af_loop_box(tree, user_initial_conditions)
    else if (ST_use_dielectric) then
       error stop "use_dielectric requires user_initial_conditions to be set"
    end if

    ! This placement is so that users can set epsilon before the coarse grid
    ! solver is constructed
    call mg_init(tree, mg)

    if (ST_use_dielectric) then
       ! To store the photon flux (single state), surface charge (multiple
       ! states) and a copy of the surface charge
       n_surface_variables = af_advance_num_steps(time_integrator) + 2
       call surface_initialize(tree, i_eps, diel, n_surface_variables)
    end if

    do n = 1, 100
       call field_compute(tree, mg, 0, time, .false.)

       if (ST_use_dielectric) then
          ! Make sure there are no refinement jumps across the dielectric
          call surface_get_refinement_links(diel, ref_links)
          call af_adjust_refinement(tree, refine_routine, ref_info, &
               refine_buffer_width, ref_links)
          call surface_update_after_refinement(tree, diel, ref_info)
       else
          call af_adjust_refinement(tree, refine_routine, ref_info, &
               refine_buffer_width)
       end if

       ! Set initial conditions on newly added boxes
       do lvl = 1, size(ref_info%lvls)
          !$omp parallel do private(id)
          do i = 1, size(ref_info%lvls(lvl)%add)
             id = ref_info%lvls(lvl)%add(i)
             call init_cond_set_box(tree%boxes(id))
             if (associated(user_initial_conditions)) &
                  call user_initial_conditions(tree%boxes(id))
          end do
          !$omp end parallel do
       end do

       if (ref_info%n_add == 0) exit
    end do

  end subroutine set_initial_conditions

  subroutine write_sim_data(my_unit)
    integer, intent(in) :: my_unit

    write(my_unit) datfile_version

    write(my_unit) it
    write(my_unit) output_cnt
    write(my_unit) time
    write(my_unit) global_time
    write(my_unit) photoi_prev_time
    write(my_unit) global_dt

    write(my_unit) ST_global_rates
    write(my_unit) ST_global_JdotE
    write(my_unit) fraction_steps_rejected
  end subroutine write_sim_data

  subroutine read_sim_data(my_unit)
    integer, intent(in) :: my_unit
    integer             :: version

    read(my_unit) version
    if (version /= datfile_version) error stop "Different datfile version"

    read(my_unit) it
    read(my_unit) output_cnt
    read(my_unit) time
    read(my_unit) global_time
    read(my_unit) photoi_prev_time
    read(my_unit) global_dt

    read(my_unit) ST_global_rates
    read(my_unit) ST_global_JdotE
    read(my_unit) fraction_steps_rejected

    dt = global_dt
  end subroutine read_sim_data

  subroutine print_program_name()
    print *, "            __ _                  _____ _                                      "
    print *, "     /\    / _(_)                / ____| |                                     "
    print *, "    /  \  | |_ ___   _____ _____| (___ | |_ _ __ ___  __ _ _ __ ___   ___ _ __ "
    print *, "   / /\ \ |  _| \ \ / / _ \______\___ \| __| '__/ _ \/ _` | '_ ` _ \ / _ \ '__|"
    print *, "  / ____ \| | | |\ V / (_) |     ____) | |_| | |  __/ (_| | | | | | |  __/ |   "
    print *, " /_/    \_\_| |_| \_/ \___/     |_____/ \__|_|  \___|\__,_|_| |_| |_|\___|_|   "
    print *, "                                                                               "
  end subroutine print_program_name

  !> Set species boundary conditions at the electrode
  subroutine set_electrode_densities(tree)
    type(af_t), intent(inout) :: tree

    call af_loop_box(tree, electrode_species_bc)
  end subroutine set_electrode_densities

  !> Set species boundary conditions at the electrode
  !> @todo Set Neumann boundary conditions in the actual flux computation
  subroutine electrode_species_bc(box)
    type(box_t), intent(inout) :: box
    real(dp)                   :: lsf_nb(2*NDIM), dens_nb(2*NDIM)
    integer                    :: nc, IJK

    if (iand(box%tag, mg_lsf_box) == 0) return

    nc = box%n_cell
    do KJI_DO(1, nc)
       if (box%cc(IJK, i_lsf) < 0) then
          ! Set all species densities to zero
          box%cc(IJK, all_densities) = 0.0_dp

#if NDIM == 1
          lsf_nb = [box%cc(i-1, i_lsf), &
               box%cc(i+1, i_lsf)]
#elif NDIM == 2
          lsf_nb = [box%cc(i-1, j, i_lsf), &
               box%cc(i+1, j, i_lsf), &
               box%cc(i, j-1, i_lsf), &
               box%cc(i, j+1, i_lsf)]
#elif NDIM == 3
          lsf_nb = [box%cc(i-1, j, k, i_lsf), &
               box%cc(i+1, j, k, i_lsf), &
               box%cc(i, j-1, k, i_lsf), &
               box%cc(i, j+1, k, i_lsf), &
               box%cc(i, j, k-1, i_lsf), &
               box%cc(i, j, k+1, i_lsf)]
#endif

          if (any(lsf_nb > 0) .and. &
               associated(bc_species, af_bc_neumann_zero)) then
             ! At the boundary of the electrode
#if NDIM == 1
             dens_nb = [box%cc(i-1, i_electron), &
                  box%cc(i+1, i_electron)]
#elif NDIM == 2
             dens_nb = [box%cc(i-1, j, i_electron), &
                  box%cc(i+1, j, i_electron), &
                  box%cc(i, j-1, i_electron), &
                  box%cc(i, j+1, i_electron)]
#elif NDIM == 3
             dens_nb = [box%cc(i-1, j, k, i_electron), &
                  box%cc(i+1, j, k, i_electron), &
                  box%cc(i, j-1, k, i_electron), &
                  box%cc(i, j+1, k, i_electron), &
                  box%cc(i, j, k-1, i_electron), &
                  box%cc(i, j, k+1, i_electron)]
#endif

             ! Set electron density to average of outside neighbors
             box%cc(IJK, i_electron) = sum(dens_nb, mask=(lsf_nb > 0)) / &
                  count(lsf_nb > 0)
             ! Set first positive ion density for charge neutrality
             box%cc(IJK, i_1pos_ion) = box%cc(IJK, i_electron)
          end if
       end if
    end do; CLOSE_DO
  end subroutine electrode_species_bc

  !> Copy current state for densities and field
  subroutine copy_current_state()
    integer :: n_states

    n_states = af_advance_num_steps(time_integrator)
    call af_tree_copy_ccs(tree, all_densities, all_densities+n_states)

    ! Copy potential
    call af_tree_copy_cc(tree, i_phi, i_phi+1)

    if (ST_use_dielectric) then
       call surface_copy_variable(diel, i_surf_dens, i_surf_dens+n_states)
    end if
  end subroutine copy_current_state

  !> Restore state for densities and field
  subroutine restore_previous_state()
    integer :: n_states

    n_states = af_advance_num_steps(time_integrator)

    call af_tree_copy_ccs(tree, all_densities+n_states, all_densities)

    ! Copy potential and compute field again
    call af_tree_copy_cc(tree, i_phi+1, i_phi)
    call field_from_potential(tree, mg)

    if (ST_use_dielectric) then
       call surface_copy_variable(diel, i_surf_dens+n_states, i_surf_dens )
    end if
  end subroutine restore_previous_state

  !> A wrapper routine to set the gas density in a box by a user-defined
  !> function
  subroutine set_gas_density_from_user_function(box, iv)
    type(box_t), intent(inout) :: box !< Box to fill values in
    integer, intent(in)        :: iv  !< Index of variable
    integer                    :: IJK, nc
    nc = box%n_cell

    do KJI_DO(0, nc+1)
       box%cc(IJK, iv) = user_gas_density(box, IJK)
    end do; CLOSE_DO
  end subroutine set_gas_density_from_user_function

  logical function axisymmetric_is_branching(tree)
    type(af_t), intent(in) :: tree
    real(dp)               :: max_field, axis_field, rmin(NDIM), rmax(NDIM)
    type(af_loc_t)         :: loc_field

    ! Get location of max(E) in full domain
    call af_tree_max_cc(tree, i_electric_fld, max_field, loc_field)

    ! Compare with maximum in region near axis
    rmin = 0.0_dp
    rmax(1) = 1e-2_dp * ST_domain_len(1)
    rmax(2:) = ST_domain_len(2:)
    call analysis_max_var_region(tree, i_electric_fld, rmin, &
         rmax, axis_field, loc_field)

    axisymmetric_is_branching = (max_field > 1.1_dp * axis_field)

    if (axisymmetric_is_branching) print *, max_field, axis_field
  end function axisymmetric_is_branching

end program streamer
