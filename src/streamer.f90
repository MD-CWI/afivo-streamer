#include "../afivo/src/cpp_macros.h"
!> Program to perform streamer simulations with AMR
program streamer

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
  use m_fluid_lfa
  use m_dt
  use m_types
  use m_user_methods
  use m_output
  use m_circuit

  implicit none

  integer, parameter        :: int8       = selected_int_kind(18)
  integer(int8)             :: t_start, t_current, count_rate
  real(dp)                  :: wc_time, inv_count_rate
  real(dp)                  :: time_last_print, time_last_output
  integer                   :: i, it, coord_type, box_bytes
  logical                   :: write_out, evolve_electrons
  real(dp)                  :: time, dt, dt_lim, photoi_prev_time
  real(dp)                  :: dt_gas_lim
  real(dp)                  :: memory_limit_GB = 16.0_dp
  type(CFG_t)               :: cfg            ! The configuration for the simulation
  type(af_t)                :: tree           ! This contains the full grid information
  type(af_t)                :: tree_copy      ! Used when reading a tree from a file
  type(ref_info_t)          :: ref_info       ! Contains info about refinement changes
  integer                   :: output_cnt = 0 ! Number of output files written
  character(len=string_len) :: restart_from_file = undefined_str
  real(dp)                  :: max_field, initial_streamer_pos
  type(af_loc_t)            :: loc_field, loc_field_initial
  real(dp), dimension(NDIM) :: loc_field_coord, loc_field_initial_coord

  call print_program_name()

  ! Parse command line configuration files and options
  call CFG_update_from_arguments(cfg)

  call CFG_add_get(cfg, "restart_from_file", restart_from_file, &
       "If set, restart simulation from a previous .dat file")
  call CFG_add_get(cfg, "memory_limit_GB", memory_limit_GB, &
       "Memory limit (GB)")

  call initialize_modules(cfg, tree, mg, restart_from_file /= undefined_str)

  call CFG_write(cfg, trim(output_name) // "_out.cfg", custom_first=.true.)

  ! Specify default methods for all the variables
  do i = n_gas_species+1, n_species
     call af_set_cc_methods(tree, species_itree(i), &
          bc_species, af_gc_interp_lim, ST_prolongation_method)
  end do

  if (.not. gas_constant_density) then
     call af_set_cc_methods(tree, i_gas_dens, &
          af_bc_neumann_zero, af_gc_interp, ST_prolongation_method)
  end if

  do i = 1, tree%n_var_cell
     if (tree%cc_write_output(i) .and. .not. &
          (tree%has_cc_method(i) .or. i == i_phi)) then
        call af_set_cc_methods(tree, i, af_bc_neumann_zero, &
             af_gc_interp, ST_prolongation_method)
     end if
  end do

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
        print *, "n_var_cell here:", tree_copy%n_var_cell
        print *, "n_var_cell file:", tree%n_var_cell
        error stop "restart_from_file: incompatible variable list"
     end if

     ! @todo more consistency checks

     ! This routine always needs to be called when using multigrid
     call mg_init(tree, mg)
  else
     time             = 0.0_dp ! Simulation time (all times are in s)
     global_time      = time
     photoi_prev_time = time   ! Time of last photoionization computation
     dt               = global_dt
     initial_streamer_pos = 0.0_dp ! Initial streamer position

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

  ! Initial wall clock time
  call system_clock(t_start, count_rate)
  inv_count_rate   = 1.0_dp / count_rate
  time_last_print  = -1e10_dp
  time_last_output = time

  do it = 1, huge(1)-1
     if (ST_use_end_time .and. time >= ST_end_time) exit

     ! Initialize starting position of streamer
     if (ST_use_end_streamer_length .and. it == ST_initial_streamer_pos_steps_wait) then
        call af_tree_max_cc(tree, i_electric_fld, max_field, loc_field_initial)
        loc_field_initial_coord = af_r_loc(tree, loc_field_initial)
     end if

     ! Check if streamer length exceeds the defined maximal streamer length
     if (ST_use_end_streamer_length .and. it > ST_initial_streamer_pos_steps_wait) then
        call af_tree_max_cc(tree, i_electric_fld, max_field, loc_field)
        loc_field_coord = af_r_loc(tree, loc_field)
        if (NORM2(loc_field_initial_coord - loc_field_coord) >= ST_end_streamer_length) exit
     end if

     ! Update wall clock time
     call system_clock(t_current)
     wc_time = (t_current - t_start) * inv_count_rate

     ! Every ST_print_status_interval, print some info about progress
     if (wc_time - time_last_print > output_status_delay) then
        call output_status(tree, time, wc_time, it, dt)
        time_last_print = wc_time
     end if

     ! Every output_dt, write output
     if (time + dt >= time_last_output + output_dt) then
        write_out        = .true.
        dt               = time_last_output + output_dt - time
        time_last_output = time_last_output + output_dt
        output_cnt       = output_cnt + 1
     else
        write_out = .false.
     end if

     evolve_electrons = .true.
     if (associated(user_evolve_electrons)) &
          evolve_electrons = user_evolve_electrons(tree, time)

     if (evolve_electrons) then
        if (photoi_enabled .and. mod(it, photoi_per_steps) == 0) then
           call photoi_set_src(tree, time - photoi_prev_time)
           photoi_prev_time = time
        end if

        if (ST_use_electrode) then
           call set_electrode_densities(tree)
        end if

        call af_advance(tree, dt, dt_lim, time, &
             species_itree(n_gas_species+1:n_species), &
             af_heuns_method, forward_euler)

        ! Make sure field is available for latest time state
        call field_compute(tree, mg, 0, time, .true.)
        if (compute_power_density) call compute_total_energy_density(tree, dt)
        if (gas_dynamics) call coupling_add_fluid_source(tree, dt)
        if (circuit_used) call circuit_update(tree, dt)
     else
        dt_lim = dt_max
     end if

     if (gas_dynamics) then
        ! Go back to time at beginning of step
        time = global_time

        call af_advance(tree, dt, dt_gas_lim, time, &
             gas_vars, time_integrator, gas_forward_euler)
        call coupling_update_gas_density(tree)
     else
        dt_gas_lim = dt_max
     end if

     ! If neither electrons or the gas is evolved, make sure time is increased
     if (.not. (evolve_electrons .or. gas_dynamics)) then
        time = time + dt
     end if

     ! dt is modified when writing output, global_dt not
     dt          = min(dt_lim, dt_gas_lim)
     global_dt   = dt_lim
     global_time = time

     if (global_dt < dt_min) then
        write(*, "(A,E12.4,A)") " Time step (dt =", global_dt, &
             ") getting too small"
        print *, "See the documentation on time integration"
        call output_status(tree, time, wc_time, it, dt)
        if (.not. write_out) then
           write_out = .true.
           output_cnt = output_cnt + 1
        end if
     end if

     if (write_out) then
        call output_write(tree, output_cnt, wc_time, write_sim_data)
     end if

     if (global_dt < dt_min) error stop "dt too small"

     if (mod(it, refine_per_steps) == 0) then
        ! Restrict species, for the ghost cells near refinement boundaries
        call af_restrict_tree(tree, species_itree(n_gas_species+1:n_species))
        call af_gc_tree(tree, species_itree(n_gas_species+1:n_species))

        if (gas_dynamics) then
           call af_restrict_tree(tree, gas_vars)
           call af_gc_tree(tree, gas_vars)
        end if

        if (associated(user_refine)) then
           call af_adjust_refinement(tree, user_refine, ref_info, &
                refine_buffer_width)
        else
           call af_adjust_refinement(tree, refine_routine, ref_info, &
                refine_buffer_width)
        end if

        if (ref_info%n_add > 0 .or. ref_info%n_rm > 0) then
           ! Compute the field on the new mesh
           call field_compute(tree, mg, 0, time, .true.)
        end if
     end if
  end do

  call output_status(tree, time, wc_time, it, dt)

contains

  subroutine initialize_modules(cfg, tree, mg, restart)
    use m_user
    use m_table_data
    use m_transport_data

    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout)  :: tree
    type(mg_t), intent(inout)  :: mg
    logical, intent(in)        :: restart

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
    call circuit_initialize(tree, cfg, restart)
    call init_cond_initialize(tree, cfg)
    call output_initialize(tree, cfg)

    call output_initial_summary()

  end subroutine initialize_modules

  subroutine set_initial_conditions(tree, mg)
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(inout) :: mg
    integer                   :: n, lvl, i, id

    call af_loop_box(tree, init_cond_set_box)
    if (associated(user_initial_conditions)) then
       call af_loop_box(tree, user_initial_conditions)
    else if (ST_use_dielectric) then
       error stop "use_dielectric requires user_initial_conditions to be set"
    end if

    ! This placement is so that users can set epsilon before the coarse grid
    ! solver is constructed
    call mg_init(tree, mg)

    do n = 1, 100
       call field_compute(tree, mg, 0, time, .false.)

       if (associated(user_refine)) then
          call af_adjust_refinement(tree, user_refine, ref_info, &
               refine_buffer_width)
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

    write(my_unit) output_cnt
    write(my_unit) time
    write(my_unit) global_time
    write(my_unit) photoi_prev_time
    write(my_unit) global_dt
  end subroutine write_sim_data

  subroutine read_sim_data(my_unit)
    integer, intent(in) :: my_unit

    read(my_unit) output_cnt
    read(my_unit) time
    read(my_unit) global_time
    read(my_unit) photoi_prev_time
    read(my_unit) global_dt

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

    if (box%tag /= mg_lsf_box) return

    nc = box%n_cell
    do KJI_DO(1, nc)
       if (box%cc(IJK, i_lsf) < 0) then
          ! Set all species densities to zero
          box%cc(IJK, species_itree(n_gas_species+1:n_species)) = 0.0_dp

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
             dens_nb = [box%cc(i-1, i_lsf), &
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

end program streamer
