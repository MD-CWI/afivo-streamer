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
  use m_fluid_lfa
  use m_dt
  use m_types
  use m_user_methods
  use m_output

  implicit none

  integer, parameter        :: int8       = selected_int_kind(18)
  integer(int8)             :: t_start, t_current, count_rate
  real(dp)                  :: wc_time, inv_count_rate
  real(dp)                  :: time_last_print, time_last_output
  integer                   :: i, it, coord_type, box_bytes
  logical                   :: write_out
  real(dp)                  :: time, dt, dt_lim, photoi_prev_time
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
  ! Circuit stuff. fine lets do it with global variables
  real(dp), dimension(2) :: r_top_electrode
  real(dp) :: old_calculated_voltage_at_top_electrode = 0.0
  real(dp) :: new_calculated_voltage_at_top_electrode = 0.0
  real(dp) :: capacitor_voltage = 0.0
  real(dp) :: capacitor_capacitance = 200e-12
  real(dp) :: top_electrode_voltage = 0.0
  real(dp) :: top_electrode_self_capacitance = 0.0
  real(dp) :: resistor = 300

  top_electrode_self_capacitance = 0.1 * capacitor_capacitance

  call print_program_name()

  ! Parse command line configuration files and options
  call CFG_update_from_arguments(cfg)

  call CFG_add_get(cfg, "restart_from_file", restart_from_file, &
       "If set, restart simulation from a previous .dat file")
  call CFG_add_get(cfg, "memory_limit_GB", memory_limit_GB, &
       "Memory limit (GB)")

  call initialize_modules(cfg, tree, mg)

  call CFG_write(cfg, trim(output_name) // "_out.cfg", custom_first=.true.)

  ! Specify default methods for all the variables
  do i = n_gas_species+1, n_species
     call af_set_cc_methods(tree, species_itree(i), &
          af_bc_neumann_zero, af_gc_interp_lim, ST_prolongation_method)
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

     if (photoi_enabled .and. mod(it, photoi_per_steps) == 0) then
        call photoi_set_src(tree, time - photoi_prev_time)
        photoi_prev_time = time
     end if

     call af_advance(tree, dt, dt_lim, time, &
          species_itree(n_gas_species+1:n_species), &
          af_heuns_method, forward_euler)

     ! Change potentials at boundaries due to circuit
# if NDIM == 2
      call add_circuit_effect(tree)
# endif
     ! Make sure field is available for latest time state
     call field_compute(tree, mg, 0, time, .true.)

     ! dt is modified when writing output, global_dt not
     global_dt   = dt_lim
     dt          = dt_lim
     global_time = time

     if (global_dt < dt_min) then
        print *, "ST_dt getting too small, instability?", global_dt
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

  subroutine initialize_modules(cfg, tree, mg)
    use m_user
    use m_table_data
    use m_transport_data

    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout)  :: tree
    type(mg_t), intent(inout)  :: mg

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
    call init_cond_initialize(cfg)
    call output_initialize(tree, cfg)

    call output_initial_summary()

  end subroutine initialize_modules

  subroutine set_initial_conditions(tree, mg)
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(inout) :: mg
    integer                   :: n

    do n = 1, 100
       call af_loop_box(tree, init_cond_set_box)

       if (associated(user_initial_conditions)) then
          call af_loop_box(tree, user_initial_conditions)
       else if (ST_use_dielectric) then
          error stop "use_dielectric requires user_initial_conditions to be set"
       end if

       ! This placement is so that users can set epsilon before the coarse grid
       ! solver is constructed
       if (n == 1) call mg_init(tree, mg)

       call field_compute(tree, mg, 0, time, .false.)

       if (associated(user_refine)) then
          call af_adjust_refinement(tree, user_refine, ref_info, &
               refine_buffer_width)
       else
          call af_adjust_refinement(tree, refine_routine, ref_info, &
               refine_buffer_width)
       end if

       if (ref_info%n_add == 0) exit
    end do

    ! Initialize capacitor voltage to be the same as applied voltage
    capacitor_voltage = -field_get_amplitude(tree, time) * ST_domain_len(2)
    print *, "Initial capacitor voltage: ", capacitor_voltage
    r_top_electrode(1) = 0
    r_top_electrode(2) = ST_domain_len(2)
    top_electrode_voltage = capacitor_voltage

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


   subroutine add_circuit_effect(tree)
      type(af_t), intent(inout)     :: tree
      real(dp) :: delta_V_electrode_dt
      real(dp) :: delta_V_electrode_capacitor
      real(dp) :: current

      call af_loop_box(tree, calculate_potential_from_box, .true.)

      ! The change in potential at the top electrode
      delta_V_electrode_dt = new_calculated_voltage_at_top_electrode &
                            - old_calculated_voltage_at_top_electrode

      ! Change the potential of the top electrode
      top_electrode_voltage = top_electrode_voltage + delta_V_electrode_dt

      print *, "Old calculated potential at top electrode: ", &
                      old_calculated_voltage_at_top_electrode
      print *, "New calculated potential at top electrode: ", &
                      new_calculated_voltage_at_top_electrode
      print *, "Delta potential at top electrode (N - O): ", &
                      delta_V_electrode_dt

      ! There is now a new potential difference between capacitor and the top electrode
      delta_V_electrode_capacitor = capacitor_voltage - top_electrode_voltage

      ! This potential difference creates a current between the top electrode and capacitor
      current = delta_V_electrode_capacitor / resistor

      ! This current will change the potential of the capacitor
      capacitor_voltage = capacitor_voltage - (current * dt) / capacitor_capacitance

      ! This current will also charge the top electrode again
      top_electrode_voltage = top_electrode_voltage + (current * dt) / top_electrode_self_capacitance
      
      ! Change the applied voltage at the top electrode for the simulation
      field_amplitude = -top_electrode_voltage / ST_domain_len(NDIM)

      old_calculated_voltage_at_top_electrode = new_calculated_voltage_at_top_electrode
      new_calculated_voltage_at_top_electrode = 0.0


   end subroutine add_circuit_effect

   subroutine calculate_potential_from_box(box)
      use m_units_constants
      type(box_t), intent(inout)     :: box
      real(dp) :: box_dr(NDIM)
      real(dp) :: box_rmin(NDIM)
      real(dp) :: cell_r_center(NDIM)
      real(dp) :: cell_volume
      integer :: n_cells
      integer :: i_cell, j_cell
      integer :: i_rhs_cell
      real(dp) :: rhs_cell
      real(dp) :: total_charge_cell
      real(dp) :: distance_cell_to_electrode
      real(dp) :: potential_from_cell
      
      n_cells = box%n_cell
      box_rmin = box%r_min
      box_dr = box%dr

      ! Loop through each cell in the box
      do i_cell = 1, n_cells
         do j_cell = 1, n_cells
            cell_r_center(1) = box_rmin(1) + i_cell * box_dr(1)
            cell_r_center(2) = box_rmin(2) + j_cell * box_dr(2)
            ! Calculate the cell volume (Cylindrical cell volume = 2 * pi * r * dr * dz)
            cell_volume = 2 * UC_pi * (cell_r_center(1) - 0.5 * box_dr(1)) * box_dr(1) * box_dr(2)
            ! Retrieve the charge density in this cell
            ! rhs = -rho / epsilon0
            i_rhs_cell = af_find_cc_variable(tree, "rhs")
# if NDIM == 2
            rhs_cell = box%cc(i_cell, j_cell, i_rhs_cell)
# endif
            ! Calculate the total charge = charge density * cell volume
            total_charge_cell = -rhs_cell * UC_eps0 * cell_volume
            ! Calculate the potential of this cell to the electrode position
            distance_cell_to_electrode = sqrt((cell_r_center(1) - r_top_electrode(1))**2 + &
                                              (cell_r_center(2) - r_top_electrode(2))**2)
            potential_from_cell = total_charge_cell / (4 * UC_pi * distance_cell_to_electrode)
            ! Increment potential at top electrode with the calculated potential of this cell
            new_calculated_voltage_at_top_electrode = new_calculated_voltage_at_top_electrode + potential_from_cell
            !print *, "Potential from cell: ", potential_from_cell
         end do
      end do

   end subroutine calculate_potential_from_box

end program streamer
