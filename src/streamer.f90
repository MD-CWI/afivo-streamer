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
  use m_units_constants

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
  ! Circuit stuff. fine lets do it with global variables
  logical :: use_circuit_andy = .false.
  integer :: circuit_andy_output_unit = 50
  real(dp), dimension(NDIM) :: r_top_electrode
  real(dp) :: old_calculated_voltage_at_top_electrode = 0.0
  real(dp) :: new_calculated_voltage_at_top_electrode = 0.0
  real(dp) :: top_electrode_voltage = 0.0
  real(dp) :: top_electrode_initial_voltage_fraction = 1
  real(dp) :: top_electrode_self_capacitance = 0.0
  real(dp) :: resistor = 300
  real(dp) :: capacitor_initial_voltage_fraction = 1
  real(dp) :: capacitor_capacitance = 200e-12
  real(dp) :: capacitor_voltage = 0.0
  real(dp), dimension(NDIM) :: r_bottom_electrode
  real(dp) :: bottom_electrode_voltage = 0.0
  real(dp) :: gnd_resistor = 10
  real(dp) :: old_calculated_voltage_at_bottom_electrode = 0.0
  real(dp) :: new_calculated_voltage_at_bottom_electrode = 0.0
  real(dp) :: delta_electrode_position = 1e-6  ! We will put the calculated potential position a bit out of the grid to avoid 1 / r influence
  logical :: use_circuit_R = .false.
  integer :: circuit_R_output_unit = 55
  real(dp) :: Q0_R = UC_elem_charge
  integer :: ix_weighting_potential_top
  integer :: ix_weighting_potential_bottom

 ! Add a variable to store the weighting potential from the top electrode in
  call af_add_cc_variable(tree, "weighting_pot_top")
  ix_weighting_potential_top = af_find_cc_variable(tree, "weighting_pot_top")
  call af_set_cc_methods(tree, ix_weighting_potential_top, af_bc_neumann_zero)

! Add a variable to store the weighting potential from the bottom electrode in
  call af_add_cc_variable(tree, "weighting_pot_bottom")
  ix_weighting_potential_bottom = af_find_cc_variable(tree, "weighting_pot_bottom")
  call af_set_cc_methods(tree, ix_weighting_potential_bottom, af_bc_neumann_zero)


  call print_program_name()

  ! Parse command line configuration files and options
  call CFG_update_from_arguments(cfg)

  call CFG_add_get(cfg, "use_circuit_andy", use_circuit_andy, &
         "If set use the circuit implementation of Andy.")
   call CFG_add_get(cfg, "use_circuit_R", use_circuit_R, &
         "If set use the circuit implementation of W. Riegler.")
  call CFG_add_get(cfg, "capacitor_capacitance", capacitor_capacitance, &
         "Capacitance of the capacitor used in the RC circuit implementation of Andy.")
  call CFG_add_get(cfg, "resistor", resistor, &
         "Resistor value for the resistor in the RC circuit implementation of Andy.")
  call CFG_add_get(cfg, "gnd_resistor", gnd_resistor, &
         "Resistor value for the connection of the bottom plate to ground.")
  call CFG_add_get(cfg, "capacitor_initial_voltage_fraction", capacitor_initial_voltage_fraction, &
         "Initial capacitor voltage as a fraction of the applied electric field for RC Andy.")
   call CFG_add_get(cfg, "top_electrode_initial_voltage_fraction", top_electrode_initial_voltage_fraction, &
        "Initial top electrode voltage as a fraction of the applied electric field for RC Andy.")

  call CFG_add_get(cfg, "restart_from_file", restart_from_file, &
       "If set, restart simulation from a previous .dat file")
  call CFG_add_get(cfg, "memory_limit_GB", memory_limit_GB, &
       "Memory limit (GB)")

  call initialize_modules(cfg, tree, mg, restart_from_file /= undefined_str)

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

     evolve_electrons = .true.
     if (associated(user_evolve_electrons)) &
          evolve_electrons = user_evolve_electrons(tree, time)

     if (evolve_electrons) then
        if (photoi_enabled .and. mod(it, photoi_per_steps) == 0) then
           call photoi_set_src(tree, time - photoi_prev_time)
           photoi_prev_time = time
        end if

        call af_advance(tree, dt, dt_lim, time, &
             species_itree(n_gas_species+1:n_species), &
             af_heuns_method, forward_euler)

#if NDIM >= 2
      if (use_circuit_andy) then
         call add_circuit_effect(tree)
      end if

      if (use_circuit_R) then
         call add_circuit_effect_R(tree)
      end if
#endif
        ! Make sure field is available for latest time state
        call field_compute(tree, mg, 0, time, .true.)

        !print *, "Top electrode voltage: ", top_electrode_voltage
        !print *, "Field amplitude: ", field_get_amplitude(tree, time)

        if (gas_dynamics) call coupling_add_fluid_source(tree, dt)
        if (circuit_used) call circuit_update(tree, dt)
     else
        dt_lim = dt_max
     end if

     if (gas_dynamics) then
        ! Go back to time at beginning of step
        time = global_time

     ! Make sure field is available for latest time state
        call af_advance(tree, dt, dt_gas_lim, time, &
             gas_vars, time_integrator, gas_forward_euler)
        call coupling_update_gas_density(tree)
     else
        dt_gas_lim = dt_max
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

   if (use_circuit_andy) then
      close(circuit_andy_output_unit)
   end if

   if (use_circuit_R) then
      close(circuit_R_output_unit)
   end if

   if (circuit_used) then
      close(circuit_jannis_output_unit)
   endif

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


    if (use_circuit_andy) then
      call initialize_circuit_andy(tree, cfg)
    end if

    if (use_circuit_R) then
      call initialize_circuit_R(tree, cfg)
    end if
    
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


   subroutine initialize_circuit_andy(tree, cfg)
      use m_field
      type(CFG_t), intent(inout) :: cfg
      type(af_t), intent(inout)  :: tree

      ! Initialize capacitor voltage to be the same as applied voltage
      capacitor_voltage = capacitor_initial_voltage_fraction * (-field_get_amplitude(tree, 0.0_dp) * ST_domain_len(NDIM))
      print *, "Initial capacitor voltage: ", capacitor_voltage
      r_top_electrode(1) = 0
      r_top_electrode(2) = 0
      r_top_electrode(NDIM) = ST_domain_len(NDIM) + delta_electrode_position
      top_electrode_voltage = top_electrode_initial_voltage_fraction * (-field_get_amplitude(tree, 0.0_dp) * ST_domain_len(NDIM))
      call field_set_voltage_externally(top_electrode_voltage)
      ! Formula for self capacitance of a infinitely flat disc with radius R; C = 8 * eps0 * R
      top_electrode_self_capacitance = 8 * UC_eps0 * ST_domain_len(1)

      r_bottom_electrode(1) = 0
      r_bottom_electrode(2) = 0
      r_bottom_electrode(NDIM) = 0 - delta_electrode_position

      ! Create an output file
      open(newunit=circuit_andy_output_unit, file=trim(output_name) // "_" // "output", action='write')
      write(circuit_andy_output_unit, *) "dt(s)           I(A)           Ve(V)           Vc(V)           Vgnd(V)           Ignd(A)"


   end subroutine initialize_circuit_andy

   subroutine add_circuit_effect(tree)
      type(af_t), intent(inout)     :: tree
      real(dp) :: delta_V_top_electrode_dt
      real(dp) :: delta_V_electrode_capacitor
      real(dp) :: current_top
      real(dp) :: delta_V_bottom_electrode_dt
      real(dp) :: current_bottom

      call af_loop_box(tree, calculate_potential_from_box, .true.)

      !!!! TOP ELECTRODE !!!!
      ! The change in potential at the top electrode
      delta_V_top_electrode_dt = new_calculated_voltage_at_top_electrode &
                            - old_calculated_voltage_at_top_electrode
      print *, "Change in potential top electrode: ", delta_V_top_electrode_dt


      ! Change the potential of the top electrode
      top_electrode_voltage = top_electrode_voltage + delta_V_top_electrode_dt

      ! There is now a new potential difference between capacitor and the top electrode
      delta_V_electrode_capacitor = capacitor_voltage - top_electrode_voltage

      ! This potential difference creates a current from capacitor to top electrode (conventional)
      current_top = delta_V_electrode_capacitor / resistor

      ! This current will change the potential of the capacitor
      capacitor_voltage = capacitor_voltage - (current_top * dt) / capacitor_capacitance

      ! This current will also charge the top electrode again
      top_electrode_voltage = top_electrode_voltage + (current_top * dt) / top_electrode_self_capacitance
      
      !!!! BOTTOM ELECTRODE !!!!
      delta_V_bottom_electrode_dt = new_calculated_voltage_at_bottom_electrode - old_calculated_voltage_at_bottom_electrode

      bottom_electrode_voltage = bottom_electrode_voltage + delta_V_bottom_electrode_dt

      ! Current from bottom electrode to ground (conventional)
      current_bottom = (bottom_electrode_voltage - 0) / gnd_resistor

      bottom_electrode_voltage = bottom_electrode_voltage - (current_bottom * dt) / top_electrode_self_capacitance

      ! Change the applied voltage at the top electrode for the simulation
      ! For simulations we still have bottom electrode  = 0 V and top electrode has the applied potential
      call field_set_voltage_externally(top_electrode_voltage - bottom_electrode_voltage)

      old_calculated_voltage_at_top_electrode = new_calculated_voltage_at_top_electrode
      new_calculated_voltage_at_top_electrode = 0.0

      old_calculated_voltage_at_bottom_electrode = new_calculated_voltage_at_bottom_electrode
      new_calculated_voltage_at_bottom_electrode = 0.0

      write(circuit_andy_output_unit, *) dt, current_top, top_electrode_voltage, capacitor_voltage, &
                                         bottom_electrode_voltage, current_bottom

   end subroutine add_circuit_effect

   subroutine calculate_potential_from_box(box)
      use m_units_constants
      use m_chemistry
      type(box_t), intent(inout)     :: box
      real(dp) :: box_dr(NDIM)
      real(dp) :: box_rmin(NDIM)
      real(dp) :: cell_r_center(NDIM)
      real(dp) :: cell_volume
      integer :: n_cells
      integer :: i_cell, j_cell
      integer :: cell_idx(NDIM)
      real(dp) :: charge_density_cell
      integer :: i_plasma_species
      real(dp) :: total_charge_cell
      real(dp) :: distance_cell_to_top_electrode
      real(dp) :: distance_cell_to_bottom_electrode
      real(dp) :: potential_from_cell_top
      real(dp) :: potential_from_cell_bottom
#if NDIM == 3
      integer :: k_cell
#endif
      n_cells = box%n_cell
      box_rmin = box%r_min
      box_dr = box%dr

#if NDIM == 2
      do i_cell = 1, n_cells
         do j_cell = 1, n_cells
            charge_density_cell = 0.0

            ! Calculate the center of each cell in the box
            cell_idx(1) = i_cell
            cell_idx(2) = j_cell
            cell_r_center = af_r_cc(box, cell_idx)

            ! Calculate the cell volume (Cylindrical cell volume = 2 * pi * r * dr * dz)
            cell_volume = 2 * UC_pi * (cell_r_center(1) - 0.5 * box_dr(1)) * box_dr(1) * box_dr(2)
            ! Retrieve the charge density in this cell
            do i_plasma_species = 1, size(charged_species_charge)
               charge_density_cell = charge_density_cell + &
                                 charged_species_charge(i_plasma_species) * &
                                 box%cc(i_cell, j_cell, charged_species_itree(i_plasma_species))
            end do
            charge_density_cell = charge_density_cell * UC_elem_charge

            total_charge_cell = charge_density_cell * cell_volume

            ! Calculate the potential of this cell to the top electrode position
            distance_cell_to_top_electrode = sqrt((cell_r_center(1) - r_top_electrode(1))**2 + &
                                                  (cell_r_center(2) - r_top_electrode(2))**2)
            potential_from_cell_top = total_charge_cell / (4 * UC_pi * UC_eps0 * distance_cell_to_top_electrode)
           
            ! Increment potential at top electrode with the calculated potential of this cell
            !$omp critical
            new_calculated_voltage_at_top_electrode = new_calculated_voltage_at_top_electrode + potential_from_cell_top
            !$omp end critical

            ! Calculate the potential of this cell to the bottom electrode position
            distance_cell_to_bottom_electrode = sqrt((cell_r_center(1) - r_bottom_electrode(1))**2 + &
                                                     (cell_r_center(2) - r_bottom_electrode(2))**2)
            potential_from_cell_bottom = total_charge_cell / (4 * UC_pi * UC_eps0 * distance_cell_to_bottom_electrode)
            
            ! Increment potential at top electrode with the calculated potential of this cell
            !$omp critical
            new_calculated_voltage_at_bottom_electrode = new_calculated_voltage_at_bottom_electrode + potential_from_cell_bottom
            !$omp end critical

         end do
      end do
#endif

#if NDIM == 3
do i_cell = 1, n_cells
   do j_cell = 1, n_cells
      do k_cell = 1, n_cells
         charge_density_cell = 0

         ! Calculate the center of each cell in the box
         cell_idx(1) = i_cell
         cell_idx(2) = j_cell
         cell_idx(3) = k_cell
         cell_r_center = af_r_cc(box, cell_idx)
         
         ! Calculate the cell volume
         cell_volume = box_dr(1) * box_dr(2) * box_dr(3)
         ! Retrieve the charge density in this cell
         do i_plasma_species = 1, size(charged_species_charge)
            charge_density_cell = charge_density_cell + &
                              charged_species_charge(i_plasma_species) * &
                              box%cc(i_cell, j_cell, k_cell, charged_species_itree(i_plasma_species))
         end do
         charge_density_cell = charge_density_cell * UC_elem_charge

         total_charge_cell = charge_density_cell * cell_volume

         ! Calculate the potential of this cell to the electrode position
         distance_cell_to_top_electrode = sqrt((cell_r_center(1) - r_top_electrode(1))**2 + &
                                               (cell_r_center(2) - r_top_electrode(2))**2 + &
                                               (cell_r_center(3) - r_top_electrode(3))**2)
         potential_from_cell_top = total_charge_cell / (4 * UC_pi * UC_eps0 * distance_cell_to_top_electrode)
         
         ! Increment potential at top electrode with the calculated potential of this cell
         !$omp critical
         new_calculated_voltage_at_top_electrode = new_calculated_voltage_at_top_electrode + potential_from_cell_top
         !$omp end critical

          ! Calculate the potential of this cell to the electrode position
         distance_cell_to_bottom_electrode = sqrt((cell_r_center(1) - r_bottom_electrode(1))**2 + &
                                                  (cell_r_center(2) - r_bottom_electrode(2))**2 + &
                                                  (cell_r_center(3) - r_bottom_electrode(3))**2)
         potential_from_cell_bottom = total_charge_cell / (4 * UC_pi * UC_eps0 * distance_cell_to_bottom_electrode)
        
         ! Increment potential at top electrode with the calculated potential of this cell
         !$omp critical
         new_calculated_voltage_at_bottom_electrode = new_calculated_voltage_at_bottom_electrode + potential_from_cell_bottom
         !$omp end critical

      end do
   end do
end do
#endif
end subroutine calculate_potential_from_box


subroutine initialize_circuit_R(tree, cfg)
   use m_field
   type(CFG_t), intent(inout) :: cfg
   type(af_t), intent(inout)  :: tree
   real(dp) :: field_voltage_temp
   real(dp) :: max_potential = 0.0
   real(dp) :: min_potential = 0.0
   integer :: i_fas_fmg

   ! Initialize capacitor voltage to be the a fraction of the applied voltage
   capacitor_voltage = capacitor_initial_voltage_fraction * (-field_get_amplitude(tree, 0.0_dp) * ST_domain_len(NDIM))
   print *, "Initial capacitor voltage: ", capacitor_voltage
   top_electrode_voltage = top_electrode_initial_voltage_fraction * (-field_get_amplitude(tree, 0.0_dp) * ST_domain_len(NDIM))
   call field_set_voltage_externally(top_electrode_voltage)
   ! Formula for self capacitance of a infinitely flat disc with radius R; C = 8 * eps0 * R
   top_electrode_self_capacitance = 8 * UC_eps0 * ST_domain_len(1)

   !!!!! Calculate weighting potential at t = 0 !!!!!
   ! TOP ELECTRODE
   ! Set RHS to 0
   call af_tree_clear_cc(tree, i_rhs)
   ! Store field voltage in temp
   field_voltage_temp = field_voltage
   ! Set field voltage to desired value
   call field_set_voltage_externally(Q0_R / top_electrode_self_capacitance)
   ! Calculate the weighting potential field using FAS FMG
   do i_fas_fmg = 1, 100
      call mg_fas_fmg(tree, mg, .false., .false.)
   end do
   ! copy cc var i_phi to i_phiR_top
   call af_tree_copy_cc(tree, i_phi, ix_weighting_potential_top)

   ! BOTTOM ELECTRODE
   ! Set RHS to 0
   call af_tree_clear_cc(tree, i_rhs)
   ! Set phi to 0
   call af_tree_clear_cc(tree, i_phi)
   ! Use specially defined boundary conditions for bottom electrode
   mg%sides_bc => weighting_field_bottom_bc
   ! Calculate the weighting potential field using FAS FMG
   do i_fas_fmg = 1, 100
      call mg_fas_fmg(tree, mg, .false., .false.)
   end do
   ! copy cc var i_phi to the weighting_potential_bottom
   call af_tree_copy_cc(tree, i_phi, ix_weighting_potential_bottom)
   ! Set normal boundary conditions back
   mg%sides_bc => field_bc_homogeneous
   ! set field voltage back from temp
   call field_set_voltage_externally(field_voltage_temp)
   ! field compute with rhs to obtain original field again
   call field_compute(tree, mg, 0, 0.0_dp, .false.)
   ! Calculate the initial contribution of the charges in the gap to the induced potential on the top electrode
   ! integral psi_V * charge_density  dV / Q0
   call af_loop_box(tree, calculate_Rterm_both_electrodes_from_box, .true.)
   ! We dont use the induced potential at the initial timestep. We need a difference in times (moving charges)
   old_calculated_voltage_at_top_electrode = new_calculated_voltage_at_top_electrode
   new_calculated_voltage_at_top_electrode = 0.0

   old_calculated_voltage_at_bottom_electrode = new_calculated_voltage_at_bottom_electrode
   new_calculated_voltage_at_bottom_electrode = 0.0
 
   call af_tree_max_cc(tree, ix_weighting_potential_bottom, max_potential)
   print *, "Max potential: ", max_potential
   print *, "Real max potential: ", Q0_R / top_electrode_self_capacitance
   call af_tree_min_cc(tree, ix_weighting_potential_bottom, min_potential)
   print *, "Min potential: ", min_potential
   print *, "Real min potential: ", 0

   ! Create an output file
   open(newunit=circuit_R_output_unit, file=trim(output_name) // "_" // "output", action='write')
   write(circuit_R_output_unit, *) "dt(s)           I(A)           Ve(V)           Vc(V)           Vgnd(V)           Ignd(A)"

end subroutine initialize_circuit_R


subroutine add_circuit_effect_R(tree)
   type(af_t), intent(inout)     :: tree
   real(dp) :: delta_V_top_electrode_dt
   real(dp) :: delta_V_electrode_capacitor
   real(dp) :: current_top
   real(dp) :: delta_V_bottom_electrode_dt
   real(dp) :: current_bottom

   ! Calculate the contribution of the charges in the gap to the induced potential on the top electrode
   ! integral psi_V * charge_density  dV / Q0
   call af_loop_box(tree, calculate_Rterm_both_electrodes_from_box, .true.)

   !!!! TOP ELECTRODE !!!!
   ! The change in potential at the top electrode
   delta_V_top_electrode_dt = new_calculated_voltage_at_top_electrode &
                         - old_calculated_voltage_at_top_electrode
   !print *, "Change in potential top electrode: ", delta_V_top_electrode_dt

   ! Change the potential of the top electrode
   top_electrode_voltage = top_electrode_voltage + delta_V_top_electrode_dt

   ! There is now a new potential difference between capacitor and the top electrode
   delta_V_electrode_capacitor = capacitor_voltage - top_electrode_voltage

   ! This potential difference creates a current from capacitor to top electrode (conventional)
   current_top = delta_V_electrode_capacitor / resistor

   ! This current will change the potential of the capacitor
   capacitor_voltage = capacitor_voltage - (current_top * dt) / capacitor_capacitance

   ! This current will also charge the top electrode again
   top_electrode_voltage = top_electrode_voltage + (current_top * dt) / top_electrode_self_capacitance
   
   !!!! BOTTOM ELECTRODE !!!!
   delta_V_bottom_electrode_dt = new_calculated_voltage_at_bottom_electrode - old_calculated_voltage_at_bottom_electrode
   !print *, "Change in potential bottom electrode: ", delta_V_bottom_electrode_dt

   bottom_electrode_voltage = bottom_electrode_voltage + delta_V_bottom_electrode_dt

   ! Current from bottom electrode to ground (conventional)
   current_bottom = (bottom_electrode_voltage - 0) / gnd_resistor

   bottom_electrode_voltage = bottom_electrode_voltage - (current_bottom * dt) / top_electrode_self_capacitance

   ! Change the applied voltage at the top electrode for the simulation
   ! For simulations we still have bottom electrode  = 0 V and top electrode has the applied potential
   call field_set_voltage_externally(top_electrode_voltage - bottom_electrode_voltage)

   old_calculated_voltage_at_top_electrode = new_calculated_voltage_at_top_electrode
   new_calculated_voltage_at_top_electrode = 0.0

   old_calculated_voltage_at_bottom_electrode = new_calculated_voltage_at_bottom_electrode
   new_calculated_voltage_at_bottom_electrode = 0.0

   write(circuit_R_output_unit, *) dt, current_top, top_electrode_voltage, capacitor_voltage, &
                                      bottom_electrode_voltage, current_bottom

end subroutine add_circuit_effect_R

subroutine calculate_Rterm_both_electrodes_from_box(box)
   use m_units_constants
   use m_chemistry
   type(box_t), intent(inout)     :: box
   real(dp) :: box_dr(NDIM)
   real(dp) :: box_rmin(NDIM)
   real(dp) :: cell_r_center(NDIM)
   real(dp) :: cell_volume
   integer :: n_cells
   integer :: i_cell, j_cell
   integer :: cell_idx(NDIM)
   real(dp) :: charge_density_cell
   integer :: i_plasma_species
   real(dp) :: total_charge_cell
   real(dp) :: potential_from_cell_top
   real(dp) :: potential_from_cell_bottom
#if NDIM == 3
   integer :: k_cell
#endif
   n_cells = box%n_cell
   box_rmin = box%r_min
   box_dr = box%dr

#if NDIM == 2
   do i_cell = 1, n_cells
      do j_cell = 1, n_cells
         charge_density_cell = 0.0

         ! Calculate the center of each cell in the box
         cell_idx(1) = i_cell
         cell_idx(2) = j_cell
         cell_r_center = af_r_cc(box, cell_idx)

         ! Calculate the cell volume (Cylindrical cell volume = 2 * pi * r * dr * dz)
         cell_volume = 2 * UC_pi * (cell_r_center(1) - 0.5 * box_dr(1)) * box_dr(1) * box_dr(2)
         ! Retrieve the charge density in this cell
         do i_plasma_species = 1, size(charged_species_charge)
            charge_density_cell = charge_density_cell + &
                              charged_species_charge(i_plasma_species) * &
                              box%cc(i_cell, j_cell, charged_species_itree(i_plasma_species))
         end do
         charge_density_cell = charge_density_cell * UC_elem_charge

         total_charge_cell = charge_density_cell * cell_volume

         ! Multiply the total charge of this cell with the weighting potential to get the induced potential
         potential_from_cell_top = total_charge_cell * box%cc(i_cell, j_cell, ix_weighting_potential_top) / Q0_R
         ! Increment potential at top electrode with the calculated potential of this cell
         !$omp critical
         new_calculated_voltage_at_top_electrode = new_calculated_voltage_at_top_electrode + potential_from_cell_top
         !$omp end critical

         ! Multiply the total charge of this cell with the weighting potential to get the induced potential
         potential_from_cell_bottom = total_charge_cell * box%cc(i_cell, j_cell, ix_weighting_potential_bottom) / Q0_R
         ! Increment potential at bottom electrode with the calculated potential of this cell
         !$omp critical
         new_calculated_voltage_at_bottom_electrode = new_calculated_voltage_at_bottom_electrode + potential_from_cell_bottom
         !$omp end critical
      end do
   end do
#endif

#if NDIM == 3
do i_cell = 1, n_cells
do j_cell = 1, n_cells
   do k_cell = 1, n_cells
      charge_density_cell = 0

      ! Calculate the center of each cell in the box
      cell_idx(1) = i_cell
      cell_idx(2) = j_cell
      cell_idx(3) = k_cell
      cell_r_center = af_r_cc(box, cell_idx)
      
      ! Calculate the cell volume
      cell_volume = box_dr(1) * box_dr(2) * box_dr(3)
      ! Retrieve the charge density in this cell
      do i_plasma_species = 1, size(charged_species_charge)
         charge_density_cell = charge_density_cell + &
                           charged_species_charge(i_plasma_species) * &
                           box%cc(i_cell, j_cell, k_cell, charged_species_itree(i_plasma_species))
      end do
      charge_density_cell = charge_density_cell * UC_elem_charge

      total_charge_cell = charge_density_cell * cell_volume

      ! Multiply the total charge of this cell with the weighting potential to get the induced potential
      potential_from_cell_top = total_charge_cell * box%cc(i_cell, j_cell, k_cell, ix_weighting_potential_top) / Q0_R
      
      ! Increment potential at top electrode with the calculated potential of this cell
      !$omp critical
      new_calculated_voltage_at_top_electrode = new_calculated_voltage_at_top_electrode + potential_from_cell_top
      !$omp end critical

      ! Multiply the total charge of this cell with the weighting potential to get the induced potential
      potential_from_cell_bottom = total_charge_cell * box%cc(i_cell, j_cell, k_cell, ix_weighting_potential_bottom) / Q0_R
      ! Increment potential at bottom electrode with the calculated potential of this cell
      !$omp critical
      new_calculated_voltage_at_bottom_electrode = new_calculated_voltage_at_bottom_electrode + potential_from_cell_bottom
      !$omp end critical
   end do
end do
end do
#endif
end subroutine calculate_Rterm_both_electrodes_from_box


!> This fills ghost cells near physical boundaries for the potential
subroutine weighting_field_bottom_bc(box, nb, iv, coords, bc_val, bc_type)
   type(box_t), intent(in) :: box
   integer, intent(in)     :: nb
   integer, intent(in)     :: iv
   real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
   real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
   integer, intent(out)    :: bc_type

   if (af_neighb_dim(nb) == NDIM) then
      if (af_neighb_low(nb)) then
         bc_type = af_bc_dirichlet
         bc_val = Q0_R / top_electrode_self_capacitance
      else
         bc_type = af_bc_dirichlet
         bc_val  = 0.0
      end if
   else
      bc_type = af_bc_neumann
      bc_val = 0.0_dp
   end if
 end subroutine weighting_field_bottom_bc

end program streamer
