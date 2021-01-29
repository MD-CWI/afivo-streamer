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
  real(dp) :: Q0_R = UC_elem_charge
  integer :: ix_weighting_potential_top
  integer :: ix_weighting_potential_bottom

  integer  :: n_streamers_circuit = 1
  character(len=string_len) :: circuit_method = "none"  ! Can be andy, riegler, pederson
  character(len=string_len) :: circuit_type = "slc" ! can be slc (small lab circuit) or marx (marx generator)
  integer :: circuit_output_unit = 55

   ! Small Lab Circuit:
   !                  (Timed switch, voltage pulse time)
   !  |----------------|------ /---R_postswitch---|
   !  |                |      |                   |
   !  |                |      |                L_postswitch
   !  |                |      |                   |
   !  |                |      |                   |
   ! R_preswitch   C_preswitch|                  \_/ (Electrode, has C_top)
   !  |                |      |            
   !  |                |      |            
   !  HV_source        |      |               ---------- (Electrode, has C_bottom)
   !  |                |      |                   |
   !  |                |      |                 R_gnd
   !  |                |      |                   |
   !  |                |      |                   |
   !  |                |      |                  L_gnd 
   !  |                |      |                   |
   !  |------GND-------|------|-------------------|
                                        
  real(dp) :: HV_source = 0.0  ! Constant Voltage Source
  real(dp) :: HV_power = 30  ! Power limit of voltage source
  real(dp) :: T_switch_on = 10e-9  ! Time the switch connects the electrode to the HV source
  real(dp) :: T_switch_off = 20e-9  ! Time the switch connects the electrode to ground
  real(dp) :: slc_pulse_period = 0  ! Complete pulse time = T_switch_on + T_switch_off
  real(dp) :: R_preswitch = 10e3  ! Resistor for charging capacitor
  real(dp) :: R_postswitch = 300  ! Current limiting resistor for electrode
  real(dp) :: R_gnd = 10  ! Resistance between "ground" electrode and actual GND

  real(dp) :: C_preswitch = 200e-12  ! Main capacitor
  real(dp) :: C_preswitch_initial_voltage_fraction = 1
  real(dp) :: C_preswitch_voltage = 0  ! Stores V(t + 1)
  real(dp) :: C_preswitch_voltage_old = 0  ! Used for when V(t) is needed.
  real(dp) :: C_preswitch_voltage_oldold = 0  ! Used for when V(t - 1) is needed.

  real(dp) :: L_postswitch = 0  ! Self inductance of wire connecting electrode
  real(dp) :: L_gnd = 0  ! Self inductance of wire connecting bottom electrode to GND

  real(dp) :: C_top = 0  ! Self inductance of top electrode
  real(dp) :: top_electrode_voltage = 0.0  ! Stores V(t + 1)
  real(dp) :: top_electrode_voltage_old = 0.0  ! Used for when V(t) is needed
  real(dp) :: top_electrode_voltage_oldold = 0.0  ! Used for when V(t - 1) is needed
  real(dp) :: top_electrode_initial_voltage_fraction = 1
  real(dp) :: top_electrode_old_calculated_voltage = 0.0
  real(dp) :: top_electrode_new_calculated_voltage = 0.0
  logical :: use_top_electrode_sphere = .false.

  real(dp) :: C_bottom = 0  ! Self inductance of bottom electrode
  real(dp) :: bottom_electrode_voltage = 0.0
  real(dp) :: bottom_electrode_old_calculated_voltage = 0.0
  real(dp) :: bottom_electrode_new_calculated_voltage = 0.0

  ! Circuit andy specific stuff
  real(dp), dimension(NDIM) :: r_bottom_electrode
  real(dp) :: delta_electrode_position = 1e-6  ! We will put the calculated potential position a bit out of the grid to avoid 1 / r influence
  real(dp), dimension(NDIM) :: r_top_electrode

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

   call CFG_add_get(cfg, "n_streamers_circuit", n_streamers_circuit, &
   "Number of streamers we pretend to simulate for power consumption of the circuit.")
   call CFG_add_get(cfg, "circuit_method", circuit_method, &
   "Method used to add circuit effect (andy, riegler, or pederson")
   call CFG_add_get(cfg, "circuit_type", circuit_type, &
   "Type of circuit (slc (small lab circuit), or marx (marx generator)")
   call CFG_add_get(cfg, "HV_source", HV_source, &
   "(slc) Fixed applied HV. (V)")
   call CFG_add_get(cfg, "HV_power", HV_power, &
   "(slc) Power limit of HV supply. (W)")
   call CFG_add_get(cfg, "T_switch_on", T_switch_on, &
   "(slc) Time the switch connects the electrode to the HV source. (s)")
   call CFG_add_get(cfg, "T_switch_off", T_switch_off, &
   "(slc) Time the switch connects the electrode to ground. (s)")
   call CFG_add_get(cfg, "R_preswitch", R_preswitch, &
   "(slc) Resistor connected to HV supply. (Ohm)")
   call CFG_add_get(cfg, "R_postswitch", R_postswitch, &
   "(slc) Resistor connected to top powered electrode after the switch. (Ohm)")
   call CFG_add_get(cfg, "R_gnd", R_gnd, &
   "(slc) Resistance of the connection between bottom grounded electrode and real GND. (Ohm)")
   call CFG_add_get(cfg, "C_preswitch", C_preswitch, &
   "(slc) Main capacitor. Connected before the switch. (F)")
   call CFG_add_get(cfg, "C_preswitch_initial_voltage_fraction", C_preswitch_initial_voltage_fraction, &
   "(slc) Fraction of the applied HV that the C_preswitch starts with.")
   call CFG_add_get(cfg, "L_postswitch", L_postswitch, &
   "(slc) Inductance connected to top powered electrode after the switch. (H)")
   call CFG_add_get(cfg, "L_gnd", L_gnd, &
   "(slc) Inductance of the connection between bottom grounded electrode and real GND. (H)")
   call CFG_add_get(cfg, "top_electrode_initial_voltage_fraction", top_electrode_initial_voltage_fraction, &
   "Fraction of the applied HV that the top powered electrode starts with.")
   call CFG_add_get(cfg, "use_top_electrode_sphere", use_top_electrode_sphere, &
   "Approximate top powered electrode as a sphere with radius the tip radius.")
   
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

#if NDIM >= 2
      if (circuit_method == "andy") then
         call calc_induced_potential_andy(tree)
      else if (circuit_method == "riegler") then
         call calc_induced_potential_riegler(tree)
      else if (circuit_method == "pederson") then
         call calc_induced_potential_pederson(tree)
      end if

      if (circuit_method /= "none") then
         if (circuit_type == "slc") then
            call add_slc_effect()
         end if
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

   if (circuit_method /= "none") then
      close(circuit_output_unit)
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


    if (circuit_method /= "none") then
      call initialize_electrode_capacitance()

      if (circuit_type == "slc") then
         call initialize_slc()
      end if

      ! Create an output file
      open(newunit=circuit_output_unit, file=trim(output_name) // "_" // "output", action='write')
      write(circuit_output_unit, *) "dt(s)\tI_HV(A)\tIe(A)\tIc(A)\tV_Cpre(V)\tV_Ctop(V)\tV_gnd(V)\tI_gnd(A)"

    endif

    if (circuit_method == "andy") then
      call initialize_circuit_andy(tree, cfg)
    else if (circuit_method == "riegler") then
      call initialize_circuit_riegler(tree, cfg)
    else if (circuit_method == "pederson") then
      call initialize_circuit_pederson(tree, cfg)
    endif
    
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

   subroutine initialize_electrode_capacitance()
      use m_field
      if (use_top_electrode_sphere .and. ST_use_electrode) then
         ! Formula for self capacitance of a sphere with radius R; C = 4 * pi * eps0 * R
         C_top = 4 * UC_pi * UC_eps0 * field_rod_radius
      else
         ! Formula for self capacitance of a infinitely flat disc with radius R; C = 8 * eps0 * R
         C_top = 8 * UC_eps0 * ST_domain_len(1)
      end if

      ! Formula for self capacitance of a infinitely flat disc with radius R; C = 8 * eps0 * R
      C_bottom = 8 * UC_eps0 * ST_domain_len(1)

   end subroutine initialize_electrode_capacitance


   subroutine initialize_slc()
      ! Initialize the starting voltage of the main capacitor and the electrode
      C_preswitch_voltage = C_preswitch_initial_voltage_fraction * HV_source
      print *, "Initial main capacitor voltage: ", C_preswitch_voltage

      top_electrode_voltage = top_electrode_initial_voltage_fraction * HV_source
      call field_set_voltage_externally(top_electrode_voltage)
      print *, "Initial electrode voltage: ", top_electrode_voltage

      slc_pulse_period = T_switch_on + T_switch_off

   end subroutine initialize_slc

   subroutine initialize_circuit_andy(tree, cfg)
      use m_field
      type(CFG_t), intent(inout) :: cfg
      type(af_t), intent(inout)  :: tree

      r_top_electrode(1) = 0
      r_top_electrode(2) = 0
      r_top_electrode(NDIM) = ST_domain_len(NDIM) + delta_electrode_position
    
      r_bottom_electrode(1) = 0
      r_bottom_electrode(2) = 0
      r_bottom_electrode(NDIM) = 0 - delta_electrode_position
   end subroutine initialize_circuit_andy

   subroutine calc_induced_potential_andy(tree)
      type(af_t), intent(inout)     :: tree

      call af_loop_box(tree, calc_induced_potential_from_box_andy, .true.)
   end subroutine calc_induced_potential_andy

   subroutine calc_induced_potential_from_box_andy(box)
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
            top_electrode_new_calculated_voltage = top_electrode_new_calculated_voltage + potential_from_cell_top
            !$omp end critical

            ! Calculate the potential of this cell to the bottom electrode position
            distance_cell_to_bottom_electrode = sqrt((cell_r_center(1) - r_bottom_electrode(1))**2 + &
                                                     (cell_r_center(2) - r_bottom_electrode(2))**2)
            potential_from_cell_bottom = total_charge_cell / (4 * UC_pi * UC_eps0 * distance_cell_to_bottom_electrode)
            
            ! Increment potential at top electrode with the calculated potential of this cell
            !$omp critical
            bottom_electrode_new_calculated_voltage = bottom_electrode_new_calculated_voltage + potential_from_cell_bottom
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
               top_electrode_new_calculated_voltage = top_electrode_new_calculated_voltage + potential_from_cell_top
               !$omp end critical

               ! Calculate the potential of this cell to the electrode position
               distance_cell_to_bottom_electrode = sqrt((cell_r_center(1) - r_bottom_electrode(1))**2 + &
                                                      (cell_r_center(2) - r_bottom_electrode(2))**2 + &
                                                      (cell_r_center(3) - r_bottom_electrode(3))**2)
               potential_from_cell_bottom = total_charge_cell / (4 * UC_pi * UC_eps0 * distance_cell_to_bottom_electrode)
            
               ! Increment potential at top electrode with the calculated potential of this cell
               !$omp critical
               bottom_electrode_new_calculated_voltage = bottom_electrode_new_calculated_voltage + potential_from_cell_bottom
               !$omp end critical

            end do
         end do
      end do
#endif
   end subroutine calc_induced_potential_from_box_andy


   subroutine initialize_circuit_riegler(tree, cfg)
      use m_field
      type(CFG_t), intent(inout) :: cfg
      type(af_t), intent(inout)  :: tree
      real(dp) :: field_voltage_temp
      real(dp) :: max_potential = 0.0
      real(dp) :: min_potential = 0.0
      integer :: i_fas_fmg

      !!!!! Calculate weighting potential at t = 0 !!!!!
      ! TOP ELECTRODE
      ! Set RHS to 0
      call af_tree_clear_cc(tree, i_rhs)
      ! Store field voltage in temp
      field_voltage_temp = field_voltage
      ! Set field voltage to desired value
      call field_set_voltage_externally(Q0_R / C_top)
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
      call af_loop_box(tree,  calc_induced_potential_from_box_riegler, .true.)
      ! We dont use the induced potential at the initial timestep. We need a difference in times (moving charges)
      top_electrode_old_calculated_voltage = top_electrode_new_calculated_voltage
      top_electrode_new_calculated_voltage = 0.0

      bottom_electrode_old_calculated_voltage = bottom_electrode_new_calculated_voltage
      bottom_electrode_new_calculated_voltage = 0.0
   
      call af_tree_max_cc(tree, ix_weighting_potential_bottom, max_potential)
      print *, "Max potential: ", max_potential
      print *, "Real max potential: ", Q0_R / C_top
      call af_tree_min_cc(tree, ix_weighting_potential_bottom, min_potential)
      print *, "Min potential: ", min_potential
      print *, "Real min potential: ", 0

   end subroutine initialize_circuit_riegler

   subroutine calc_induced_potential_riegler(tree)
      type(af_t), intent(inout)     :: tree
      call af_loop_box(tree, calc_induced_potential_from_box_riegler, .true.)
   end subroutine calc_induced_potential_riegler

   subroutine calc_induced_potential_from_box_riegler(box)
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

            total_charge_cell = charge_density_cell * cell_volume * n_streamers_circuit

            ! Multiply the total charge of this cell with the weighting potential to get the induced potential
            potential_from_cell_top = total_charge_cell * box%cc(i_cell, j_cell, ix_weighting_potential_top) / Q0_R
            ! Increment potential at top electrode with the calculated potential of this cell
            !$omp critical
            top_electrode_new_calculated_voltage = top_electrode_new_calculated_voltage + potential_from_cell_top
            !$omp end critical

            ! Multiply the total charge of this cell with the weighting potential to get the induced potential
            potential_from_cell_bottom = total_charge_cell * box%cc(i_cell, j_cell, ix_weighting_potential_bottom) / Q0_R
            ! Increment potential at bottom electrode with the calculated potential of this cell
            !$omp critical
            bottom_electrode_new_calculated_voltage = bottom_electrode_new_calculated_voltage + potential_from_cell_bottom
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

               total_charge_cell = charge_density_cell * cell_volume * n_streamers_circuit

               ! Multiply the total charge of this cell with the weighting potential to get the induced potential
               potential_from_cell_top = total_charge_cell * box%cc(i_cell, j_cell, k_cell, ix_weighting_potential_top) / Q0_R
               
               ! Increment potential at top electrode with the calculated potential of this cell
               !$omp critical
               top_electrode_new_calculated_voltage = top_electrode_new_calculated_voltage + potential_from_cell_top
               !$omp end critical

               ! Multiply the total charge of this cell with the weighting potential to get the induced potential
               potential_from_cell_bottom = total_charge_cell * box%cc(i_cell, j_cell, k_cell, ix_weighting_potential_bottom) / Q0_R
               ! Increment potential at bottom electrode with the calculated potential of this cell
               !$omp critical
               bottom_electrode_new_calculated_voltage = bottom_electrode_new_calculated_voltage + potential_from_cell_bottom
               !$omp end critical
            end do
         end do
      end do
#endif
   end subroutine calc_induced_potential_from_box_riegler


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
            bc_val = Q0_R / C_top
         else
            bc_type = af_bc_dirichlet
            bc_val  = 0.0
         end if
      else
         bc_type = af_bc_neumann
         bc_val = 0.0_dp
      end if
   end subroutine weighting_field_bottom_bc


   subroutine initialize_circuit_pederson(tree, cfg)
      use m_field
      type(CFG_t), intent(inout) :: cfg
      type(af_t), intent(inout)  :: tree
      real(dp) :: field_voltage_temp
      real(dp) :: max_potential = 0.0
      real(dp) :: min_potential = 0.0
      integer :: i_fas_fmg 
   
   end subroutine initialize_circuit_pederson

   subroutine calc_induced_potential_pederson(tree)
      type(af_t), intent(inout)     :: tree
      call af_loop_box(tree, calc_induced_potential_from_box_pederson, .true.)
   end subroutine calc_induced_potential_pederson

   subroutine calc_induced_potential_from_box_pederson(box)
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

   end subroutine calc_induced_potential_from_box_pederson

   subroutine add_slc_effect()
      real(dp) :: delta_V_top_electrode_dt
      real(dp) :: delta_V_electrode_capacitor
      real(dp) :: delta_V_bottom_electrode_dt
      real(dp) :: time_within_slc_pulse_period
      real(dp) :: I_HV
      real(dp) :: Ie
      real(dp) :: Ic
      real(dp) :: I_gnd

      C_preswitch_voltage_oldold = C_preswitch_voltage_old
      C_preswitch_voltage_old = C_preswitch_voltage

      top_electrode_voltage_oldold = top_electrode_voltage_old
      top_electrode_voltage_old = top_electrode_voltage

      time_within_slc_pulse_period = time - floor(time / slc_pulse_period) * slc_pulse_period

      !!!! TOP ELECTRODE !!!!
      ! The change in potential at the top electrode
      delta_V_top_electrode_dt = top_electrode_new_calculated_voltage &
                                 - top_electrode_old_calculated_voltage
      !print *, "Change in potential top electrode: ", delta_V_top_electrode_dt


      ! Change the potential of the top electrode
      top_electrode_voltage_old = top_electrode_voltage_old + delta_V_top_electrode_dt

      if (time_within_slc_pulse_period <= T_switch_on) then
         ! Electrode is connected via the switch to the main capacitor and HV source

         top_electrode_voltage = (1 / ((C_top * R_postswitch / dt) +  (L_postswitch * C_top / dt**2))) * &
         (&
            ((C_top * R_postswitch / dt) + (2 * L_postswitch * C_top / dt**2) -1) * top_electrode_voltage_old &
            - (L_postswitch * C_top / dt**2) * top_electrode_voltage_oldold &
            + C_preswitch_voltage_old &
         )

         C_preswitch_voltage = (dt / (C_preswitch * R_preswitch)) * &
         (&
            - (1 - C_preswitch * R_preswitch) * C_preswitch_voltage_old &
            - (C_top * R_preswitch / dt) * top_electrode_voltage &
            + (C_top * R_preswitch / dt) * top_electrode_voltage_old &
            + HV_source &
         )


         Ie = C_top * (top_electrode_voltage_old - top_electrode_voltage_oldold) / dt
         I_HV = (C_top * (top_electrode_voltage_old - top_electrode_voltage_oldold) / dt) &
               + (C_preswitch * (C_preswitch_voltage_old - C_preswitch_voltage_oldold) / dt)
         Ic = C_preswitch * (C_preswitch_voltage_old - C_preswitch_voltage_oldold) / dt
         ! There is now a new potential difference between capacitor and the top electrode
         !delta_V_electrode_capacitor = C_preswitch_voltage - top_electrode_voltage
      
         ! This potential difference creates a current from capacitor to top electrode (conventional)
         !current_top = delta_V_electrode_capacitor / R_postswitch

         ! This current will change the potential of the capacitor
         !C_preswitch_voltage = C_preswitch_voltage - (current_top * dt) / C_preswitch

         ! This current will also charge the top electrode again
         !top_electrode_voltage = top_electrode_voltage + (current_top * dt) / C_top
      

      else
         ! Electrode is connected via R_postswitch to ground

         ! C_preswitch is being charged by HV_source through R_preswitch

      end if


      !!!! BOTTOM ELECTRODE !!!!
      delta_V_bottom_electrode_dt = bottom_electrode_new_calculated_voltage - bottom_electrode_old_calculated_voltage

      bottom_electrode_voltage = bottom_electrode_voltage + delta_V_bottom_electrode_dt

      ! Current from bottom electrode to ground (conventional)
      I_gnd = (bottom_electrode_voltage - 0) / R_gnd

      bottom_electrode_voltage = bottom_electrode_voltage - (I_gnd * dt) / C_bottom

      ! Change the applied voltage at the top electrode for the simulation
      ! For simulations we still have bottom electrode  = 0 V and top electrode has the applied potential
      call field_set_voltage_externally(top_electrode_voltage - bottom_electrode_voltage)

      top_electrode_old_calculated_voltage = top_electrode_new_calculated_voltage
      top_electrode_new_calculated_voltage = 0.0

      bottom_electrode_old_calculated_voltage = bottom_electrode_new_calculated_voltage
      bottom_electrode_new_calculated_voltage = 0.0

      write(circuit_output_unit, *) dt, I_HV, Ie, Ic, C_preswitch_voltage, top_electrode_voltage, &
                                    bottom_electrode_voltage, I_gnd

   end subroutine add_slc_effect

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
