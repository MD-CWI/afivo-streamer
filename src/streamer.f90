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
  use m_advance
  use m_types
  use m_user_methods
  use m_output

  implicit none

  integer, parameter        :: int8       = selected_int_kind(18)
  integer(int8)             :: t_start, t_current, count_rate
  real(dp)                  :: wc_time, inv_count_rate
  real(dp)                  :: time_last_print, time_last_output
  integer                   :: i, s, it, coord_type, box_bytes
  logical                   :: write_out
  real(dp)                  :: time, dt, photoi_prev_time
  real(dp)                  :: memory_limit_GB = 16.0_dp
  type(CFG_t)               :: cfg            ! The configuration for the simulation
  type(af_t)                :: tree           ! This contains the full grid information
  type(af_t)                :: tree_copy      ! Used when reading a tree from a file
  type(mg_t)                :: mg             ! Multigrid option structure
  type(ref_info_t)          :: ref_info       ! Contains info about refinement changes
  integer                   :: output_cnt = 0 ! Number of output files written
  character(len=string_len) :: restart_from_file = undefined_str
  real(dp)                  :: max_field, initial_streamer_pos
  type(af_loc_t)            :: loc_field, loc_field_initial
  real(dp), dimension(NDIM) :: loc_field_coord, loc_field_initial_coord

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
     do s = 0, advance_num_states-1
        call af_set_cc_methods(tree, species_itree(i)+s, &
             af_bc_neumann_zero, af_gc_interp_lim, ST_prolongation_method)
     end do
  end do

  if (.not. gas_constant_density) then
     call af_set_cc_methods(tree, i_gas_dens, &
          af_bc_neumann_zero, af_gc_interp, ST_prolongation_method)
  end if

  if (photoi_enabled) then
     call photoi_set_methods(tree)
  end if

  do i = 1, tree%n_var_cell
     if (tree%cc_write_output(i) .and. .not. tree%has_cc_method(i)) then
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
     dt               = advance_max_dt
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

     call advance(tree, mg, dt, time)
     dt = advance_max_dt
     global_time = time

     if (advance_max_dt < dt_min) then
        print *, "ST_dt getting too small, instability?", advance_max_dt
        call output_status(tree, time, wc_time, it, dt)
        if (.not. write_out) then
           write_out = .true.
           output_cnt = output_cnt + 1
        end if
     end if

     if (write_out) then
        call output_write(tree, output_cnt, wc_time, write_sim_data)
     end if

     if (advance_max_dt < dt_min) error stop "dt too small"

     if (mod(it, refine_per_steps) == 0) then
        ! Restrict species, for the ghost cells near refinement boundaries
        do i = n_gas_species+1, n_species
           call af_restrict_tree(tree, species_itree(i))
        end do

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
    use m_advance_base

    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout)  :: tree
    type(mg_t), intent(inout)  :: mg

    call user_initialize(cfg, tree)
    call advance_base_initialize(cfg)
    call table_data_initialize(cfg)
    call gas_initialize(tree, cfg)
    call transport_data_initialize(cfg)
    call chemistry_initialize(tree, cfg)
    call ST_initialize(tree, cfg, NDIM)
    call photoi_initialize(tree, cfg)

    call dt_initialize(cfg)
    call advance_set_max_dt(dt_min)
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
  end subroutine set_initial_conditions

  subroutine write_sim_data(my_unit)
    integer, intent(in) :: my_unit

    write(my_unit) output_cnt
    write(my_unit) time
    write(my_unit) global_time
    write(my_unit) photoi_prev_time
    write(my_unit) advance_max_dt
  end subroutine write_sim_data

  subroutine read_sim_data(my_unit)
    integer, intent(in) :: my_unit
    real(dp)            :: tmp

    read(my_unit) output_cnt
    read(my_unit) time
    read(my_unit) global_time
    read(my_unit) photoi_prev_time
    read(my_unit) tmp

    call advance_set_max_dt(tmp)
    dt = tmp
  end subroutine read_sim_data

end program streamer
