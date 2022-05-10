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
  use m_dielectric

  implicit none

  integer, parameter        :: int8       = selected_int_kind(18)
  integer(int8)             :: t_start, t_current, count_rate
  real(dp)                  :: wc_time, inv_count_rate
  real(dp)                  :: time_last_print, time_last_output
  integer                   :: i, it, coord_type, box_bytes
  integer, allocatable      :: ref_links(:, :)
  logical                   :: write_out, evolve_electrons
  real(dp)                  :: time, dt, dt_lim, photoi_prev_time
  real(dp)                  :: dt_gas_lim
  real(dp)                  :: memory_limit_GB = 16.0_dp
  type(af_t)                :: tree_copy      ! Used when reading a tree from a file
  type(ref_info_t)          :: ref_info       ! Contains info about refinement changes
  integer                   :: output_cnt = 0 ! Number of output files written
  character(len=string_len) :: restart_from_file = undefined_str
  real(dp)                  :: max_field, initial_streamer_pos
  type(af_loc_t)            :: loc_field, loc_field_initial
  real(dp), dimension(NDIM) :: loc_field_coord, loc_field_initial_coord
  real(dp)                  :: t_field_off = -1
  character(len=string_len) :: electron_bc = "standard"
  real(dp)                  :: output_dt_interpulse_fac = 1.0
  real(dp)                  :: breakdown_field_Td

  !> The configuration for the simulation
  type(CFG_t) :: cfg
  !> This contains the full grid information
  type(af_t)  :: tree

  call print_program_name()

  ! Parse command line configuration files and options
  call CFG_update_from_arguments(cfg)

  call CFG_add_get(cfg, "restart_from_file", restart_from_file, &
       "If set, restart simulation from a previous .dat file")
  call CFG_add_get(cfg, "memory_limit_GB", memory_limit_GB, &
       "Memory limit (GB)")

  call CFG_add_get(cfg, "t_field_off", t_field_off, &
       "Time after which the applied field is turned off.")
  call CFG_add_get(cfg, "electron_bc", electron_bc, &
       "BC used for electron: standard (means using bc_species variable), dirichlet_custom, &
        neuman_custom.")

  call CFG_add_get(cfg, "output_dt_interpulse_fac", output_dt_interpulse_fac, &
   "Factor which is multiplied with output_dt resulting in an output_dt when applied voltage is 0.")

  call initialize_modules(cfg, tree, mg, restart_from_file /= undefined_str)

  call CFG_write(cfg, trim(output_name) // "_out.cfg", custom_first=.true.)

  call chemistry_get_breakdown_field(breakdown_field_Td, 1.0e3_dp)
  write(*, '(A,E12.4)') " Estimated breakdown field (Td): ", breakdown_field_Td

  ! Specify default methods for all the variables
  do i = n_gas_species+1, n_species
      if (i /= ix_electron) then
         call af_set_cc_methods(tree, species_itree(i), &
               bc_species, af_gc_interp_lim, ST_prolongation_method)
      end if
  end do

  if (electron_bc == "standard") then
      print *, "Using standard electron_BC"
      call af_set_cc_methods(tree, species_itree(ix_electron), &
               bc_species, af_gc_interp_lim, ST_prolongation_method)
  else if (electron_bc == "dirichlet_custom") then
      print *, "Using custom dirichlet 0 electron_BC"
      call af_set_cc_methods(tree, species_itree(ix_electron), &
               rb=af_gc_interp_lim, prolong=ST_prolongation_method, bc_custom=dirichlet_custom)
  else if (electron_bc == "neuman_custom") then
      print *, "Using custom neuman 0 electron_BC"
      call af_set_cc_methods(tree, species_itree(ix_electron), &
               rb=af_gc_interp_lim, prolong=ST_prolongation_method, bc_custom=neuman_custom)
  else if (electron_bc == "robin_custom") then
      print *, "Using custom robin electron_BC"
      call af_set_cc_methods(tree, species_itree(ix_electron), &
               rb=af_gc_interp_lim, prolong=ST_prolongation_method, bc_custom=robin_custom)
  else if (electron_bc == "outflow") then
      print *, "Using custom outflow electron_BC"
      call af_set_cc_methods(tree, species_itree(ix_electron), &
               rb=af_gc_interp_lim, prolong=ST_prolongation_method, bc_custom=outflow_custom)
  end if

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

     if (ST_use_dielectric) error stop "Restarting not support with dielectric"

     ! @todo more consistency checks

     ! This routine always needs to be called when using multigrid
     call mg_init(tree, mg)

     ! We want to set the initial seed again
     call set_initial_conditions_restart(tree, mg)
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
  call af_stencil_print_info(tree)

  ! Initial wall clock time
  call system_clock(t_start, count_rate)
  inv_count_rate   = 1.0_dp / count_rate
  time_last_print  = -1e10_dp
  time_last_output = time

  do it = 1, huge(1)-1
     if (ST_use_end_time .and. time >= ST_end_time) exit

     ! If the applied voltage is 0. turn on the drt limit flux fix
     if (ST_use_drt_fix_voltage_off) then
      if (field_voltage == 0) then
         ST_drt_limit_flux = .true.
      else
         ST_drt_limit_flux = .false.
      end if
     end if

     ! Turn off the applied field
     if (t_field_off >= 0 .and. time >= t_field_off) then
        call field_set_voltage_externally(0.0_dp)
     end if
     
     if (associated(user_generic_method)) then
        call user_generic_method(tree, time)
     end if

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
     if (field_voltage == 0) then
         if (time + dt >= time_last_output + (output_dt_interpulse_fac * output_dt)) then
            write_out        = .true.
            dt               = time_last_output + (output_dt_interpulse_fac * output_dt) - time
            time_last_output = time_last_output + (output_dt_interpulse_fac * output_dt)
            output_cnt       = output_cnt + 1
         else
            write_out = .false.
         end if
      else if (time + dt >= time_last_output + output_dt) then
         write_out        = .true.
         dt               = time_last_output + output_dt - time
         time_last_output = time_last_output + output_dt
         output_cnt       = output_cnt + 1
      else
         write_out = .false.
      end if

     !if (time + dt >= time_last_output + output_dt) then
     !   write_out        = .true.
     !   dt               = time_last_output + output_dt - time
     !   time_last_output = time_last_output + output_dt
     !   output_cnt       = output_cnt + 1
     !else
     !   write_out = .false.
     !end if

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

        if (ST_use_dielectric) then
           call af_advance(tree, dt, dt_lim, time, &
                species_itree(n_gas_species+1:n_species), &
                time_integrator, forward_euler, &
                dielectric_combine_substeps)
        else
           call af_advance(tree, dt, dt_lim, time, &
                species_itree(n_gas_species+1:n_species), &
             time_integrator, forward_euler)
        end if

        ! Make sure field is available for latest time state
        call field_compute(tree, mg, 0, time, .true.)
        if (compute_power_density) call compute_total_energy_density(tree, dt)
        if (gas_dynamics) then
            !print *, "Dt vector: ", dt
            call coupling_add_fluid_source(tree, dt)
        end if
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
       ! To store surface charge and the photon flux
       n_surface_variables = af_advance_num_steps(time_integrator) + 1
       call surface_initialize(tree, i_eps, diel, n_surface_variables)
    end if

    do n = 1, 100

       ! This placement is so that users can set epsilon before the coarse grid
       ! solver is constructed
       if ((n == 1) .and. (mg%initialized .eqv. .false.)) call mg_init(tree, mg)

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

  subroutine set_initial_conditions_restart(tree, mg)
   type(af_t), intent(inout) :: tree
    type(mg_t), intent(inout) :: mg
    integer                   :: n

    do n = 1, 100
       call af_loop_box(tree, init_cond_restart_set_box)

       if (associated(user_initial_conditions)) then
          call af_loop_box(tree, user_initial_conditions)
       else if (ST_use_dielectric) then
          error stop "use_dielectric requires user_initial_conditions to be set"
       end if

       ! This placement is so that users can set epsilon before the coarse grid
       ! solver is constructed
       if ((n == 1) .and. (mg%initialized .eqv. .false.)) call mg_init(tree, mg)

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
  end subroutine set_initial_conditions_restart


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


   !> With this method we can set ghost cells manually
  subroutine custom_boundary_method(box, nb, iv, n_gc, cc)
   type(box_t), intent(inout) :: box     !< Box that needs b.c.
   integer, intent(in)     :: nb      !< Direction
   integer, intent(in)     :: iv      !< Index of variable
   integer, intent(in)     :: n_gc !< Number of ghost cells to fill
   !> If n_gc > 1, fill ghost cell values in this array instead of box%cc
   real(dp), intent(inout), optional :: &
        cc(DTIMES(1-n_gc:box%n_cell+n_gc))

   integer :: lo(NDIM), hi(NDIM)
   integer :: lo_inside(NDIM), hi_inside(NDIM)
   integer :: nc

   nc = box%n_cell

#if NDIM == 2


   if (n_gc >= 3) error stop "not implemented"


   !if (iv == ix_electron) then
      if (n_gc == 2) then
         ! #define DTIMES(TXT) TXT, TXT
         !call af_get_index_bc_outside(nb, box%n_cell, 2, lo, hi)
         !print *, "2 ghost cells: ", iv, nb, iv, n_gc, lo, hi
         !cc(DSLICE(lo, hi)) = 200e3
         select case (nb)
            case (af_neighb_lowx)
               cc(0, 1:nc) = box%cc(1, 1:nc, iv)
               cc(-1, 1:nc) = 0
            case (af_neighb_highx)
               cc(nc + 1, 1:nc) = box%cc(nc, 1:nc, iv)
               cc(nc + 2, 1:nc) = 0
            case (af_neighb_lowy)
               cc(1:nc, 0) = box%cc(1:nc, 1, iv)
               cc(1:nc, -1) = 0
            case (af_neighb_highy)
               cc(1:nc, nc + 1) = box%cc(1:nc, nc, iv)
               cc(1:nc, nc + 2) = 0
         end select
      else if (n_gc == 1) then

         ! Get ghost cell index range
         !call af_get_index_bc_outside(nb, box%n_cell, 1, lo, hi)
         !call af_get_index_bc_inside(nb, box%n_cell, 1, lo_inside, hi_inside)
         select case (nb)
            case (af_neighb_lowx)
               !print *, "af_neigb_lowx: ", nb
               box%cc(0, 1:nc, iv)    = box%cc(1, 1:nc, iv)
               !box%cc(DSLICE(lo,hi), iv) = box%cc(DSLICE(lo_inside,hi_inside), iv)
               !box%cc(DSLICE(lo_inside,hi_inside), iv) = -4
            case (af_neighb_highx)
               !print *, "af_neigb_highx: ", nb
               box%cc(nc+1, 1:nc, iv) = box%cc(nc, 1:nc, iv)
               !box%cc(DSLICE(lo,hi), iv) = box%cc(DSLICE(lo_inside,hi_inside), iv)
               !box%cc(DSLICE(lo_inside,hi_inside), iv) = -5
            case (af_neighb_lowy)
               !print *, "af_neigb_lowy: ", nb
               box%cc(1:nc, 0, iv)    = box%cc(1:nc, 1, iv)
               !box%cc(DSLICE(lo,hi), iv) = -1 * box%cc(DSLICE(lo_inside,hi_inside), iv)
               !box%cc(DSLICE(lo_inside,hi_inside), iv) = -6
            case (af_neighb_highy)
               !print *, "af_neigb_highy: ", nb
               !print *, "1 ghost cell: ", lo, hi
               !print *, it, box%cc(DSLICE(lo,hi), iv)
               !print *, it, box%cc(DSLICE(lo_inside,hi_inside), iv)
               box%cc(1:nc, nc+1, iv) = box%cc(1:nc, nc, iv)

               !box%cc(DSLICE(lo,hi), iv) = -1 * box%cc(DSLICE(lo_inside,hi_inside), iv)
               !box%cc(DSLICE(lo_inside,hi_inside), iv) = -7
         end select

         ! Set all ghost cells to a scalar value
         ! #define DSLICE(lo,hi) lo(1):hi(1), lo(2):hi(2)
      end if
   !end if
#endif
  end subroutine custom_boundary_method



!> With this method we can set ghost cells manually
  subroutine dirichlet_custom(box, nb, iv, n_gc, cc)
   type(box_t), intent(inout) :: box     !< Box that needs b.c.
   integer, intent(in)     :: nb      !< Direction
   integer, intent(in)     :: iv      !< Index of variable
   integer, intent(in)     :: n_gc !< Number of ghost cells to fill
   !> If n_gc > 1, fill ghost cell values in this array instead of box%cc
   real(dp), intent(inout), optional :: &
        cc(DTIMES(1-n_gc:box%n_cell+n_gc))
   integer :: nc

   nc = box%n_cell

#if NDIM == 2
   if (n_gc >= 3) error stop "not implemented"

   if (n_gc == 2) then
      select case (nb)
         case (af_neighb_lowx)
            !Hema thinks: cc(0, 1:nc) = -box%cc(1, 1:nc, iv)
            cc(0, 1:nc) = -box%cc(1, 1:nc, iv)
            !cc(0, 1:nc) = box%cc(1, 1:nc, iv)
            cc(-1, 1:nc) = 0
         case (af_neighb_highx)
            !Hema thinks: cc(nc+1, 1:nc) = -box%cc(nc, 1:nc, iv)
            cc(nc + 1, 1:nc) = -box%cc(nc, 1:nc, iv)
            !cc(nc + 1, 1:nc) = box%cc(nc, 1:nc, iv)
            cc(nc + 2, 1:nc) = 0
         case (af_neighb_lowy)
            cc(1:nc, 0) = -box%cc(1:nc, 1, iv)
            cc(1:nc, -1) = 0
         case (af_neighb_highy)
            cc(1:nc, nc + 1) = -box%cc(1:nc, nc, iv)
            cc(1:nc, nc + 2) = 0
      end select
   else if (n_gc == 1) then
      select case (nb)
         case (af_neighb_lowx)
            !box%cc(0, 1:nc, iv) = box%cc(1, 1:nc, iv)
            box%cc(0, 1:nc, iv) = - box%cc(1, 1:nc, iv)
         case (af_neighb_highx)
            !box%cc(nc + 1, 1:nc, iv) = box%cc(nc, 1:nc, iv)
            box%cc(nc + 1, 1:nc, iv) = -box%cc(nc, 1:nc, iv)
         case (af_neighb_lowy)
            box%cc(1:nc, 0, iv) = -box%cc(1:nc, 1, iv)
         case (af_neighb_highy)
            box%cc(1:nc, nc + 1, iv) = -box%cc(1:nc, nc, iv)
      end select
   end if
#endif
  end subroutine dirichlet_custom



!> With this method we can set ghost cells manually
  subroutine neuman_custom(box, nb, iv, n_gc, cc)
   type(box_t), intent(inout) :: box     !< Box that needs b.c.
   integer, intent(in)     :: nb      !< Direction
   integer, intent(in)     :: iv      !< Index of variable
   integer, intent(in)     :: n_gc !< Number of ghost cells to fill
   !> If n_gc > 1, fill ghost cell values in this array instead of box%cc
   real(dp), intent(inout), optional :: &
        cc(DTIMES(1-n_gc:box%n_cell+n_gc))
   integer :: nc

   nc = box%n_cell

#if NDIM == 2
   if (n_gc >= 3) error stop "not implemented"

   if (n_gc == 2) then
      select case (nb)
         case (af_neighb_lowx)
            cc(0, 1:nc) = box%cc(1, 1:nc, iv)
            cc(-1, 1:nc) = 0
         case (af_neighb_highx)
            cc(nc + 1, 1:nc) = box%cc(nc, 1:nc, iv)
            cc(nc + 2, 1:nc) = 0
         case (af_neighb_lowy)
            cc(1:nc, 0) = box%cc(1:nc, 1, iv)
            cc(1:nc, -1) = 0
         case (af_neighb_highy)
            cc(1:nc, nc + 1) = box%cc(1:nc, nc, iv)
            cc(1:nc, nc + 2) = 0
      end select
   else if (n_gc == 1) then
      select case (nb)
         case (af_neighb_lowx)
            box%cc(0, 1:nc, iv) = box%cc(1, 1:nc, iv)
         case (af_neighb_highx)
            box%cc(nc + 1, 1:nc, iv) = box%cc(nc, 1:nc, iv)
         case (af_neighb_lowy)
            box%cc(1:nc, 0, iv) = box%cc(1:nc, 1, iv)
         case (af_neighb_highy)
            box%cc(1:nc, nc + 1, iv) = box%cc(1:nc, nc, iv)
      end select
   end if
#endif
  end subroutine neuman_custom


!> With this method we can set ghost cells manually
subroutine robin_custom(box, nb, iv, n_gc, cc)
   use m_transport_data
   use m_lookup_table
   type(box_t), intent(inout) :: box     !< Box that needs b.c.
   integer, intent(in)     :: nb      !< Direction
   integer, intent(in)     :: iv      !< Index of variable
   integer, intent(in)     :: n_gc !< Number of ghost cells to fill
   !> If n_gc > 1, fill ghost cell values in this array instead of box%cc
   real(dp), intent(inout), optional :: &
         cc(DTIMES(1-n_gc:box%n_cell+n_gc))
   integer :: nc
   real(dp) :: mobility_face(box%n_cell)  ! Mobility at each boundary face
   real(dp) :: diffusion_face(box%n_cell)  ! Diffusion at each boundary face
   real(dp) :: electric_field_y_face(box%n_cell)  ! Electric field at each boundary face
   real(dp) :: Td_electric_field_y_face(box%n_cell)  ! Electric field in Td at each boundary face
   real(dp) :: n_face(box%n_cell) ! density at boundary faces (we use upwinding)
   integer :: ix_electric_field_fc
   real(dp) :: robin_alpha(box%n_cell)
   real(dp) :: robin_beta(box%n_cell)
   real(dp) :: robin_gamma(box%n_cell)
   integer :: idx_gamma

   ! This implements a robin BC: alpha * n + beta * grad(n) = gamma
   ! We choose alpha and beta so that our robin BC prescribes the total flux J at the boundary
   ! alpha = +- mu * E  ; beta = -D  => J = gamma
   ! Outflow BC for electrons will have: gamma = (a - 1) * mu * E * n
   ! with a = 1 if E . surface_normal < 0 and 0 otherwise.
   ! In simpler terms: If the flow direction is AWAY from the boundary, gamma = 0  (dirichlet 0)
   !                   If the flow direction is TOWARDS the boundary, gamma = +- mu * E * n (+- based on charge of species, neuman 0)
   ! Normally we would need 'n' on the face of the boundary, but this is not stored explicitly.
   ! We need 'n' in alpha and in gamma, but for gamma it is only in the case if the flow is towards
   ! the boundary. For the gamma we can argue that an upwind scheme is fine and we use the cc 'n' value.
   ! For the alpha, I am not too sure, but either way both alpha and gamma should contain the same values 
   ! to get the cancelling of the advection term to obtain the neuman 0 condition so we will also use the upwind value

   nc = box%n_cell

   ix_electric_field_fc = af_find_fc_variable(tree, "field")


#if NDIM == 2
   if (n_gc >= 3) error stop "not implemented"

   if (n_gc == 2) then
      select case (nb)
         case (af_neighb_lowx)
            cc(0, 1:nc) = box%cc(1, 1:nc, iv)
            cc(-1, 1:nc) = 0
         case (af_neighb_highx)
            cc(nc + 1, 1:nc) = box%cc(nc, 1:nc, iv)
            cc(nc + 2, 1:nc) = 0
         case (af_neighb_lowy)
            print *, "Low y 2 GC"
            electric_field_y_face = box%fc(1:nc, 1, 2, ix_electric_field_fc)
            Td_electric_field_y_face = abs(electric_field_y_face) * SI_to_Townsend / gas_number_density

            mobility_face = LT_get_col(td_tbl, td_mobility, Td_electric_field_y_face) / gas_number_density
            diffusion_face = LT_get_col(td_tbl, td_diffusion, Td_electric_field_y_face) / gas_number_density
            !print * , "Mobility: ", mobility_face
            robin_alpha = -mobility_face * electric_field_y_face
            robin_beta = -diffusion_face


            if (time == 0) then
               !print *, "Time 0: ", time
               n_face = box%cc(1:nc, 1, iv)
            else
               do idx_gamma = 1, nc
                  if (electric_field_y_face(idx_gamma) > 0) then
                     ! E field is pointing up from the bottom boundary
                     ! This means electrons are flowing towards boundary
                     n_face(idx_gamma) = box%cc(idx_gamma, 1, iv)
                  else
                     n_face(idx_gamma) = box%cc(idx_gamma, 0, iv)
                  end if
               end do
            end if

            
            do idx_gamma = 1, nc
               if (electric_field_y_face(idx_gamma) > 0) then
                  ! E field is pointing up from the bottom boundary
                  ! This means electrons are flowing towards boundary
                  robin_gamma(idx_gamma) = -mobility_face(idx_gamma) * electric_field_y_face(idx_gamma) &
                                          * n_face(idx_gamma)
               else
                  robin_gamma(idx_gamma) = 0
               end if
            end do

            !print *, "n_face: ", n_face
            !print *, "Gamma: ", robin_gamma
            !print *, "Alpha: ", robin_alpha
            !print *, "E field face: ", electric_field_y_face
            !print *, "dr: ", box%dr
            !print *, "Beta: ", robin_beta

            cc(1:nc, 0) = (((robin_alpha * box%dr(2)) / robin_beta) * n_face) &
                           + cc(1:nc, 1) &
                           - ((robin_gamma * box%dr(2)) / robin_beta)
            
            print *, "Ghost: ", cc(1:nc, 0)

            cc(1:nc, -1) = 0
         case (af_neighb_highy)
            cc(1:nc, nc + 1) = box%cc(1:nc, nc, iv)
            cc(1:nc, nc + 2) = 0
      end select
   else if (n_gc == 1) then
      select case (nb)
         case (af_neighb_lowx)
            box%cc(0, 1:nc, iv) = box%cc(1, 1:nc, iv)
         case (af_neighb_highx)
            box%cc(nc + 1, 1:nc, iv) = box%cc(nc, 1:nc, iv)
         case (af_neighb_lowy)
            print *, "Low y"
            electric_field_y_face = box%fc(1:nc, 1, 2, ix_electric_field_fc)
            Td_electric_field_y_face = abs(electric_field_y_face) * SI_to_Townsend / gas_number_density

            mobility_face = LT_get_col(td_tbl, td_mobility, Td_electric_field_y_face) / gas_number_density
            !print * , "Mobility: ", mobility_face
            diffusion_face = LT_get_col(td_tbl, td_diffusion, Td_electric_field_y_face) / gas_number_density
            
            robin_alpha = -mobility_face * electric_field_y_face
            robin_beta = -diffusion_face


            if (time == 0) then
               !print *, "Time 0: ", time
               n_face = box%cc(1:nc, 1, iv)
            else
               do idx_gamma = 1, nc
                  if (electric_field_y_face(idx_gamma) > 0) then
                     ! E field is pointing up from the bottom boundary
                     ! This means electrons are flowing towards boundary
                     n_face(idx_gamma) = box%cc(idx_gamma, 1, iv)
                  else
                     n_face(idx_gamma) = box%cc(idx_gamma, 0, iv)
                  end if
               end do
            end if


            do idx_gamma = 1, nc
               if (electric_field_y_face(idx_gamma) > 0) then
                  ! E field is pointing up from the bottom boundary
                  ! This means electrons are flowing towards boundary
                  robin_gamma(idx_gamma) = -mobility_face(idx_gamma) * electric_field_y_face(idx_gamma) &
                                          * n_face(idx_gamma)
               else 
                  robin_gamma(idx_gamma) = 0
               end if
            end do

            print *, box%lvl, box%in_use, box%r_min
            !print *, "n_face: ", n_face
            !print *, "Gamma: ", robin_gamma
            !print *, "Alpha: ", robin_alpha
            !print *, "E field face: ", electric_field_y_face
            print *, "dr: ", box%dr
            !print *, "Beta: ", robin_beta

            box%cc(1:nc, 0, iv) = (((robin_alpha * box%dr(2)) / robin_beta) * n_face) &
                                    + box%cc(1:nc, 1, iv) &
                                    - ((robin_gamma * box%dr(2)) / robin_beta)

            !print *, "Ghost: ", box%cc(1:nc, 0, iv)
            !print *, "First factor: ", (robin_alpha * box%dr(2)) / robin_beta

         case (af_neighb_highy)
            box%cc(1:nc, nc + 1, iv) = box%cc(1:nc, nc, iv)
      end select
   end if
#endif
end subroutine robin_custom


!> With this method we can set ghost cells manually
subroutine outflow_custom(box, nb, iv, n_gc, cc)
   type(box_t), intent(inout) :: box     !< Box that needs b.c.
   integer, intent(in)     :: nb      !< Direction
   integer, intent(in)     :: iv      !< Index of variable
   integer, intent(in)     :: n_gc !< Number of ghost cells to fill
   !> If n_gc > 1, fill ghost cell values in this array instead of box%cc
   real(dp), intent(inout), optional :: &
        cc(DTIMES(1-n_gc:box%n_cell+n_gc))
   integer :: nc
   real(dp) :: electric_field_y_face(box%n_cell)  ! Electric field at each boundary face
   real(dp) :: electric_field_x_face(box%n_cell)  ! Electric field at each boundary face
   integer :: ix_electric_field_fc
   integer :: idx_gamma



   nc = box%n_cell
   ix_electric_field_fc = af_find_fc_variable(tree, "field")


#if NDIM == 2
   if (n_gc >= 3) error stop "not implemented"

   if (n_gc == 2) then
      select case (nb)
         case (af_neighb_lowx)
            electric_field_x_face = box%fc(1,1:nc,1,ix_electric_field_fc)
            do idx_gamma = 1,nc
               if (electric_field_x_face(idx_gamma) > 0) then
                  cc(idx_gamma,0) = box%cc(1,idx_gamma,iv)
               else
                  cc(idx_gamma,0) = -box%cc(1,idx_gamma,iv)
               end if
            end do
            !cc(0, 1:nc) = box%cc(1, 1:nc, iv)
            cc(-1, 1:nc) = 0
         case (af_neighb_highx)
            electric_field_x_face = box%fc(nc+1,1:nc,1,ix_electric_field_fc)
            do idx_gamma = 1,nc
               if (electric_field_x_face(idx_gamma) > 0) then
                  cc(idx_gamma,nc+1) = box%cc(nc,idx_gamma,iv)
               else
                  cc(idx_gamma,nc+1) = -box%cc(nc,idx_gamma,iv)
               end if
            end do
            !cc(nc + 1, 1:nc) = box%cc(nc, 1:nc, iv)
            cc(nc + 2, 1:nc) = 0
         case (af_neighb_lowy)

            electric_field_y_face = box%fc(1:nc, 1, 2, ix_electric_field_fc)

            do idx_gamma = 1, nc
               if (electric_field_y_face(idx_gamma) > 0) then
                  ! E field is pointing up from the bottom boundary
                  ! This means electrons are flowing towards boundary -> Neuman 0
                  cc(idx_gamma, 0) = box%cc(idx_gamma, 1, iv)
               else
                  ! Electrons flowing away from bottom boundary -> Dirichlet 0
                  cc(idx_gamma, 0) = -box%cc(idx_gamma, 1, iv)
               end if
            end do

            cc(1:nc, -1) = 0
         case (af_neighb_highy)

            electric_field_y_face = box%fc(1:nc, nc + 1, 2, ix_electric_field_fc)

            do idx_gamma = 1, nc
               if (electric_field_y_face(idx_gamma) > 0) then
                  ! E field is pointing up towards the top boundary
                  ! This means electrons are flowing away from the top boundary -> Dirichlet 0
                  cc(idx_gamma, nc + 1) = -box%cc(idx_gamma, nc, iv)
               else
                  ! Electrons flowing towards the top boundary -> Neuman 0
                  cc(idx_gamma, nc + 1) = box%cc(idx_gamma, nc, iv)
               end if
            end do

            cc(1:nc, nc + 2) = 0
      end select
   else if (n_gc == 1) then
      select case (nb)
         case (af_neighb_lowx)
            box%cc(0, 1:nc, iv) = box%cc(1, 1:nc, iv)
         case (af_neighb_highx)
            box%cc(nc + 1, 1:nc, iv) = box%cc(nc, 1:nc, iv)
         case (af_neighb_lowy)

            electric_field_y_face = box%fc(1:nc, 1, 2, ix_electric_field_fc)

            do idx_gamma = 1, nc
               if (electric_field_y_face(idx_gamma) > 0) then
                  ! E field is pointing up from the bottom boundary
                  ! This means electrons are flowing towards boundary -> Neuman 0
                  box%cc(idx_gamma, 0, iv) = box%cc(idx_gamma, 1, iv)
               else
                  ! Electrons flowing away from bottom boundary -> Dirichlet 0
                  box%cc(idx_gamma, 0, iv) = -box%cc(idx_gamma, 1, iv)
               end if
            end do

         case (af_neighb_highy)

            electric_field_y_face = box%fc(1:nc, nc + 1, 2, ix_electric_field_fc)

            do idx_gamma = 1, nc
               if (electric_field_y_face(idx_gamma) > 0) then
                  ! E field is pointing up towards the top boundary
                  ! This means electrons are flowing away from the top boundary -> Dirichlet 0
                  box%cc(idx_gamma, nc + 1, iv) = -box%cc(idx_gamma, nc, iv)
               else
                  ! Electrons flowing towards the top boundary -> Neuman 0
                  box%cc(idx_gamma, nc + 1, iv) = box%cc(idx_gamma, nc, iv)
               end if
            end do

      end select
   end if
#endif
  end subroutine outflow_custom
  
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
          box%cc(IJK, species_itree(n_gas_species+1:n_species)) = 0.0_dp
          !box%cc(IJK, gas_vars(i_rho)) = 0

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

end program streamer
