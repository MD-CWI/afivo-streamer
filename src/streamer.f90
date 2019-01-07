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
  real(dp)                  :: wc_time, inv_count_rate, time_last_print
  integer                   :: i, s, it
  logical                   :: write_out
  real(dp)                  :: time, dt, photoi_prev_time
  type(CFG_t)               :: cfg            ! The configuration for the simulation
  type(af_t)                :: tree           ! This contains the full grid information
  type(mg_t)                :: mg             ! Multigrid option struct
  type(ref_info_t)          :: ref_info
  integer                   :: output_cnt = 0 ! Number of output files written

  call CFG_update_from_arguments(cfg)

  call initialize_modules(cfg, tree, mg)

  call CFG_write(cfg, trim(output_name) // "_out.cfg")

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Specify default methods for all the variables
  do i = 1, n_species
     do s = 0, advance_num_states-1
        call af_set_cc_methods(tree, species_ix(i)+s, &
             af_bc_neumann_zero, af_gc_interp_lim)
     end do
  end do

  if (photoi_enabled) then
     call photoi_set_methods(tree)
  end if

  do i = 1, tree%n_var_cell
     if (tree%cc_write_output(i) .and. .not. tree%has_cc_method(i)) then
        call af_set_cc_methods(tree, i, af_bc_neumann_zero, af_gc_interp)
     end if
  end do

  output_cnt       = 0      ! Number of output files written
  time             = 0.0_dp ! Simulation time (all times are in s)
  global_time      = time
  photoi_prev_time = 0.0_dp ! Time of last photoionization computation

  ! Set up the initial conditions
  call set_initial_conditions(tree, mg)

  print *, "Number of threads", af_get_max_threads()
  call af_print_info(tree)

  ! Initial wall clock time
  call system_clock(t_start, count_rate)
  inv_count_rate = 1.0_dp / count_rate
  time_last_print = -1e10_dp

  do it = 1, huge(1)-1
     if (time >= ST_end_time) exit

     dt = advance_max_dt

     ! Update wall clock time
     call system_clock(t_current)
     wc_time = (t_current - t_start) * inv_count_rate

     ! Every ST_print_status_interval, print some info about progress
     if (wc_time - time_last_print > output_status_delay) then
        call output_status(tree, time, wc_time, it, dt)
        time_last_print = wc_time
     end if

     ! Every output_dt, write output
     if (output_cnt * output_dt <= time + dt) then
        write_out  = .true.
        dt         = output_cnt * output_dt - time
        output_cnt = output_cnt + 1
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
        if (.not. write_out) then
           write_out = .true.
           output_cnt = output_cnt + 1
        end if
     end if

     if (write_out) then
        call output_write(tree, output_cnt, wc_time)
     end if

     if (advance_max_dt < dt_min) error stop "dt too small"

     if (mod(it, refine_per_steps) == 0) then
        call af_gc_tree(tree, species_ix(1:n_species))

        if (associated(user_refine)) then
           call af_adjust_refinement(tree, user_refine, ref_info, &
                refine_buffer_width, .true.)
        else
           call af_adjust_refinement(tree, refine_routine, ref_info, &
                refine_buffer_width, .true.)
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

    call advance_base_initialize(cfg)
    call table_data_initialize(cfg)
    call gas_initialize(cfg)
    call transport_data_initialize(cfg)
    call chemistry_initialize(tree, cfg)
    call ST_initialize(tree, cfg, NDIM)
    call photoi_initialize(tree, cfg)

    call dt_initialize(cfg)
    call advance_initialize()
    call refine_initialize(cfg)
    call field_initialize(tree, cfg, mg)
    call init_cond_initialize(cfg)
    call output_initialize(tree, cfg)
    call user_initialize(cfg, tree)

  end subroutine initialize_modules

  !> Initialize the AMR tree
  subroutine init_tree(tree)
    type(af_t), intent(inout) :: tree

    ! Variables used below to initialize tree
    real(dp)                  :: dr
    integer                   :: id
    integer                   :: ix_list(NDIM, 1) ! Spatial indices of initial boxes
    integer                   :: nb_list(2*NDIM, 1) ! Index of neighbors

    dr = ST_domain_len / ST_box_size

    ! Initialize tree
    if (ST_cylindrical) then
       call af_init(tree, ST_box_size, dr, coarsen_to=2, coord=af_cyl)
    else
       call af_init(tree, ST_box_size, dr, coarsen_to=2)
    end if

    ! Set up geometry
    id             = 1          ! One box ...
    ix_list(:, id) = 1          ! With index 1,1 ...

    if (ST_periodic) then
       if (ST_cylindrical) error stop "Cannot have periodic cylindrical domain"
       nb_list(:, 1)               = -1
       nb_list(af_neighb_lowx, 1)  = 1
       nb_list(af_neighb_highx, 1) = 1
#if NDIM == 3
       nb_list(af_neighb_lowy, 1)  = 1
       nb_list(af_neighb_highy, 1) = 1
#endif
       call af_set_base(tree, 1, ix_list, nb_list)
    else
       ! Create the base mesh
       call af_set_base(tree, 1, ix_list)
    end if

  end subroutine init_tree

  subroutine set_initial_conditions(tree, mg)
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(in)    :: mg

    do
       call af_loop_box(tree, init_cond_set_box)

       if (associated(user_initial_conditions)) then
          call af_loop_box(tree, user_initial_conditions)
       else if (ST_use_dielectric) then
          error stop "use_dielectric requires user_initial_conditions to be set"
       end if

       call field_compute(tree, mg, 0, time, .false.)

       if (associated(user_refine)) then
          call af_adjust_refinement(tree, user_refine, ref_info, &
               refine_buffer_width, .true.)
       else
          call af_adjust_refinement(tree, refine_routine, ref_info, &
               refine_buffer_width, .true.)
       end if

       if (ref_info%n_add == 0) exit
    end do
  end subroutine set_initial_conditions

end program streamer
