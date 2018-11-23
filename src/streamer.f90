#include "../afivo/src/cpp_macros.h"
!> Program to perform NDIMd streamer simulations in Cartesian and cylindrical coordinates
program streamer

  use m_af_all
  use m_streamer
  use m_field
  use m_init_cond
  use m_refine
  use m_photoi
  use m_chemistry
  use m_gas
  use m_fluid_lfa

  implicit none

  integer, parameter        :: int8 = selected_int_kind(18)
  integer(int8)             :: t_start, t_current, count_rate
  real(dp)                  :: wc_time, inv_count_rate, time_last_print
  integer                   :: i, it
  character(len=string_len) :: fname
  character(len=name_len)   :: prolong_method
  character(len=name_len)   :: time_integrator
  logical                   :: write_out
  real(dp)                  :: dt_prev, photoi_prev_time
  type(CFG_t)               :: cfg  ! The configuration for the simulation
  type(af_t)                :: tree ! This contains the full grid information
  type(mg_t)                :: mg   ! Multigrid option struct
  type(ref_info_t)          :: ref_info

  ! Method used to prolong (interpolate) densities
  procedure(af_subr_prolong), pointer :: prolong_density => null()

  integer :: output_cnt = 0 ! Number of output files written

  call CFG_update_from_arguments(cfg)

  call gas_init(cfg)
  call chemistry_init(tree, cfg)
  call ST_initialize(tree, cfg, NDIM)
  call photoi_initialize(tree, cfg)

  call ST_load_transport_data(cfg)
  call refine_initialize(cfg)
  call field_initialize(cfg, mg)
  call init_cond_initialize(cfg, NDIM)

  time_integrator = "heuns_method"
  call CFG_add_get(cfg, "time_integrator", time_integrator, &
       "Time integrator (forward_euler, heuns_method)")

  prolong_method = "limit"
  call CFG_add_get(cfg, "prolong_density", prolong_method, &
       "Density prolongation method (limit, linear, linear_cons, sparse)")
  select case (prolong_method)
  case ("limit")
     prolong_density => af_prolong_limit
  case ("linear")
     prolong_density => af_prolong_linear
  case ("linear_cons")
     prolong_density => af_prolong_linear_cons
  case ("sparse")
     prolong_density => af_prolong_sparse
  case default
     error stop "Unknown prolong_density method"
  end select

  call check_path_writable(ST_output_dir)

  fname = trim(ST_output_dir) // "/" // trim(ST_simulation_name) // "_out.cfg"
  call CFG_write(cfg, trim(fname))

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the multigrid options. First define the variables to use
  mg%i_phi = i_phi
  mg%i_tmp = i_electric_fld
  mg%i_rhs = i_rhs

  ! This automatically handles cylindrical symmetry
  mg%box_op => mg_auto_op
  mg%box_gsrb => mg_auto_gsrb
  mg%box_corr => mg_auto_corr

  ! This routine always needs to be called when using multigrid
  call mg_init_mg(mg)

  do i = 1, n_species
     call af_set_cc_methods(tree, species_ix(i), &
          af_bc_neumann_zero, af_gc_interp_lim, prolong_density)
  end do

  call af_set_cc_methods(tree, i_phi, mg%sides_bc, mg%sides_rb)
  call af_set_cc_methods(tree, i_electric_fld, &
       af_bc_neumann_zero, af_gc_interp)

  if (ST_output_src_term) then
     call af_set_cc_methods(tree, i_src, &
          af_bc_neumann_zero, af_gc_interp)
  end if

  if (photoi_enabled) then
     call photoi_set_methods(tree)
  end if

  do i = 1, tree%n_var_cell
     if (tree%cc_write_output(i) .and. .not. tree%has_cc_method(i)) then
        call af_set_cc_methods(tree, i, af_bc_neumann_zero, af_gc_interp)
     end if
  end do

  output_cnt       = 0 ! Number of output files written
  ST_time          = 0 ! Simulation time (all times are in s)
  photoi_prev_time = 0 ! Time of last photoionization computation

  ! Set up the initial conditions
  do
     call af_loop_box(tree, init_cond_set_box)
     call field_compute(tree, mg, .false.)
     call af_adjust_refinement(tree, refine_routine, ref_info, &
          refine_buffer_width, .true.)
     if (ref_info%n_add == 0) exit
  end do

  call init_cond_stochastic_density(tree)

  print *, "Number of threads", af_get_max_threads()
  call af_print_info(tree)

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
        call photoi_set_src(tree, ST_time - photoi_prev_time)
        photoi_prev_time = ST_time
     end if

     dt_matrix(1:ST_dt_num_cond, :) = ST_dt_max      ! Maximum time step

     select case (time_integrator)
     case ("forward_euler")
        call forward_euler(tree, ST_dt, 0, 0, .true.)
        call restrict_flux_species(tree, 0)
        call field_compute(tree, mg, .true.)
        ST_time = ST_time + ST_dt
     case ("rk2")
        call forward_euler(tree, 0.5_dp * ST_dt, 0, 1, .false.)
        call restrict_flux_species(tree, 1)
        ST_time = ST_time + 0.5_dp * ST_dt
        call field_compute(tree, mg, .true.)
        call forward_euler(tree, ST_dt, 1, 0, .true.)
        call restrict_flux_species(tree, 0)
        call field_compute(tree, mg, .true.)
     case ("heuns_method")
        call forward_euler(tree, ST_dt, 0, 1, .false.)
        call restrict_flux_species(tree, 1)
        ST_time = ST_time + ST_dt
        call field_compute(tree, mg, .true.)
        call forward_euler(tree, ST_dt, 1, 1, .true.)
        call combine_substeps(tree, species_ix(1:n_species), &
             [0, 1], [0.5_dp, 0.5_dp], 0)
        call restrict_flux_species(tree, 0)
        call field_compute(tree, mg, .true.)
     end select

     ! Determine next time step
     ST_dt   = min(2 * dt_prev, ST_dt_safety_factor * &
          minval(dt_matrix(1:ST_dt_num_cond, :)))
     dt_prev = ST_dt

     if (ST_dt < ST_dt_min .and. .not. write_out) then
        print *, "ST_dt getting too small, instability?", ST_dt
        write_out = .true.
        output_cnt = output_cnt + 1
     end if

     if (write_out) then
        if (ST_silo_write) then
           ! Because the mesh could have changed
           if (photoi_enabled) call photoi_set_src(tree, ST_dt)

           do i = 1, tree%n_var_cell
              if (tree%cc_write_output(i)) then
                 call af_restrict_tree(tree, i)
                 call af_gc_tree(tree, i)
              end if
           end do

           write(fname, "(A,I6.6)") trim(ST_simulation_name) // "_", output_cnt
           call af_write_silo(tree, fname, output_cnt, ST_time, dir=ST_output_dir)
        end if

        if (ST_datfile_write) then
           call af_write_tree(tree, fname, ST_output_dir)
        end if

        write(fname, "(A,I6.6)") trim(ST_simulation_name) // "_log.txt"
        call write_log_file(tree, fname, output_cnt, ST_output_dir)

        if (ST_plane_write) then
           write(fname, "(A,I6.6)") trim(ST_simulation_name) // &
                "_plane_", output_cnt
           call af_write_plane(tree, fname, [ST_plane_ivar], &
                ST_plane_rmin * ST_domain_len, &
                ST_plane_rmax * ST_domain_len, &
                ST_plane_npixels, ST_output_dir)
        end if

        if (ST_lineout_write) then
           write(fname, "(A,I6.6)") trim(ST_simulation_name) // &
                "_line_", output_cnt
           call af_write_line(tree, trim(fname), &
                [i_electron, i_1pos_ion, i_phi, i_electric_fld], &
                r_min = ST_lineout_rmin(1:NDIM) * ST_domain_len, &
                r_max = ST_lineout_rmax(1:NDIM) * ST_domain_len, &
                n_points=ST_lineout_npoints, dir=ST_output_dir)
        end if
     end if

     if (ST_dt < ST_dt_min) error stop "dt too small"

     if (mod(it, refine_per_steps) == 0) then
        call af_gc_tree(tree, species_ix(1:n_species))
        call af_adjust_refinement(tree, refine_routine, ref_info, &
             refine_buffer_width, .true.)

        if (ref_info%n_add > 0 .or. ref_info%n_rm > 0) then
           ! Compute the field on the new mesh
           call field_compute(tree, mg, .true.)
        end if
     end if
  end do

  call print_status()
  call af_destroy(tree)

contains

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

  !> Restrict species for which we compute fluxes near refinement boundaries,
  !> for the coarse grid ghost cells
  subroutine restrict_flux_species(tree, s_out)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: s_out
    integer                   :: lvl, i, id, p_id

    !$omp parallel private(lvl, i, id, p_id)
    do lvl = tree%lowest_lvl, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          p_id = tree%boxes(id)%parent
          if (p_id > af_no_box .and. &
               any(tree%boxes(id)%neighbors == af_no_box)) then
             call af_restrict_box_vars(tree%boxes(id), tree%boxes(p_id), &
                  flux_species + s_out)
          end if
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine restrict_flux_species

  subroutine write_log_file(tree, filename, out_cnt, dir)
    type(af_t), intent(in)      :: tree
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: out_cnt
    character(len=*), intent(in) :: dir
    character(len=string_len)       :: fname
    character(len=30), save      :: fmt
    integer, parameter           :: my_unit = 123
    real(dp)                     :: velocity, dt
    real(dp), save               :: prev_pos(NDIM) = 0
    real(dp)                     :: sum_elec, sum_pos_ion
    real(dp)                     :: max_elec, max_field, max_Er
    type(af_loc_t)              :: loc_elec, loc_field, loc_Er

    call af_prepend_directory(filename, dir, fname)

    call af_tree_sum_cc(tree, i_electron, sum_elec)
    call af_tree_sum_cc(tree, i_1pos_ion, sum_pos_ion)
    call af_tree_max_cc(tree, i_electron, max_elec, loc_elec)
    call af_tree_max_cc(tree, i_electric_fld, max_field, loc_field)
    call af_tree_max_fc(tree, 1, electric_fld, max_Er, loc_Er)

    dt = ST_dt_safety_factor * minval(dt_matrix(1:ST_dt_num_cond, :))

    if (out_cnt == 1) then
       open(my_unit, file=trim(fname), action="write")
#if NDIM == 2
       write(my_unit, *) "# it time dt v sum(n_e) sum(n_i) ", &
            "max(E) x y max(n_e) x y max(E_r) x y wc_time n_cells min(dx) highest(lvl)"
       fmt = "(I6,15E16.8,I12,1E16.8,I3)"
#elif NDIM == 3
       write(my_unit, *) "# it time dt v sum(n_e) sum(n_i) ", &
            "max(E) x y z max(n_e) x y z wc_time n_cells min(dx) highest(lvl)"
       fmt = "(I6,14E16.8,I12,1E16.8,I3)"
#endif
       close(my_unit)

       ! Start with velocity zero
       prev_pos = af_r_loc(tree, loc_field)
    end if

    velocity = norm2(af_r_loc(tree, loc_field) - prev_pos) / ST_dt_output
    prev_pos = af_r_loc(tree, loc_field)

    open(my_unit, file=trim(fname), action="write", &
         position="append")
#if NDIM == 2
    write(my_unit, fmt) out_cnt, ST_time, dt, velocity, sum_elec, &
         sum_pos_ion, max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec), max_Er, af_r_loc(tree, loc_Er), &
         wc_time, af_num_cells_used(tree), af_min_dr(tree),tree%highest_lvl
#elif NDIM == 3
    write(my_unit, fmt) out_cnt, ST_time, dt, velocity, sum_elec, &
         sum_pos_ion, max_field, af_r_loc(tree, loc_field), max_elec, &
         af_r_loc(tree, loc_elec), wc_time, af_num_cells_used(tree), af_min_dr(tree),tree%highest_lvl
#endif
    close(my_unit)

  end subroutine write_log_file

  subroutine print_status()
    write(*, "(F7.2,A,I0,A,E10.3,A,E10.3,A,E10.3,A,E10.3)") &
         100 * ST_time / ST_end_time, "% it: ", it, &
         " t:", ST_time, " dt:", ST_dt, " wc:", wc_time, &
         " ncell:", real(af_num_cells_used(tree), dp)
  end subroutine print_status

  subroutine check_path_writable(pathname)
    character(len=*), intent(in) :: pathname
    integer                      :: my_unit, iostate
    open(newunit=my_unit, file=trim(pathname)//"/DUMMY", iostat=iostate)
    if (iostate /= 0) then
       print *, "Output directory: " // trim(pathname)
       error stop "Directory not writable (does it exist?)"
    else
       close(my_unit, status='delete')
    end if
  end subroutine check_path_writable

  subroutine combine_substeps(tree, ivs, in_steps, coeffs, out_step)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: ivs(:)
    integer, intent(in)       :: in_steps(:)
    real(dp), intent(in)      :: coeffs(:)
    integer, intent(in)       :: out_step
    integer                   :: lvl, i, id, n, nc
    real(dp), allocatable     :: tmp(DTIMES(:), :)

    nc = tree%n_cell

    !$omp parallel private(lvl, i, id, tmp, n)
    allocate(tmp(DTIMES(0:nc+1), size(ivs)))

    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)

          tmp = 0.0_dp
          do n = 1, size(in_steps)
             tmp = tmp + coeffs(n) * &
                  tree%boxes(id)%cc(DTIMES(:), ivs+in_steps(n))
          end do
          tree%boxes(id)%cc(DTIMES(:), ivs+out_step) = tmp
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine combine_substeps

end program streamer
