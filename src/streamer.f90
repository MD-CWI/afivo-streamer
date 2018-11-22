#include "../afivo/src/cpp_macros.h"
!> Program to perform NDIMd streamer simulations in Cartesian and cylindrical coordinates
program streamer

  use m_af_all
  use m_streamer
  use m_field
  use m_init_cond
  use m_photoi
  use m_chemistry
  use m_gas

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
          ST_refine_buffer_width, .true.)
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

     ST_dt_matrix(1:ST_dt_num_cond, :) = ST_dt_max      ! Maximum time step

     select case (time_integrator)
     case ("forward_euler")
        call forward_euler(tree, ST_dt, 0, 0)
        call restrict_flux_species(tree, 0)
        call field_compute(tree, mg, .true.)
        ST_time = ST_time + ST_dt
     case ("rk2")
        call forward_euler(tree, 0.5_dp * ST_dt, 0, 1)
        call restrict_flux_species(tree, 1)
        ST_time = ST_time + 0.5_dp * ST_dt
        call field_compute(tree, mg, .true.)
        call forward_euler(tree, ST_dt, 1, 0)
        call restrict_flux_species(tree, 0)
        call field_compute(tree, mg, .true.)
     case ("heuns_method")
        call forward_euler(tree, ST_dt, 0, 1)
        call restrict_flux_species(tree, 1)
        ST_time = ST_time + ST_dt
        call field_compute(tree, mg, .true.)
        call forward_euler(tree, ST_dt, 1, 1)
        call combine_substeps(tree, species_ix(1:n_species), [0, 1], [0.5_dp, 0.5_dp], 0)
        call restrict_flux_species(tree, 0)
        call field_compute(tree, mg, .true.)
     end select

     ! Determine next time step
     ST_dt   = min(2 * dt_prev, ST_dt_safety_factor * &
          minval(ST_dt_matrix(1:ST_dt_num_cond, :)))
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

     if (mod(it, ST_refine_per_steps) == 0) then
        call af_gc_tree(tree, species_ix(1:n_species))
        call af_adjust_refinement(tree, refine_routine, ref_info, &
             ST_refine_buffer_width, .true.)

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
       nb_list(:, 1)                = -1
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

  ! This routine sets the cell refinement flags for box
  subroutine refine_routine(box, cell_flags)
    use m_geometry
    use m_init_cond
    type(box_t), intent(in) :: box
    ! Refinement flags for the cells of the box
    integer, intent(out)     :: &
         cell_flags(DTIMES(box%n_cell))
    integer                  :: IJK, n, nc
    real(dp)                 :: dx, dx2
    real(dp)                 :: alpha, adx, fld, elec_dens
    real(dp)                 :: dist, rmin(NDIM), rmax(NDIM)

    nc      = box%n_cell
    dx      = box%dr
    dx2     = box%dr**2

    do KJI_DO(1,nc)
       fld   = box%cc(IJK, i_electric_fld) * SI_to_Townsend
       alpha = LT_get_col(ST_td_tbl, i_alpha_eff, ST_refine_adx_fac * fld) * &
            gas_number_density / ST_refine_adx_fac
       adx   = box%dr * alpha
       elec_dens = box%cc(IJK, i_electron)

       if (adx  > ST_refine_adx .and. elec_dens > ST_refine_min_dens) then
          cell_flags(IJK) = af_do_ref
       else if (adx < 0.125_dp * ST_refine_adx .and. dx < ST_derefine_dx) then
          cell_flags(IJK) = af_rm_ref
       else
          cell_flags(IJK) = af_keep_ref
       end if

       ! Refine around the initial conditions
       if (ST_time < ST_refine_init_time) then
          do n = 1, init_conds%n_cond
             dist = GM_dist_line(af_r_cc(box, [IJK]), &
                  init_conds%seed_r0(:, n), &
                  init_conds%seed_r1(:, n), NDIM)
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

    ! Make sure we don't have or get a too fine or too coarse grid
    if (dx > ST_refine_max_dx) then
       cell_flags = af_do_ref
    else if (dx < 2 * ST_refine_min_dx) then
       where(cell_flags == af_do_ref) cell_flags = af_keep_ref
    end if

  end subroutine refine_routine

  subroutine forward_euler(tree, dt, s_in, s_out)
    use m_chemistry
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt
    integer, intent(in)       :: s_in
    integer, intent(in)       :: s_out
    integer                   :: lvl, i, id, p_id

    ! First calculate fluxes
    !$omp parallel private(lvl, i, id)
    do lvl = tree%lowest_lvl, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call fluxes_elec(tree%boxes, id, dt, s_in)
       end do
       !$omp end do
    end do
    !$omp end parallel

    call af_consistent_fluxes(tree, [flux_elec])

    ! Update the solution
    !$omp parallel private(lvl, i, id, p_id)
    do lvl = tree%lowest_lvl, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call update_solution(tree%boxes(id), dt, s_in, s_out)

          ! This can be important for setting ghost cells
          p_id = tree%boxes(id)%parent
          if (p_id > af_no_box) then
             call af_restrict_box_vars(tree%boxes(id), tree%boxes(p_id), &
                  flux_species)
          end if
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine forward_euler

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

  !> Compute the electron fluxes due to drift and diffusion
  subroutine fluxes_elec(boxes, id, dt, s_in)
    use m_flux_schemes
    use m_units_constants
    type(box_t), intent(inout) :: boxes(:)
    integer, intent(in)        :: id
    real(dp), intent(in)       :: dt
    integer, intent(in)        :: s_in
    real(dp)                   :: inv_dr, fld, Td
    ! Velocities at cell faces
    real(dp), allocatable      :: v(DTIMES(:), :)
    ! Diffusion coefficients at cell faces
    real(dp), allocatable      :: dc(DTIMES(:), :)
    ! Maximal fluxes at cell faces
    real(dp), allocatable      :: fmax(DTIMES(:), :)
    ! Cell-centered densities
    real(dp), allocatable      :: cc(DTIMES(:))
    real(dp)                   :: mu, fld_face, drt_fac, tmp
    real(dp)                   :: nsmall, N_inv
    integer                    :: nc, n, m
#if NDIM == 3
    integer                    :: l
#endif

    nc      = boxes(id)%n_cell
    inv_dr  = 1/boxes(id)%dr
    drt_fac = UC_eps0 / max(1e-100_dp, UC_elem_charge * dt)
    nsmall  = 1.0_dp ! A small density
    N_inv   = 1/gas_number_density

    allocate(v(DTIMES(1:nc+1), NDIM))
    allocate(dc(DTIMES(1:nc+1), NDIM))
    allocate(cc(DTIMES(-1:nc+2)))
    allocate(fmax(DTIMES(1:nc+1), NDIM))
    if (ST_drt_limit_flux) then
       fmax = 0.0_dp
    end if

    ! Fill ghost cells on the sides of boxes (no corners)
    call af_gc_box(boxes, id, i_electron+s_in, af_gc_interp_lim, &
         af_bc_neumann_zero, .false.)

    ! Fill cc with interior data plus a second layer of ghost cells
    call af_gc2_box(boxes, id, i_electron+s_in, af_gc2_prolong_linear, &
         af_bc2_neumann_zero, cc, nc)

    ! We use the average field to compute the mobility and diffusion coefficient
    ! at the interface
    do n = 1, nc+1
       do m = 1, nc
#if NDIM == 2
          fld = 0.5_dp * (boxes(id)%cc(n-1, m, i_electric_fld) + &
               boxes(id)%cc(n, m, i_electric_fld))
          Td = fld * SI_to_Townsend
          mu = LT_get_col(ST_td_tbl, i_mobility, Td) * N_inv
          fld_face = boxes(id)%fc(n, m, 1, electric_fld)
          v(n, m, 1)  = mu * fld_face
          dc(n, m, 1) = LT_get_col(ST_td_tbl, i_diffusion, Td) * N_inv

          if (ST_drt_limit_flux) then
             tmp = abs(cc(n-1, m) - cc(n, m))/max(cc(n-1, m), cc(n, m), nsmall)
             tmp = max(fld, tmp * inv_dr * dc(n, m, 1) / mu)
             fmax(n, m, 1) = drt_fac * tmp
          end if

          if (abs(fld_face) > ST_diffusion_field_limit .and. &
               (cc(n-1, m) - cc(n, m)) * fld_face > 0.0_dp) then
             dc(n, m, 1) = 0.0_dp
          end if

          fld = 0.5_dp * (boxes(id)%cc(m, n-1, i_electric_fld) + &
               boxes(id)%cc(m, n, i_electric_fld))
          Td = fld * SI_to_Townsend
          mu = LT_get_col(ST_td_tbl, i_mobility, Td) * N_inv
          fld_face = boxes(id)%fc(m, n, 2, electric_fld)
          v(m, n, 2)  = -mu * fld_face
          dc(m, n, 2) = LT_get_col(ST_td_tbl, i_diffusion, Td) * N_inv

          if (ST_drt_limit_flux) then
             tmp = abs(cc(m, n-1) - cc(m, n))/max(cc(m, n-1), cc(m, n), nsmall)
             tmp = max(fld, tmp * inv_dr * dc(m, n, 2) / mu)
             fmax(m, n, 2) = drt_fac * tmp
          end if

          if (abs(fld_face) > ST_diffusion_field_limit .and. &
               (cc(m, n-1) - cc(m, n)) * fld_face > 0.0_dp) then
             dc(m, n, 2) = 0.0_dp
          end if
#elif NDIM == 3
          do l = 1, nc
             fld = 0.5_dp * (&
                  boxes(id)%cc(n-1, m, l, i_electric_fld) + &
                  boxes(id)%cc(n, m, l, i_electric_fld))
             Td = fld * SI_to_Townsend
             mu = LT_get_col(ST_td_tbl, i_mobility, Td) * N_inv
             fld_face = boxes(id)%fc(n, m, l, 1, electric_fld)
             v(n, m, l, 1)  = -mu * fld_face
             dc(n, m, l, 1) = LT_get_col(ST_td_tbl, i_diffusion, Td) * N_inv

             if (ST_drt_limit_flux) then
                tmp = abs(cc(n-1, m, l) - cc(n, m, l)) / &
                     max(cc(n-1, m, l), cc(n, m, l), nsmall)
                tmp = max(fld, tmp * inv_dr * dc(n, m, l, 1) / mu)
                fmax(n, m, l, 1) = drt_fac * tmp
             end if

             if (abs(fld_face) > ST_diffusion_field_limit .and. &
                  (cc(n-1, m, l) - cc(n, m, l)) * fld_face > 0.0_dp) then
                dc(n, m, l, 1) = 0.0_dp
             end if

             fld = 0.5_dp * (&
                  boxes(id)%cc(m, n-1, l, i_electric_fld) + &
                  boxes(id)%cc(m, n, l, i_electric_fld))
             Td = fld * SI_to_Townsend
             mu = LT_get_col(ST_td_tbl, i_mobility, Td) * N_inv
             fld_face = boxes(id)%fc(m, n, l, 2, electric_fld)
             v(m, n, l, 2)  = -mu * fld_face
             dc(m, n, l, 2) = LT_get_col(ST_td_tbl, i_diffusion, Td) * N_inv

             if (ST_drt_limit_flux) then
                tmp = abs(cc(m, n-1, l) - cc(m, n, l)) / &
                     max(cc(m, n-1, l), cc(m, n, l), nsmall)
                tmp = max(fld, tmp * inv_dr * dc(m, n, l, 2) / mu)
                fmax(m, n, l, 2) = drt_fac * tmp
             end if

             if (abs(fld_face) > ST_diffusion_field_limit .and. &
                  (cc(m, n-1, l) - cc(m, n, l)) * fld_face > 0.0_dp) then
                dc(m, n, l, 2) = 0.0_dp
             end if

             fld = 0.5_dp * (&
                  boxes(id)%cc(m, l, n-1, i_electric_fld) + &
                  boxes(id)%cc(m, l, n, i_electric_fld))
             Td = fld * SI_to_Townsend
             mu = LT_get_col(ST_td_tbl, i_mobility, Td) * N_inv
             fld_face = boxes(id)%fc(m, l, n, 3, electric_fld)
             v(m, l, n, 3)  = -mu * fld_face
             dc(m, l, n, 3) = LT_get_col(ST_td_tbl, i_diffusion, Td) * N_inv

             if (ST_drt_limit_flux) then
                tmp = abs(cc(m, l, n-1) - cc(m, l, n)) / &
                     max(cc(m, l, n-1), cc(m, l, n), nsmall)
                tmp = max(fld, tmp * inv_dr * dc(m, l, n, 3) / mu)
                fmax(m, l, n, 3) = drt_fac * tmp
             end if

             if (abs(fld_face) > ST_diffusion_field_limit .and. &
                  (cc(m, l, n-1) - cc(m, l, n)) * fld_face > 0.0_dp) then
                dc(m, l, n, 3) = 0.0_dp
             end if
          end do
#endif
       end do
    end do

    if (ST_max_velocity > 0.0_dp) then
       where (abs(v) > ST_max_velocity)
          v = sign(ST_max_velocity, v)
       end where
    end if

    ! if (set_dt) then
    !    ! CFL condition
    !    if (ST_max_velocity > 0) then
    !       tmp_vec = abs(fld_vec * max(mu, abs(ST_ion_mobility)))
    !       where (tmp_vec > ST_max_velocity) tmp_vec = ST_max_velocity
    !       dt_cfl = 1.0_dp/sum(tmp_vec * inv_dr)
    !    else
    !       dt_cfl = 1.0_dp/sum(abs(fld_vec * &
    !            max(mu, abs(ST_ion_mobility))) * inv_dr)
    !    end if

    !    ! Diffusion condition
    !    dt_dif = box%dr**2 / max(2 * NDIM * max(diff, ST_ion_diffusion), &
    !         epsilon(1.0_dp))

    !    ! Dielectric relaxation time
    !    if (ST_drt_limit_flux .and. fld < ST_drt_limit_flux_Emax) then
    !       dt_drt = UC_eps0 / max(UC_elem_charge * abs(ST_ion_mobility) * &
    !            box%cc(IJK, i_electron), epsilon(1.0_dp))
    !    else
    !       dt_drt = UC_eps0 / (UC_elem_charge * (mu + abs(ST_ion_mobility)) * &
    !            max(box%cc(IJK, i_electron), epsilon(1.0_dp)))
    !    end if

    !    ST_dt_matrix(ST_ix_drt, tid) = min(ST_dt_matrix(ST_ix_drt, tid), &
    !         dt_drt)

    !    ! Take the combined CFL-diffusion condition with Courant number 0.5.
    !    ! The 0.5 is emperical, to have good accuracy (and TVD/positivity) in
    !    ! combination with the explicit trapezoidal rule
    !    dt_cfl = 0.5_dp/(1/dt_cfl + 1/dt_dif)

    !    ST_dt_matrix(ST_ix_cfl, tid) = min(ST_dt_matrix(ST_ix_cfl, tid), &
    !         dt_cfl)
    !    ST_dt_matrix(ST_ix_diff, tid) = min(ST_dt_matrix(ST_ix_diff, tid), &
    !         dt_dif)
    ! end if

#if NDIM == 2
    call flux_koren_2d(cc, v, nc, 2)
    call flux_diff_2d(cc, dc, inv_dr, nc, 2)
#elif NDIM == 3
    call flux_koren_3d(cc, v, nc, 2)
    call flux_diff_3d(cc, dc, inv_dr, nc, 2)
#endif

    boxes(id)%fc(DTIMES(:), :, flux_elec) = v + dc

    if (ST_drt_limit_flux) then
       where (abs(boxes(id)%fc(DTIMES(:), :, flux_elec)) > fmax)
          boxes(id)%fc(DTIMES(:), :, flux_elec) = &
               sign(fmax, boxes(id)%fc(DTIMES(:), :, flux_elec))
       end where
    end if

    !     if (ST_update_ions) then
    !        ! Use a constant diffusion coefficient for ions
    !        dc = ST_ion_diffusion

    !        ! Use a constant mobility for ions
    !        v(DTIMES(:), 1:NDIM) = ST_ion_mobility * &
    !             boxes(id)%fc(DTIMES(:), 1:NDIM, electric_fld)

    !        ! Fill ghost cells on the sides of boxes (no corners)
    !        call af_gc_box(boxes, id, i_pos_ion, af_gc_interp_lim, &
    !             af_bc_neumann_zero, .false.)

    !        ! Fill cc with interior data plus a second layer of ghost cells
    !        call af_gc2_box(boxes, id, i_pos_ion, af_gc2_prolong_linear, &
    !             af_bc2_neumann_zero, cc, nc)

    ! #if NDIM == 2
    !        call flux_koren_2d(cc, v, nc, 2)
    !        call flux_diff_2d(cc, dc, inv_dr, nc, 2)
    ! #elif NDIM == 3
    !        call flux_koren_3d(cc, v, nc, 2)
    !        call flux_diff_3d(cc, dc, inv_dr, nc, 2)
    ! #endif

    !        boxes(id)%fc(DTIMES(:), :, flux_ion) = v + dc
    !     end if
  end subroutine fluxes_elec

  !> Advance solution over dt based on the fluxes / source term, using forward Euler
  subroutine update_solution(box, dt, s_in, s_out)
    use omp_lib
    use m_units_constants
    use m_gas
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: dt
    integer, intent(in)        :: s_in
    integer, intent(in)        :: s_out
    real(dp)                   :: inv_dr
    real(dp), allocatable      :: rates(:, :)
    real(dp), allocatable      :: derivs(:, :)
    real(dp), allocatable      :: dens(:, :)
    real(dp), allocatable      :: fields(:)
#if NDIM == 2
    real(dp)                   :: rfac(2)
    integer                    :: ioff
#endif
    integer                    :: IJK, ix, nc, n_cells, n, iv

    nc     = box%n_cell
    n_cells = box%n_cell**NDIM
    inv_dr = 1/box%dr
#if NDIM == 2
    ioff   = (box%ix(1)-1) * nc
#endif

    allocate(rates(n_cells, n_reactions))
    allocate(derivs(n_cells, n_species))

    fields = SI_to_Townsend * &
         reshape(box%cc(DTIMES(1:nc), i_electric_fld), [n_cells])

    dens = reshape(box%cc(DTIMES(1:nc), species_ix(1:n_species)+s_in), &
         [n_cells, n_species])

    call get_rates(fields, rates, n_cells)
    call get_derivatives(dens, rates, derivs, n_cells)

    ix = 0
    do KJI_DO(1,nc)
       ix = ix + 1
       ! Contribution of flux
#if NDIM == 2
       if (ST_cylindrical) then
          ! Weighting of flux contribution for cylindrical coordinates
          rfac(:) = [i+ioff-1, i+ioff] / (i+ioff-0.5_dp)
       else
          rfac(:) = 1.0_dp
       end if

       derivs(ix, i_electron) = derivs(ix, i_electron) + inv_dr * ( &
            box%fc(i, j, 2, flux_elec) - &
            box%fc(i, j+1, 2, flux_elec) + &
            rfac(1) * box%fc(i, j, 1, flux_elec) - &
            rfac(2) * box%fc(i+1, j, 1, flux_elec))
#elif NDIM == 3
       derivs(ix, i_electron) = derivs(ix, i_electron) + inv_dr * ( &
            sum(box%fc(i, j, k, 1:3, flux_elec)) - &
            box%fc(i+1, j, k, 1, flux_elec) - &
            box%fc(i, j+1, k, 2, flux_elec) - &
            box%fc(i, j, k+1, 3, flux_elec))
#endif

       if (photoi_enabled) then
          derivs(ix, i_electron) = derivs(ix, i_electron) + &
               box%cc(IJK, i_photo)
          derivs(ix, i_1pos_ion) = derivs(ix, i_1pos_ion) + &
               box%cc(IJK, i_photo)
       end if

       do n = 1, n_species
          iv = species_ix(n)
          box%cc(IJK, iv+s_out) = box%cc(IJK, iv+s_in) + dt * derivs(ix, n)
       end do
    end do; CLOSE_DO

  end subroutine update_solution

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

    dt = ST_dt_safety_factor * minval(ST_dt_matrix(1:ST_dt_num_cond, :))

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
