!> Program that solves a 2d streamer
!
! Author: Jannis Teunissen
! License: GPLv3
program streamer_2d

  use m_a2_types
  use m_a2_core
  use m_a2_ghostcell
  use m_a2_utils
  use m_a2_restrict
  use m_a2_multigrid
  use m_a2_output
  use m_write_silo
  use m_streamer

  implicit none

  integer                :: i, n
  character(len=ST_slen) :: fname
  logical                :: write_out

  type(a2_t)             :: tree ! This contains the full grid information
  type(mg2_t)            :: mg   ! Multigrid option struct
  type(ref_info_t)       :: ref_info

  call ST_create_config()
  call ST_read_config_files()
  call ST_load_config()

  ! Initialize the transport coefficients
  call ST_load_transport_data()

  ! Set the initial conditions from the configuration
  call ST_get_init_cond(2)

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the multigrid options. First define the variables to use
  mg%i_phi        = i_phi
  mg%i_tmp        = i_electric_fld
  mg%i_rhs        = i_rhs

  ! The number of cycles at the lowest level
  mg%n_cycle_base = 8

  ! Routines to use for ...
  mg%sides_bc    => sides_bc_potential ! Filling ghost cell on physical boundaries
  mg%box_op      => mg2_auto_op
  mg%box_corr    => mg2_auto_corr
  mg%box_gsrb    => mg2_auto_gsrb

  ! This routine always needs to be called when using multigrid
  call mg2_init_mg(mg)

  ST_out_cnt = 0 ! Number of output files written
  ST_time    = 0 ! Simulation time (all times are in s)

  ! Set up the initial conditions
  do
     call a2_loop_box(tree, set_initial_condition)
     call compute_electric_field(tree, n_fmg_cycles, .false.)
     call a2_adjust_refinement(tree, refine_routine, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  call a2_print_info(tree)

  if (ST_photoi_enabled) &
       call set_photoionization(tree, ST_photoi_eta, ST_photoi_num_photons)

  do
     ! Get a new time step, which is at most dt_amr
     ST_dt = get_max_dt(tree)

     if (ST_dt < 1e-14) then
        print *, "ST_dt getting too small, instability?"
        ST_time = ST_end_time + 1.0_dp
     end if

     ! Every ST_dt_out, write output
     if (ST_out_cnt * ST_dt_out <= ST_time) then
        write_out = .true.
        ST_out_cnt = ST_out_cnt + 1
        write(fname, "(A,I6.6)") trim(ST_simulation_name) // "_", ST_out_cnt
     else
        write_out = .false.
     end if

     if (ST_time > ST_end_time) exit

     ! We perform n_steps between mesh-refinements
     do n = 1, ST_refine_per_steps

        if (ST_photoi_enabled) &
             call set_photoionization(tree, ST_photoi_eta, ST_photoi_num_photons, ST_dt)

        ! Copy previous solution
        call a2_tree_copy_cc(tree, i_electron, i_electron_old)
        call a2_tree_copy_cc(tree, i_pos_ion, i_pos_ion_old)

        ! Two forward Euler steps over ST_dt
        do i = 1, 2
           ST_time = ST_time + ST_dt

           ! First calculate fluxes
           call a2_loop_boxes_arg(tree, fluxes_koren, [ST_dt], .true.)
           call a2_consistent_fluxes(tree, [flux_elec])

           ! Update the solution
           call a2_loop_box_arg(tree, update_solution, [ST_dt], .true.)

           ! Restrict the electron and ion densities to lower levels
           call a2_restrict_tree(tree, i_electron)
           call a2_restrict_tree(tree, i_pos_ion)

           ! Fill ghost cells
           call a2_gc_tree(tree, i_electron, a2_gc_interp_lim, a2_bc_neumann_zero)
           call a2_gc_tree(tree, i_pos_ion, a2_gc_interp_lim, a2_bc_neumann_zero)

           ! Compute new field on first iteration
           if (i == 1) call compute_electric_field(tree, n_fmg_cycles, .true.)
        end do

        ST_time = ST_time - ST_dt        ! Go back one time step

        ! Take average of phi_old and phi (explicit trapezoidal rule)
        call a2_loop_box(tree, average_density)

        ! Compute field with new density
        call compute_electric_field(tree, n_fmg_cycles, .true.)
     end do

     if (write_out) call a2_write_silo(tree, fname, ST_out_cnt, &
          ST_time, dir=ST_output_dir)

     call a2_adjust_refinement(tree, refine_routine, ref_info)

     if (ref_info%n_add > 0 .or. ref_info%n_rm > 0) then
        ! For boxes which just have been refined, set data on their children
        call prolong_to_new_boxes(tree, ref_info)

        ! Compute the field on the new mesh
        call compute_electric_field(tree, n_fmg_cycles, .true.)
     end if

  end do

  call a2_destroy(tree)

contains

  !> Initialize the AMR tree
  subroutine init_tree(tree)
    type(a2_t), intent(inout) :: tree

    ! Variables used below to initialize tree
    real(dp)                  :: dr
    integer                   :: id
    integer                   :: ix_list(2, 1) ! Spatial indices of initial boxes
    integer                   :: nb_list(4, 1) ! Neighbors of initial boxes
    integer                   :: n_boxes_init = 1000

    dr = ST_domain_len / ST_box_size

    ! Initialize tree
    call a2_init(tree, ST_box_size, n_var_cell, n_var_face, dr, &
         coarsen_to=2, n_boxes=n_boxes_init, cc_names=ST_cc_names)

    ! Set up geometry
    id             = 1          ! One box ...
    ix_list(:, id) = [1,1]      ! With index 1,1 ...
    nb_list(:, id) = -1         ! And neighbors -1 (physical boundary)

    ! Create the base mesh
    call a2_set_base(tree, ix_list, nb_list)

  end subroutine init_tree

  ! This routine sets the refinement flag for boxes(id)
  subroutine refine_routine(boxes, id, refine_flag)
    use m_geometry
    type(box2_t), intent(in) :: boxes(:) ! List of all boxes
    integer, intent(in)      :: id       ! Index of box to look at
    integer, intent(inout)   :: refine_flag ! Refinement flag for the box
    integer                  :: n, nc
    real(dp)                 :: cphi, dx, dx2
    real(dp)                 :: alpha, adx, max_fld
    real(dp)                 :: boxlen, dist

    nc      = boxes(id)%n_cell
    dx      = boxes(id)%dr
    dx2     = boxes(id)%dr**2

    max_fld = maxval(boxes(id)%cc(1:nc, 1:nc, i_electric_fld))
    alpha   = LT_get_col(ST_td_tbl, i_alpha, max_fld)

    ! The refinement is based on the ionization length times the grid spacing
    adx     = boxes(id)%dr * alpha

    ! Increase the refinement in fields below ST_refine_adx_field
    adx     = adx * max(1.0_dp, (2 * ST_refine_adx_field) / &
         (max_fld + ST_refine_adx_field))

    ! The refinement is also based on the intensity of the source term. Here we
    ! estimate the curvature of phi (given by dx**2 * Laplacian(phi))
    cphi    = dx2 * maxval(abs(boxes(id)%cc(1:nc, 1:nc, i_rhs)))

    if (adx / ST_refine_adx + cphi / ST_refine_cphi > 1) then
       refine_flag = af_do_ref
    else if (adx < 0.125_dp * ST_refine_adx .and. cphi < 0.0625_dp * ST_refine_cphi &
         .and. dx < ST_derefine_dx) then
       refine_flag = af_rm_ref
    end if

    ! Refine around the initial conditions
    if (ST_time < ST_refine_init_time) then
       boxlen = boxes(id)%n_cell * boxes(id)%dr

       do n = 1, ST_init_cond%n_cond
          dist = GM_dist_line(a2_r_center(boxes(id)), &
               ST_init_cond%seed_r0(:, n), &
               ST_init_cond%seed_r1(:, n), 2)
          if (dist - 2 * ST_init_cond%seed_width(n) < boxlen &
               .and. boxes(id)%dr > ST_refine_init_fac * &
               ST_init_cond%seed_width(n)) then
             refine_flag = af_do_ref
          end if
       end do
    end if

    ! Make sure we don't have or get a too fine or too coarse grid
    if (dx > ST_refine_max_dx) then
       refine_flag = af_do_ref
    else if (dx < ST_refine_min_dx) then
       refine_flag = af_rm_ref
    else if (dx < 2 * ST_refine_min_dx .and. refine_flag == af_do_ref) then
       refine_flag = af_keep_ref
    else if (dx > 0.5_dp * ST_refine_max_dx .and. refine_flag == af_rm_ref) then
       refine_flag = af_keep_ref
    end if

  end subroutine refine_routine

  !> Sets the initial condition
  subroutine set_initial_condition(box)
    use m_geometry
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, n, nc
    real(dp)                    :: xy(2)
    real(dp)                    :: density

    nc = box%n_cell
    box%cc(:, :, i_electron) = ST_init_cond%bg_dens

    do j = 0, nc+1
       do i = 0, nc+1
          xy   = a2_r_cc(box, [i,j])

          do n = 1, ST_init_cond%n_cond
             density = ST_init_cond%seed_density(n) * &
                       GM_density_line(xy, ST_init_cond%seed_r0(:, n), &
                       ST_init_cond%seed_r1(:, n), 2, &
                       ST_init_cond%seed_width(n), &
                       ST_init_cond%seed_falloff(n))
             box%cc(i, j, i_electron) = box%cc(i, j, i_electron) + density
          end do
       end do
    end do

    box%cc(:, :, i_pos_ion) = box%cc(:, :, i_electron)
    box%cc(:, :, i_phi)     = 0     ! Inital potential set to zero

  end subroutine set_initial_condition

  !> Get maximum time step based on e.g. CFL criteria
  real(dp) function get_max_dt(tree)
    type(a2_t), intent(in) :: tree
    real(dp), parameter    :: UC_eps0        = 8.8541878176d-12
    real(dp), parameter    :: UC_elem_charge = 1.6022d-19
    real(dp)               :: max_fld, min_fld, max_dns, dr_min
    real(dp)               :: mobility, diff_coeff, alpha, max_mobility
    real(dp)               :: dt_cfl, dt_dif, dt_drt, dt_alpha

    call a2_tree_max_cc(tree, i_electric_fld, max_fld)
    call a2_tree_min_cc(tree, i_electric_fld, min_fld)
    call a2_tree_max_cc(tree, i_electron, max_dns)

    dr_min       = a2_min_dr(tree)
    mobility     = LT_get_col(ST_td_tbl, i_mobility, max_fld)
    max_mobility = LT_get_col(ST_td_tbl, i_mobility, min_fld)
    diff_coeff   = LT_get_col(ST_td_tbl, i_diffusion, max_fld)
    alpha        = LT_get_col(ST_td_tbl, i_alpha, max_fld)

    ! CFL condition
    dt_cfl = dr_min / (mobility * max_fld) ! Factor ~ sqrt(0.5)

    ! Diffusion condition
    dt_dif = 0.25_dp * dr_min**2 / diff_coeff

    ! Dielectric relaxation time
    dt_drt = UC_eps0 / (UC_elem_charge * max_mobility * max_dns)

    ! Ionization limit
    dt_alpha =  1 / max(mobility * max_fld * alpha, epsilon(1.0_dp))

    get_max_dt = 0.5_dp * min(1/(1/dt_cfl + 1/dt_dif), dt_alpha, ST_dt_max)
  end function get_max_dt

  !> Compute electric field on the tree. First perform multigrid to get electric
  ! potential, then take numerical gradient to geld field.
  subroutine compute_electric_field(tree, n_cycles, have_guess)
    use m_units_constants
    type(a2_t), intent(inout) :: tree
    integer, intent(in)       :: n_cycles
    logical, intent(in)       :: have_guess
    real(dp), parameter       :: fac = UC_elem_charge / UC_eps0
    integer                   :: lvl, i, id, nc

    nc = tree%n_cell

    ! Set the source term (rhs)
    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          tree%boxes(id)%cc(:, :, i_rhs) = fac * (&
               tree%boxes(id)%cc(:, :, i_electron) - &
               tree%boxes(id)%cc(:, :, i_pos_ion))
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    call ST_set_voltage(ST_time)

    ! Perform n_cycles fmg cycles (logicals: store residual, first call)
    do i = 1, n_cycles
       call mg2_fas_fmg(tree, mg, .false., have_guess .or. i > 1)
    end do

    ! Compute field from potential
    call a2_loop_box(tree, electric_field_from_potential)

    ! Set the field norm also in ghost cells
    call a2_gc_tree(tree, i_electric_fld, a2_gc_interp, a2_bc_neumann_zero)
  end subroutine compute_electric_field

  !> Compute electric field from electrical potential
  subroutine electric_field_from_potential(box)
    type(box2_t), intent(inout) :: box
    integer                     :: nc
    real(dp)                    :: inv_dr

    nc     = box%n_cell
    inv_dr = 1 / box%dr

    box%fx(:, :, electric_fld) = inv_dr * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi))
    box%fy(:, :, electric_fld) = inv_dr * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi))

    box%cc(1:nc, 1:nc, i_electric_fld) = 0.5_dp * sqrt(&
         (box%fx(1:nc, 1:nc, electric_fld) + box%fx(2:nc+1, 1:nc, electric_fld))**2 + &
         (box%fy(1:nc, 1:nc, electric_fld) + box%fy(1:nc, 2:nc+1, electric_fld))**2)
  end subroutine electric_field_from_potential

  !> This fills ghost cells near physical boundaries for the potential
  subroutine sides_bc_potential(box, nb, iv, bc_type)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: nb ! Direction for the boundary condition
    integer, intent(in)         :: iv ! Index of variable
    integer, intent(out)        :: bc_type ! Type of boundary condition
    integer                     :: nc

    nc = box%n_cell

    select case (nb)
    case (a2_neighb_lowx)
       bc_type = af_bc_neumann
       box%cc(   0, 1:nc, iv) = 0
    case (a2_neighb_highx)
       bc_type = af_bc_neumann
       box%cc(nc+1, 1:nc, iv) = 0
    case (a2_neighb_lowy)
       bc_type = af_bc_dirichlet
       box%cc(1:nc,    0, iv) = 0
    case (a2_neighb_highy)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, nc+1, iv) = ST_applied_voltage
    end select
  end subroutine sides_bc_potential

  !> Compute the electron fluxes due to drift and diffusion
  subroutine fluxes_koren(boxes, id, dt_vec)
    use m_units_constants
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp), intent(in)        :: dt_vec(:)
    real(dp)                    :: inv_dr, tmp, gradp, gradc, gradn
    real(dp)                    :: mobility, diff_coeff, v_drift
    real(dp)                    :: fld_avg, fld
    real(dp)                    :: gc_data(boxes(id)%n_cell, a2_num_neighbors)
    integer                     :: i, j, nc
    type(LT_loc_t)              :: loc

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    call a2_gc2_box(boxes, id, i_electron, a2_gc2_prolong_linear, &
         a2_bc2_neumann_zero, gc_data, nc)

    ! x-fluxes interior, advective part with flux limiter
    do j = 1, nc
       do i = 1, nc+1
          fld_avg   = 0.5_dp * (boxes(id)%cc(i, j, i_electric_fld) + &
               boxes(id)%cc(i-1, j, i_electric_fld))
          loc        = LT_get_loc(ST_td_tbl, fld_avg)
          mobility   = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)
          diff_coeff = LT_get_col_at_loc(ST_td_tbl, i_diffusion, loc)
          fld        = boxes(id)%fx(i, j, electric_fld)
          v_drift    = -mobility * fld
          gradc      = boxes(id)%cc(i  , j, i_electron) - &
                       boxes(id)%cc(i-1, j, i_electron)

          if (v_drift < 0.0_dp) then
             if (i == nc+1) then
                tmp = gc_data(j, a2_neighb_highx)
             else
                tmp = boxes(id)%cc(i+1, j, i_electron)
             end if
             gradn = tmp - boxes(id)%cc(i, j, i_electron)
             boxes(id)%fx(i, j, flux_elec) = v_drift * &
                  (boxes(id)%cc(i, j, i_electron) - koren_mlim(gradc, gradn))
          else                  ! v_drift > 0
             if (i == 1) then
                tmp = gc_data(j, a2_neighb_lowx)
             else
                tmp = boxes(id)%cc(i-2, j, i_electron)
             end if
             gradp = boxes(id)%cc(i-1, j, i_electron) - tmp
             boxes(id)%fx(i, j, flux_elec) = v_drift * &
                  (boxes(id)%cc(i-1, j, i_electron) + koren_mlim(gradc, gradp))
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be
          ! scaled by 1/dx
          boxes(id)%fx(i, j, flux_elec) = boxes(id)%fx(i, j, flux_elec) - &
               diff_coeff * gradc * inv_dr
       end do
    end do

    ! y-fluxes interior, advective part with flux limiter
    do j = 1, nc+1
       do i = 1, nc
          fld_avg    = 0.5_dp * (boxes(id)%cc(i, j, i_electric_fld) + &
               boxes(id)%cc(i, j-1, i_electric_fld))
          loc        = LT_get_loc(ST_td_tbl, fld_avg)
          mobility   = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)
          diff_coeff = LT_get_col_at_loc(ST_td_tbl, i_diffusion, loc)
          fld        = boxes(id)%fy(i, j, electric_fld)
          v_drift    = -mobility * fld
          gradc      = boxes(id)%cc(i, j  , i_electron) - &
                       boxes(id)%cc(i, j-1, i_electron)

          if (v_drift < 0.0_dp) then
             if (j == nc+1) then
                tmp = gc_data(i, a2_neighb_highy)
             else
                tmp = boxes(id)%cc(i, j+1, i_electron)
             end if
             gradn = tmp - boxes(id)%cc(i, j, i_electron)
             boxes(id)%fy(i, j, flux_elec) = v_drift * &
                  (boxes(id)%cc(i, j, i_electron) - koren_mlim(gradc, gradn))
          else                  ! v_drift > 0
             if (j == 1) then
                tmp = gc_data(i, a2_neighb_lowy)
             else
                tmp = boxes(id)%cc(i, j-2, i_electron)
             end if
             gradp = boxes(id)%cc(i, j-1, i_electron) - tmp
             boxes(id)%fy(i, j, flux_elec) = v_drift * &
                  (boxes(id)%cc(i, j-1, i_electron) + koren_mlim(gradc, gradp))
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be
          ! scaled by 1/dx
          boxes(id)%fy(i, j, flux_elec) = boxes(id)%fy(i, j, flux_elec) - &
               diff_coeff * gradc * inv_dr
       end do
    end do

  end subroutine fluxes_koren

  !> Take average of new and old electron/ion density for explicit trapezoidal rule
  subroutine average_density(box)
    type(box2_t), intent(inout) :: box
    box%cc(:, :, i_electron) = 0.5_dp *  &
         (box%cc(:, :, i_electron) + box%cc(:, :, i_electron_old))
    box%cc(:, :, i_pos_ion) = 0.5_dp * &
         (box%cc(:, :, i_pos_ion)  + box%cc(:, :, i_pos_ion_old))
  end subroutine average_density

  !> Advance solution over dt based on the fluxes / source term, using forward Euler
  subroutine update_solution(box, dt)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr, src, fld
    real(dp)                    :: alpha, eta, sflux, mu
    integer                     :: i, j, nc, ioff
    type(LT_loc_t)              :: loc

    nc     = box%n_cell
    inv_dr = 1/box%dr
    ioff   = (box%ix(1)-1) * nc

    do j = 1, nc
       do i = 1, nc
          ! Weighting of flux contribution for cylindrical coordinates
          fld   = box%cc(i,j, i_electric_fld)
          loc   = LT_get_loc(ST_td_tbl, fld)
          alpha = LT_get_col_at_loc(ST_td_tbl, i_alpha, loc)
          eta   = LT_get_col_at_loc(ST_td_tbl, i_eta, loc)
          mu    = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)

          ! Contribution of flux
          sflux = (box%fy(i, j, flux_elec) - box%fy(i, j+1, flux_elec) + &
                   box%fx(i, j, flux_elec) - box%fx(i+1, j, flux_elec)) * inv_dr * dt(1)

          ! Source term
          src = fld * mu * box%cc(i, j, i_electron) * (alpha - eta) * dt(1)

          if (ST_photoi_enabled) src = src + box%cc(i,j, i_photo) * dt(1)

          ! Add flux and source term
          box%cc(i, j, i_electron) = box%cc(i, j, i_electron) + sflux + src

          ! Add source term
          box%cc(i, j, i_pos_ion)  = box%cc(i, j, i_pos_ion) + src

       end do
    end do
  end subroutine update_solution

  !> Sets the photoionization
  subroutine set_photoionization(tree, eta, num_photons, dt)
    use m_units_constants

    type(a2_t), intent(inout) :: tree
    real(dp), intent(in)      :: eta
    real(dp), intent(in), optional :: dt
    integer, intent(in)       :: num_photons
    real(dp), parameter       :: p_quench = 30.0D0 * UC_torr_to_bar
    real(dp)                  :: quench_fac

    ! Compute quench factor, because some excited species will be quenched by
    ! collisions, preventing the emission of a UV photon
    quench_fac = p_quench / (ST_gas_pressure + p_quench)

    ! Set photon production rate per cell, which is proportional to the
    ! ionization rate.
    call a2_loop_box_arg(tree, set_photoionization_rate, [eta * quench_fac], .true.)

    call photoi_set_src_2d(tree, ST_photoi_tbl, ST_rng, num_photons, &
         i_photo, i_photo, 0.25e-3_dp, .false., .false., 1e-9_dp, dt)

  end subroutine set_photoionization

  !> Sets the photoionization_rate
  subroutine set_photoionization_rate(box, coeff)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: coeff(:)
    integer                     :: i, j, nc
    real(dp)                    :: fld, alpha, mobility, dr, tmp
    type(LT_loc_t)              :: loc

    nc = box%n_cell

    do j = 1, nc
       do i = 1, nc
          dr       = box%dr
          fld      = box%cc(i, j, i_electric_fld)
          loc      = LT_get_loc(ST_td_tbl, fld)
          alpha    = LT_get_col_at_loc(ST_td_tbl, i_alpha, loc)
          mobility = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)

          tmp = fld * mobility * alpha * box%cc(i, j, i_electron) * coeff(1)
          if (tmp < 0) tmp = 0
          box%cc(i, j, i_photo) = tmp
       end do
    end do
  end subroutine set_photoionization_rate

  !> For each box that gets refined, set data on its children using this routine
  subroutine prolong_to_new_boxes(tree, ref_info)
    use m_a2_prolong
    type(a2_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id, p_id

    do lvl = 1, tree%highest_lvl
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          p_id = tree%boxes(id)%parent
          call a2_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_electron)
          call a2_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_pos_ion)
          call a2_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_phi)
       end do

       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a2_gc_box(tree%boxes, id, i_electron, &
               a2_gc_interp_lim, a2_bc_neumann_zero)
          call a2_gc_box(tree%boxes, id, i_pos_ion, &
               a2_gc_interp_lim, a2_bc_neumann_zero)
          call a2_gc_box(tree%boxes, id, i_phi, &
               mg2_sides_rb, mg%sides_bc)
       end do
    end do
  end subroutine prolong_to_new_boxes

end program streamer_2d
