!> Program to perform axisymmetric streamer simulations
program streamer_cyl

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
  mg%n_cycle_base = 4

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
           call a2_loop_boxes(tree, fluxes_koren, .true.)
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
          ST_time, dir=ST_output_dir, ixs_fc=[electric_fld])

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
    integer                   :: n_boxes_init = 1000

    dr = ST_domain_len / ST_box_size

    ! Initialize tree
    call a2_init(tree, ST_box_size, n_var_cell, n_var_face, dr, &
         coarsen_to=2, n_boxes=n_boxes_init, coord=af_cyl, cc_names=ST_cc_names)

    ! Set up geometry
    id             = 1          ! One box ...
    ix_list(:, id) = [1,1]      ! With index 1,1 ...

    ! Create the base mesh
    call a2_set_base(tree, 1, ix_list)

  end subroutine init_tree

  ! This routine sets the cell refinement flags for box
  subroutine refine_routine(box, cell_flags)
    use m_geometry
    type(box2_t), intent(in) :: box
    ! Refinement flags for the cells of the box
    integer, intent(out)     :: cell_flags(box%n_cell, box%n_cell)
    integer                  :: i, j, n, nc
    real(dp)                 :: cphi, dx, dx2
    real(dp)                 :: alpha, adx, fld
    real(dp)                 :: dist

    nc      = box%n_cell
    dx      = box%dr
    dx2     = box%dr**2

    do j = 1, nc
       do i = 1, nc
          fld   = box%cc(i, j, i_electric_fld)
          alpha = LT_get_col(ST_td_tbl, i_alpha, fld)
          ! The refinement is based on the ionization length
          adx   = box%dr * alpha

          ! The refinement is also based on the intensity of the source term.
          ! Here we estimate the curvature of phi (given by dx**2 *
          ! Laplacian(phi))
          cphi = dx2 * abs(box%cc(i, j, i_rhs))

          if (adx / ST_refine_adx + cphi / ST_refine_cphi > 1) then
             cell_flags(i, j) = af_do_ref
          else if (adx < 0.125_dp * ST_refine_adx .and. &
               cphi < 0.0625_dp * ST_refine_cphi &
               .and. dx < ST_derefine_dx) then
             cell_flags(i, j) = af_rm_ref
          else
             cell_flags(i, j) = af_keep_ref
          end if

          ! Refine around the initial conditions
          if (ST_time < ST_refine_init_time) then
             do n = 1, ST_init_cond%n_cond
                dist = GM_dist_line(a2_r_cc(box, [i, j]), &
                     ST_init_cond%seed_r0(:, n), &
                     ST_init_cond%seed_r1(:, n), 2)
                if (dist - ST_init_cond%seed_width(n) < 2 * dx &
                     .and. box%dr > ST_refine_init_fac * &
                     ST_init_cond%seed_width(n)) then
                   cell_flags(i, j) = af_do_ref
                end if
             end do
          end if

       end do
    end do

    ! Make sure we don't have or get a too fine or too coarse grid
    if (dx > ST_refine_max_dx) then
       cell_flags = af_do_ref
    else if (dx < 2 * ST_refine_min_dx) then
       where(cell_flags == af_do_ref) cell_flags = af_keep_ref
    end if

  end subroutine refine_routine

  !> Sets the initial condition
  subroutine set_initial_condition(box)
    use m_geometry
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, n, nc
    real(dp)                    :: xy(2)
    real(dp)                    :: density

    nc                       = box%n_cell
    box%cc(:, :, i_electron) = ST_init_cond%background_density
    box%cc(:, :, i_pos_ion)  = ST_init_cond%background_density
    box%cc(:, :, i_phi)      = 0 ! Inital potential set to zero

    do j = 0, nc+1
       do i = 0, nc+1
          xy   = a2_r_cc(box, [i,j])

          do n = 1, ST_init_cond%n_cond
             density = ST_init_cond%seed_density(n) * &
                  GM_density_line(xy, ST_init_cond%seed_r0(:, n), &
                  ST_init_cond%seed_r1(:, n), 2, &
                  ST_init_cond%seed_width(n), &
                  ST_init_cond%seed_falloff(n))

             ! Add electrons and/or ions depending on the seed charge type
             ! (positive, negative or neutral)
             if (ST_init_cond%seed_charge_type(n) <= 0) then
                box%cc(i, j, i_electron) = box%cc(i, j, i_electron) + density
             end if

             if (ST_init_cond%seed_charge_type(n) >= 0) then
                box%cc(i, j, i_pos_ion) = box%cc(i, j, i_pos_ion) + density
             end if

          end do
       end do
    end do

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
    dt_cfl = sqrt(0.5_dp) * dr_min / (mobility * max_fld)

    ! Diffusion condition
    dt_dif = 0.25_dp * dr_min**2 / diff_coeff

    ! Dielectric relaxation time
    dt_drt = UC_eps0 / (UC_elem_charge * max_mobility * &
         max(epsilon(1.0_dp), max_dns))

    ! Ionization limit
    dt_alpha =  1 / max(mobility * max_fld * alpha, epsilon(1.0_dp))

    get_max_dt = 0.9_dp * min(1/(1/dt_cfl + 1/dt_dif), &
         dt_drt, dt_alpha, ST_dt_max)

    print *, dt_cfl, dt_dif, dt_drt, dt_alpha
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

    box%fc(:, 1:nc, 1, electric_fld) = inv_dr * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi))
    box%fc(1:nc, :, 2, electric_fld) = inv_dr * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi))

    box%cc(1:nc, 1:nc, i_electric_fld) = 0.5_dp * sqrt(&
         (box%fc(1:nc, 1:nc, 1, electric_fld) + &
         box%fc(2:nc+1, 1:nc, 1, electric_fld))**2 + &
         (box%fc(1:nc, 1:nc, 2, electric_fld) + &
         box%fc(1:nc, 2:nc+1, 2, electric_fld))**2)
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
    case (a2_neighb_highx)             ! Neumann
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
  subroutine fluxes_koren(boxes, id)
    use m_units_constants
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp)                    :: inv_dr, gradp, gradc, gradn
    real(dp)                    :: mobility, diff_coeff, v_drift
    real(dp)                    :: fld
    real(dp)                    :: cc(-1:boxes(id)%n_cell+2, -1:boxes(id)%n_cell+2)
    integer                     :: i, j, nc, dim, dix(2)
    type(LT_loc_t)              :: loc

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    call a2_gc2_box(boxes, id, i_electron, a2_gc2_prolong_linear, &
         a2_bc2_neumann_zero, cc, nc)

    do dim = 1, 2
       dix(:) = 0
       dix(dim) = 1

       do j = 1, nc+dix(2)
          do i = 1, nc+dix(1)
             fld        = boxes(id)%fc(i, j, dim, electric_fld)
             loc        = LT_get_loc(ST_td_tbl, fld)
             mobility   = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)
             diff_coeff = LT_get_col_at_loc(ST_td_tbl, i_diffusion, loc)
             v_drift    = -mobility * fld
             gradc      = cc(i, j) - cc(i-dix(1), j-dix(2))

             if (v_drift < 0.0_dp) then
                gradn = cc(i+dix(1), j+dix(2)) - cc(i, j)
                boxes(id)%fc(i, j, dim, flux_elec) = v_drift * &
                     (cc(i, j) - koren_mlim(gradc, gradn))
             else                  ! v_drift > 0
                gradp = cc(i-dix(1), j-dix(2)) - cc(i-2*dix(1), j-2*dix(2))
                boxes(id)%fc(i, j, dim, flux_elec) = v_drift * &
                     (cc(i-dix(1), j-dix(2)) + koren_mlim(gradc, gradp))
             end if

             ! Diffusive part with 2-nd order explicit method. dif_f has to be
             ! scaled by 1/dx
             boxes(id)%fc(i, j, dim, flux_elec) = &
                  boxes(id)%fc(i, j, dim, flux_elec) - &
                  diff_coeff * gradc * inv_dr
          end do
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
    real(dp)                    :: alpha, eta, sflux, mu, rfac(2)
    integer                     :: i, j, nc, ioff
    type(LT_loc_t)              :: loc

    nc     = box%n_cell
    inv_dr = 1/box%dr
    ioff   = (box%ix(1)-1) * nc

    do j = 1, nc
       do i = 1, nc
          ! Weighting of flux contribution for cylindrical coordinates
          rfac  = [i+ioff-1, i+ioff] / (i+ioff-0.5_dp)
          fld   = box%cc(i,j, i_electric_fld)
          loc   = LT_get_loc(ST_td_tbl, fld)
          alpha = LT_get_col_at_loc(ST_td_tbl, i_alpha, loc)
          eta   = LT_get_col_at_loc(ST_td_tbl, i_eta, loc)
          mu    = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)

          ! Contribution of flux
          sflux = (box%fc(i, j, 2, flux_elec) - box%fc(i, j+1, 2, flux_elec) + &
               rfac(1) * box%fc(i, j, 1, flux_elec) - &
               rfac(2) * box%fc(i+1, j, 1, flux_elec)) * inv_dr * dt(1)

          ! Source term
          src = fld * mu * abs(box%cc(i, j, i_electron)) * (alpha - eta) *dt(1)

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
    real(dp), parameter       :: p_quench = 30.0e-3_dp
    real(dp)                  :: quench_fac

    ! Compute quench factor, because some excited species will be quenched by
    ! collisions, preventing the emission of a UV photon
    quench_fac = p_quench / (ST_gas_pressure + p_quench)

    ! Set photon production rate per cell, which is proportional to the
    ! ionization rate.
    call a2_loop_box_arg(tree, set_photoionization_rate, [eta * quench_fac], .true.)

    call photoi_set_src_2d(tree, ST_photoi_tbl, ST_rng, num_photons, &
         i_photo, i_photo, 0.25e-3_dp, .true., .true., 1e-9_dp, dt)

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
          call a2_prolong_quadratic(tree%boxes(p_id), tree%boxes(id), i_electron)
          call a2_prolong_quadratic(tree%boxes(p_id), tree%boxes(id), i_pos_ion)
          call a2_prolong_quadratic(tree%boxes(p_id), tree%boxes(id), i_phi)
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

end program streamer_cyl
