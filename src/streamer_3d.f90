!> Program to perform 3d streamer simulations
program streamer_3d

  use m_a3_types
  use m_a3_core
  use m_a3_ghostcell
  use m_a3_utils
  use m_a3_restrict
  use m_a3_multigrid
  use m_a3_output
  use m_write_silo
  use m_streamer

  implicit none

  integer                :: i, n
  character(len=ST_slen) :: fname
  logical                :: write_out

  type(a3_t)             :: tree ! This contains the full grid information
  type(mg3_t)            :: mg   ! Multigrid option struct
  type(ref_info_t)       :: ref_info

  call ST_create_config()
  call ST_read_config_files()
  call ST_load_config()

  ! Initialize the transport coefficients
  call ST_load_transport_data()

  ! Set the initial conditions from the configuration
  call ST_get_init_cond(3)

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
  mg%box_op      => mg3_auto_op
  mg%box_corr    => mg3_auto_corr
  mg%box_gsrb    => mg3_auto_gsrb

  ! This routine always needs to be called when using multigrid
  call mg3_init_mg(mg)

  ST_out_cnt = 0 ! Number of output files written
  ST_time    = 0 ! Simulation time (all times are in s)

  ! Set up the initial conditions
  do
     call a3_loop_box(tree, set_initial_condition)
     call compute_electric_field(tree, n_fmg_cycles, .false.)
     call a3_adjust_refinement(tree, refine_routine, ref_info, 4)
     if (ref_info%n_add == 0 .and. ref_info%n_rm == 0) exit !!
  end do

  call a3_print_info(tree)

  if (ST_photoi_enabled) &
       call set_photoionization(tree, ST_photoi_eta, ST_photoi_num_photons)

  do
     ! Get a new time step, which is at most dt_max
     call a3_reduction_vec(tree, get_max_dt, get_min, &
          [ST_dt_max, ST_dt_max, ST_dt_max], ST_dt_vec, ST_dt_num_cond)
     ST_dt = minval(ST_dt_vec)

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
        call a3_tree_copy_cc(tree, i_electron, i_electron_old)
        call a3_tree_copy_cc(tree, i_pos_ion, i_pos_ion_old)

        ! Two forward Euler steps over ST_dt
        do i = 1, 2
           ST_time = ST_time + ST_dt

           ! First calculate fluxes
           call a3_loop_boxes_arg(tree, fluxes_koren, [ST_dt], .true.)
           call a3_consistent_fluxes(tree, [flux_elec])

           ! Update the solution
           call a3_loop_box_arg(tree, update_solution, [ST_dt], .true.)

           ! Restrict the electron and ion densities to lower levels
           call a3_restrict_tree(tree, i_electron)

           ! Fill ghost cells
           call a3_gc_tree(tree, i_electron, a3_gc_interp_lim, a3_bc_neumann_zero)

           ! Compute new field on first iteration
           if (i == 1) call compute_electric_field(tree, n_fmg_cycles, .true.)
        end do

        ST_time = ST_time - ST_dt        ! Go back one time step

        ! Take average of phi_old and phi (explicit trapezoidal rule)
        call a3_loop_box(tree, average_density)

        ! Compute field with new density
        call compute_electric_field(tree, n_fmg_cycles, .true.)
     end do

     if (write_out) call a3_write_silo(tree, fname, ST_out_cnt, ST_time, &
          ixs_cc=[i_electron, i_pos_ion, i_electric_fld, i_photo, i_rhs], &
          dir=ST_output_dir)

     call a3_restrict_tree(tree, i_pos_ion)
     call a3_gc_tree(tree, i_pos_ion, a3_gc_interp_lim, a3_bc_neumann_zero)
     call a3_adjust_refinement(tree, refine_routine, ref_info, 4)

     if (ref_info%n_add > 0 .or. ref_info%n_rm > 0) then
        ! For boxes which just have been refined, set data on their children
        call prolong_to_new_boxes(tree, ref_info)

        ! Compute the field on the new mesh
        call compute_electric_field(tree, n_fmg_cycles, .true.)
     end if

  end do

  call a3_destroy(tree)

contains

  !> Initialize the AMR tree
  subroutine init_tree(tree)
    type(a3_t), intent(inout) :: tree

    ! Variables used below to initialize tree
    real(dp)                  :: dr
    integer                   :: id
    integer                   :: ix_list(3, 1) ! Spatial indices of initial boxes
    integer                   :: nb_list(6, 1) ! Neighbors of initial boxes
    integer                   :: n_boxes_init = 1000

    dr = ST_domain_len / ST_box_size

    ! Initialize tree
    call a3_init(tree, ST_box_size, n_var_cell, n_var_face, dr, &
         coarsen_to=2, n_boxes=n_boxes_init, cc_names=ST_cc_names)

    ! Set up geometry
    id             = 1          ! One box ...
    ix_list(:, id) = [1,1,1]    ! With index 1,1,1 ...
    nb_list(:, id) = -1         ! And neighbors -1 (physical boundary)

    ! Create the base mesh
    call a3_set_base(tree, ix_list, nb_list)

  end subroutine init_tree

  ! This routine sets the cell refinement flags for box
  subroutine refine_routine(box, cell_flags)
    use m_geometry
    type(box3_t), intent(in) :: box
    ! Refinement flags for the cells of the box
    integer, intent(out)     :: cell_flags(box%n_cell, box%n_cell, box%n_cell)
    integer                  :: i, j, k, n, nc
    real(dp)                 :: cphi, dx, dx2
    real(dp)                 :: alpha, adx, fld
    real(dp)                 :: dist

    nc      = box%n_cell
    dx      = box%dr
    dx2     = box%dr**2

    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             fld   = box%cc(i, j, k, i_electric_fld)
             alpha = LT_get_col(ST_td_tbl, i_alpha, fld)
             ! The refinement is based on the ionization length
             adx   = box%dr * alpha

             ! The refinement is also based on the intensity of the source term.
             ! Here we estimate the curvature of phi (given by dx**2 *
             ! Laplacian(phi))
             cphi = dx2 * abs(box%cc(i, j, k, i_rhs))

             if (adx / ST_refine_adx + cphi / ST_refine_cphi > 1) then
                cell_flags(i, j, k) = af_do_ref
             else if (adx < 0.125_dp * ST_refine_adx .and. &
                  cphi < 0.0625_dp * ST_refine_cphi &
                  .and. dx < ST_derefine_dx) then
                cell_flags(i, j, k) = af_rm_ref
             else
                cell_flags(i, j, k) = af_keep_ref
             end if

             ! Refine around the initial conditions
             if (ST_time < ST_refine_init_time) then
                do n = 1, ST_init_cond%n_cond
                   dist = GM_dist_line(a3_r_cc(box, [i, j, k]), &
                        ST_init_cond%seed_r0(:, n), &
                        ST_init_cond%seed_r1(:, n), 3)
                   if (dist - ST_init_cond%seed_width(n) < 2 * dx &
                        .and. box%dr > ST_refine_init_fac * &
                        ST_init_cond%seed_width(n)) then
                      cell_flags(i, j, k) = af_do_ref
                   end if
                end do
             end if

          end do
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
    type(box3_t), intent(inout) :: box
    integer                     :: i, j, k, n, nc
    real(dp)                    :: xy(3)
    real(dp)                    :: density

    nc                          = box%n_cell
    box%cc(:, :, :, i_electron) = ST_init_cond%background_density
    box%cc(:, :, :, i_pos_ion)  = ST_init_cond%background_density
    box%cc(:, :, :, i_phi)      = 0 ! Inital potential set to zero

    do k = 0, nc+1
       do j = 0, nc+1
          do i = 0, nc+1
             xy   = a3_r_cc(box, [i,j,k])

             do n = 1, ST_init_cond%n_cond
                density = ST_init_cond%seed_density(n) * &
                     GM_density_line(xy, ST_init_cond%seed_r0(:, n), &
                     ST_init_cond%seed_r1(:, n), 3, &
                     ST_init_cond%seed_width(n), &
                     ST_init_cond%seed_falloff(n))

                ! Add electrons and/or ions depending on the seed charge type
                ! (positive, negative or neutral)
                if (ST_init_cond%seed_charge_type(n) <= 0) then
                   box%cc(i, j, k, i_electron) = &
                        box%cc(i, j, k, i_electron) + density
                end if

                if (ST_init_cond%seed_charge_type(n) >= 0) then
                   box%cc(i, j, k, i_pos_ion) = &
                        box%cc(i, j, k, i_pos_ion) + density
                end if
             end do
          end do
       end do
    end do

  end subroutine set_initial_condition

  ! Get maximum time step based on e.g. CFL criteria
  function get_max_dt(box, n_cond) result(dt_vec)
    use m_units_constants
    type(box3_t), intent(in) :: box
    integer, intent(in)      :: n_cond
    integer                  :: i, j, k, nc
    real(dp)                 :: fld(3), fld_norm, mobility, diffusion_c
    real(dp)                 :: dt_vec(n_cond)

    nc = box%n_cell
    dt_vec = ST_dt_max

    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             fld(1) = 0.5_dp * (box%fx(i, j, k, electric_fld) + &
                  box%fx(i+1, j, k, electric_fld))
             fld(2) = 0.5_dp * (box%fy(i, j, k, electric_fld) + &
                  box%fy(i, j+1, k, electric_fld))
             fld(3) = 0.5_dp * (box%fz(i, j, k, electric_fld) + &
                  box%fz(i, j, k+1, electric_fld))

             fld_norm = box%cc(i, j, k, i_electric_fld)
             mobility = LT_get_col(ST_td_tbl, i_mobility, fld_norm)
             diffusion_c = LT_get_col(ST_td_tbl, i_diffusion, fld_norm)

             ! The 0.5 is here because of the explicit trapezoidal rule
             dt_vec(ST_ix_cfl) = min(dt_vec(ST_ix_cfl), &
                  0.5_dp/sum(abs(fld * mobility) / box%dr))

             ! Dielectric relaxation time
             dt_vec(ST_ix_drt) = min(dt_vec(ST_ix_drt), &
                  UC_eps0 / (UC_elem_charge * mobility * &
                  max(box%cc(i, j, k, i_electron), epsilon(1.0_dp))))

             ! Diffusion condition
             dt_vec(ST_ix_diff) = min(dt_vec(ST_ix_diff), &
                  0.25_dp * box%dr**2 / max(diffusion_c, epsilon(1.0_dp)))
          end do
       end do
    end do

  end function get_max_dt

  function get_min(a, b, n) result(min_vec)
    integer, intent(in)  :: n
    real(dp), intent(in) :: a(n), b(n)
    real(dp)             :: min_vec(n)

    min_vec = min(a, b)
  end function get_min

  ! !> Get maximum time step based on e.g. CFL criteria
  ! real(dp) function get_max_dt(tree)
  !   use m_units_constants
  !   type(a3_t), intent(in) :: tree
  !   real(dp)               :: max_fld, min_fld, max_dns, dr_min
  !   real(dp)               :: mobility, diff_coeff, alpha, max_mobility
  !   real(dp)               :: dt_cfl, dt_dif, dt_drt, dt_alpha

  !   call a3_tree_max_cc(tree, i_electric_fld, max_fld)
  !   call a3_tree_min_cc(tree, i_electric_fld, min_fld)
  !   call a3_tree_max_cc(tree, i_electron, max_dns)

  !   dr_min       = a3_min_dr(tree)
  !   mobility     = LT_get_col(ST_td_tbl, i_mobility, max_fld)
  !   max_mobility = LT_get_col(ST_td_tbl, i_mobility, min_fld)
  !   diff_coeff   = LT_get_col(ST_td_tbl, i_diffusion, max_fld)
  !   alpha        = LT_get_col(ST_td_tbl, i_alpha, max_fld)

  !   ! CFL condition
  !   dt_cfl = 0.5_dp * dr_min / (mobility * max_fld)

  !   ! Diffusion condition
  !   dt_dif = 0.25_dp * dr_min**2 / diff_coeff

  !   ! Dielectric relaxation time
  !   dt_drt = UC_eps0 / (UC_elem_charge * max_mobility * max_dns)

  !   ! Ionization limit
  !   dt_alpha =  1 / max(mobility * max_fld * alpha, epsilon(1.0_dp))

  !   get_max_dt = 0.95_dp * min(dt_cfl, dt_dif, dt_alpha, ST_dt_max)
  ! end function get_max_dt

  !> Compute electric field on the tree. First perform multigrid to get electric
  ! potential, then take numerical gradient to geld field.
  subroutine compute_electric_field(tree, n_cycles, have_guess)
    use m_units_constants
    type(a3_t), intent(inout) :: tree
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
          tree%boxes(id)%cc(:, :, :, i_rhs) = fac * (&
               tree%boxes(id)%cc(:, :, :, i_electron) - &
               tree%boxes(id)%cc(:, :, :, i_pos_ion))
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    call ST_set_voltage(ST_time)

    ! Perform n_cycles fmg cycles (logicals: store residual, first call)
    do i = 1, n_cycles
       call mg3_fas_fmg(tree, mg, .false., have_guess .or. i > 1)
    end do

    ! Compute field from potential
    call a3_loop_box(tree, electric_field_from_potential)

    ! Set the field norm also in ghost cells
    call a3_gc_tree(tree, i_electric_fld, a3_gc_interp, a3_bc_neumann_zero)
  end subroutine compute_electric_field

  !> Compute electric field from electrical potential
  subroutine electric_field_from_potential(box)
    type(box3_t), intent(inout) :: box
    integer                     :: nc
    real(dp)                    :: inv_dr

    nc     = box%n_cell
    inv_dr = 1 / box%dr

    box%fx(:, :, :, electric_fld) = inv_dr * &
         (box%cc(0:nc, 1:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, 1:nc, i_phi))
    box%fy(:, :, :, electric_fld) = inv_dr * &
         (box%cc(1:nc, 0:nc, 1:nc, i_phi) - box%cc(1:nc, 1:nc+1, 1:nc, i_phi))
    box%fz(:, :, :, electric_fld) = inv_dr * &
         (box%cc(1:nc, 1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc, 1:nc+1, i_phi))

    box%cc(1:nc, 1:nc, 1:nc, i_electric_fld) = 0.5_dp * sqrt(&
         (box%fx(1:nc, 1:nc, 1:nc, electric_fld) + box%fx(2:nc+1, 1:nc, 1:nc, electric_fld))**2 + &
         (box%fy(1:nc, 1:nc, 1:nc, electric_fld) + box%fy(1:nc, 2:nc+1, 1:nc, electric_fld))**2 + &
         (box%fz(1:nc, 1:nc, 1:nc, electric_fld) + box%fz(1:nc, 1:nc, 2:nc+1, electric_fld))**2)
  end subroutine electric_field_from_potential

  !> This fills ghost cells near physical boundaries for the potential
  subroutine sides_bc_potential(box, nb, iv, bc_type)
    type(box3_t), intent(inout) :: box
    integer, intent(in)         :: nb ! Direction for the boundary condition
    integer, intent(in)         :: iv ! Index of variable
    integer, intent(out)        :: bc_type ! Type of boundary condition
    integer                     :: nc

    nc = box%n_cell

    select case (nb)
    case (a3_neighb_lowx)
       bc_type = af_bc_neumann
       box%cc(   0, 1:nc, 1:nc, iv) = 0
    case (a3_neighb_highx)
       bc_type = af_bc_neumann
       box%cc(nc+1, 1:nc, 1:nc, iv) = 0
    case (a3_neighb_lowy)
       bc_type = af_bc_neumann
       box%cc(1:nc,    0, 1:nc, iv) = 0
    case (a3_neighb_highy)
       bc_type = af_bc_neumann
       box%cc(1:nc, nc+1, 1:nc, iv) = 0
    case (a3_neighb_lowz)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, 1:nc,    0, iv) = 0
    case (a3_neighb_highz)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, 1:nc, nc+1, iv) = ST_applied_voltage
    end select
  end subroutine sides_bc_potential

  !> Compute the electron fluxes due to drift and diffusion
  subroutine fluxes_koren(boxes, id, dt_vec)
    use m_units_constants
    type(box3_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp), intent(in)        :: dt_vec(:)
    real(dp)                    :: inv_dr, tmp, gradp, gradc, gradn
    real(dp)                    :: mobility, diff_coeff, v_drift
    real(dp)                    :: fld_avg, fld
    real(dp)                    :: gc_data(boxes(id)%n_cell, &
         boxes(id)%n_cell, a3_num_neighbors)
    integer                     :: i, j, k, nc
    type(LT_loc_t)              :: loc

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    call a3_gc2_box(boxes, id, i_electron, a3_gc2_prolong_linear, &
         a3_bc2_neumann_zero, gc_data, nc)

    ! x-fluxes interior, advective part with flux limiter
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc+1
             fld_avg    = 0.5_dp * (boxes(id)%cc(i, j, k, i_electric_fld) + &
                  boxes(id)%cc(i-1, j, k, i_electric_fld))
             loc        = LT_get_loc(ST_td_tbl, fld_avg)
             mobility   = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)
             diff_coeff = LT_get_col_at_loc(ST_td_tbl, i_diffusion, loc)
             fld        = boxes(id)%fx(i, j, k, electric_fld)
             v_drift    = -mobility * fld
             gradc      = boxes(id)%cc(i, j, k, i_electron) - &
                  boxes(id)%cc(i-1, j, k, i_electron)

             if (v_drift < 0.0_dp) then
                if (i == nc+1) then
                   tmp = gc_data(j, k, a3_neighb_highx)
                else
                   tmp = boxes(id)%cc(i+1, j, k, i_electron)
                end if
                gradn = tmp - boxes(id)%cc(i, j, k, i_electron)
                boxes(id)%fx(i, j, k, flux_elec) = v_drift * &
                     (boxes(id)%cc(i, j, k, i_electron) - koren_mlim(gradc, gradn))
             else                  ! v_drift > 0
                if (i == 1) then
                   tmp = gc_data(j, k, a3_neighb_lowx)
                else
                   tmp = boxes(id)%cc(i-2, j, k, i_electron)
                end if
                gradp = boxes(id)%cc(i-1, j, k, i_electron) - tmp
                boxes(id)%fx(i, j, k, flux_elec) = v_drift * &
                     (boxes(id)%cc(i-1, j, k, i_electron) + koren_mlim(gradc, gradp))
             end if

             ! Diffusive part with 2-nd order explicit method. dif_f has to be
             ! scaled by 1/dx
             boxes(id)%fx(i, j, k, flux_elec) = boxes(id)%fx(i, j, k, flux_elec) - &
                  diff_coeff * gradc * inv_dr
          end do
       end do
    end do

    ! y-fluxes interior, advective part with flux limiter
    do k = 1, nc
       do j = 1, nc+1
          do i = 1, nc
             fld_avg    = 0.5_dp * (boxes(id)%cc(i, j, k, i_electric_fld) + &
                  boxes(id)%cc(i, j-1, k, i_electric_fld))
             loc        = LT_get_loc(ST_td_tbl, fld_avg)
             mobility   = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)
             diff_coeff = LT_get_col_at_loc(ST_td_tbl, i_diffusion, loc)
             fld        = boxes(id)%fy(i, j, k, electric_fld)
             v_drift    = -mobility * fld
             gradc      = boxes(id)%cc(i, j, k, i_electron) - &
                  boxes(id)%cc(i, j-1, k, i_electron)

             if (v_drift < 0.0_dp) then
                if (j == nc+1) then
                   tmp = gc_data(i, k, a3_neighb_highy)
                else
                   tmp = boxes(id)%cc(i, j+1, k, i_electron)
                end if
                gradn = tmp - boxes(id)%cc(i, j, k, i_electron)
                boxes(id)%fy(i, j, k, flux_elec) = v_drift * &
                     (boxes(id)%cc(i, j, k, i_electron) - koren_mlim(gradc, gradn))
             else                  ! v_drift > 0
                if (j == 1) then
                   tmp = gc_data(i, k, a3_neighb_lowy)
                else
                   tmp = boxes(id)%cc(i, j-2, k, i_electron)
                end if
                gradp = boxes(id)%cc(i, j-1, k, i_electron) - tmp
                boxes(id)%fy(i, j, k, flux_elec) = v_drift * &
                     (boxes(id)%cc(i, j-1, k, i_electron) + koren_mlim(gradc, gradp))
             end if

             ! Diffusive part with 2-nd order explicit method. dif_f has to be
             ! scaled by 1/dx
             boxes(id)%fy(i, j, k, flux_elec) = boxes(id)%fy(i, j, k, flux_elec) - &
                  diff_coeff * gradc * inv_dr
          end do
       end do
    end do

    ! z-fluxes interior, advective part with flux limiter
    do k = 1, nc+1
       do j = 1, nc
          do i = 1, nc
             fld_avg    = 0.5_dp * (boxes(id)%cc(i, j, k, i_electric_fld) + &
                  boxes(id)%cc(i, j, k-1, i_electric_fld))
             loc        = LT_get_loc(ST_td_tbl, fld_avg)
             mobility   = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)
             diff_coeff = LT_get_col_at_loc(ST_td_tbl, i_diffusion, loc)
             fld        = boxes(id)%fz(i, j, k, electric_fld)
             v_drift    = -mobility * fld
             gradc      = boxes(id)%cc(i, j, k, i_electron) - &
                  boxes(id)%cc(i, j, k-1, i_electron)

             if (v_drift < 0.0_dp) then
                if (k == nc+1) then
                   tmp = gc_data(i, j, a3_neighb_highz)
                else
                   tmp = boxes(id)%cc(i, j, k+1, i_electron)
                end if
                gradn = tmp - boxes(id)%cc(i, j, k, i_electron)
                boxes(id)%fz(i, j, k, flux_elec) = v_drift * &
                     (boxes(id)%cc(i, j, k, i_electron) - koren_mlim(gradc, gradn))
             else                  ! v_drift > 0
                if (k == 1) then
                   tmp = gc_data(i, j, a3_neighb_lowz)
                else
                   tmp = boxes(id)%cc(i, j, k-2, i_electron)
                end if
                gradp = boxes(id)%cc(i, j, k-1, i_electron) - tmp
                boxes(id)%fz(i, j, k, flux_elec) = v_drift * &
                     (boxes(id)%cc(i, j, k-1, i_electron) + koren_mlim(gradc, gradp))
             end if

             ! Diffusive part with 2-nd order explicit method. dif_f has to be
             ! scaled by 1/dx
             boxes(id)%fz(i, j, k, flux_elec) = boxes(id)%fz(i, j, k, flux_elec) - &
                  diff_coeff * gradc * inv_dr
          end do
       end do
    end do

  end subroutine fluxes_koren

  !> Take average of new and old electron/ion density for explicit trapezoidal rule
  subroutine average_density(box)
    type(box3_t), intent(inout) :: box
    box%cc(:, :, :, i_electron) = 0.5_dp * &
         (box%cc(:, :, :, i_electron) + box%cc(:, :, :, i_electron_old))
    box%cc(:, :, :, i_pos_ion) = 0.5_dp * &
         (box%cc(:, :, :, i_pos_ion)  + box%cc(:, :, :, i_pos_ion_old))
  end subroutine average_density

  !> Advance solution over dt based on the fluxes / source term, using forward Euler
  subroutine update_solution(box, dt)
    type(box3_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr, src, fld
    real(dp)                    :: alpha, eta, sflux, mu
    integer                     :: i, j, k, nc, ioff
    type(LT_loc_t)              :: loc

    nc     = box%n_cell
    inv_dr = 1/box%dr
    ioff   = (box%ix(1)-1) * nc

    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             ! Weighting of flux contribution for cylindrical coordinates
             fld   = box%cc(i,j, k, i_electric_fld)
             loc   = LT_get_loc(ST_td_tbl, fld)
             alpha = LT_get_col_at_loc(ST_td_tbl, i_alpha, loc)
             eta   = LT_get_col_at_loc(ST_td_tbl, i_eta, loc)
             mu    = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)

             ! Contribution of flux
             sflux = &
                  box%fy(i, j, k, flux_elec) - box%fy(i, j+1, k, flux_elec) + &
                  box%fx(i, j, k, flux_elec) - box%fx(i+1, j, k, flux_elec) + &
                  box%fz(i, j, k, flux_elec) - box%fz(i, j, k+1, flux_elec)
             sflux = sflux * inv_dr * dt(1)

             ! Source term
             src = fld * mu * box%cc(i, j, k, i_electron) * (alpha - eta) * dt(1)

             if (ST_photoi_enabled) src = src + box%cc(i, j, k, i_photo) * dt(1)

             ! Add flux and source term
             box%cc(i, j, k, i_electron) = box%cc(i, j, k, i_electron) + sflux + src

             ! Add source term
             box%cc(i, j, k, i_pos_ion) = box%cc(i, j, k, i_pos_ion) + src
          end do
       end do
    end do

  end subroutine update_solution

  !> Sets the photoionization
  subroutine set_photoionization(tree, eta, num_photons, dt)
    use m_units_constants

    type(a3_t), intent(inout) :: tree
    real(dp), intent(in)      :: eta
    real(dp), intent(in), optional :: dt
    integer, intent(in)       :: num_photons
    real(dp), parameter       :: p_quench = 30e-3_dp ! 30 mbar
    real(dp)                  :: quench_fac

    ! Compute quench factor, because some excited species will be quenched by
    ! collisions, preventing the emission of a UV photon
    quench_fac = p_quench / (ST_gas_pressure + p_quench)

    ! Set photon production rate per cell, which is proportional to the
    ! ionization rate.
    call a3_loop_box_arg(tree, set_photoionization_rate, [eta * quench_fac], .true.)

    call photoi_set_src_3d(tree, ST_photoi_tbl, ST_rng, num_photons, &
         i_photo, i_photo, 0.25e-3_dp, .true., 1e-9_dp, dt)

  end subroutine set_photoionization

  !> Sets the photoionization_rate
  subroutine set_photoionization_rate(box, coeff)
    type(box3_t), intent(inout) :: box
    real(dp), intent(in)        :: coeff(:)
    integer                     :: i, j, k, nc
    real(dp)                    :: fld, alpha, mobility
    type(LT_loc_t)              :: loc

    nc = box%n_cell

    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             fld      = box%cc(i, j, k, i_electric_fld)
             loc      = LT_get_loc(ST_td_tbl, fld)
             alpha    = LT_get_col_at_loc(ST_td_tbl, i_alpha, loc)
             mobility = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)

             box%cc(i, j, k, i_photo) = fld * mobility * alpha * &
                  box%cc(i, j, k, i_electron) * coeff(1)
          end do
       end do
    end do
  end subroutine set_photoionization_rate

  !> For each box that gets refined, set data on its children using this routine
  subroutine prolong_to_new_boxes(tree, ref_info)
    use m_a3_prolong
    type(a3_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id, p_id

    !$omp parallel private(lvl, i, id, p_id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          p_id = tree%boxes(id)%parent
          call a3_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_electron)
          call a3_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_pos_ion)
          call a3_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_phi)
       end do
       !$omp end do

       !$omp do
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a3_gc_box(tree%boxes, id, i_electron, &
               a3_gc_interp_lim, a3_bc_neumann_zero)
          call a3_gc_box(tree%boxes, id, i_pos_ion, &
               a3_gc_interp_lim, a3_bc_neumann_zero)
          call a3_gc_box(tree%boxes, id, i_phi, &
               mg3_sides_rb, mg%sides_bc)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine prolong_to_new_boxes

end program streamer_3d
