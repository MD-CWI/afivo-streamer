program streamer_2d

  use m_afivo_2d
  use m_mg_2d
  use m_write_silo
  use m_lookup_table
  use m_config
  use m_random

  implicit none

  integer, parameter :: dp       = kind(0.0d0)
  integer, parameter :: name_len = 200

  ! Indices of cell-centered variables
  integer, parameter :: n_var_cell = 9
  integer, parameter :: i_elec     = 1 ! Electron density
  integer, parameter :: i_pion     = 2 ! Positive ion density
  integer, parameter :: i_elec_old = 3 ! For time-stepping scheme
  integer, parameter :: i_pion_old = 4 ! For time-stepping scheme
  integer, parameter :: i_phi      = 5 ! Electrical potential
  integer, parameter :: i_fld      = 6 ! Electric field norm
  integer, parameter :: i_rhs      = 7 ! Source term Poisson
  integer, parameter :: i_res      = 8 ! Residual (multigrid)
  integer, parameter :: i_pho      = 9 ! Phototionization rate
  character(len=10)  :: cc_names(n_var_cell) = &
       [character(len=10) :: "elec", "pion", "elec_old", &
       "pion_old", "phi", "fld", "rhs", "res", "pho"]

  ! Indices of face-centered variables
  integer, parameter :: n_var_face = 2
  integer, parameter :: f_elec     = 1 ! Electron flux
  integer, parameter :: f_fld      = 2 ! Electric field vector

  ! Indices of transport data
  integer, parameter :: n_var_td = 4
  integer, parameter :: i_mobility = 1
  integer, parameter :: i_diffusion = 2
  integer, parameter :: i_alpha = 3
  integer, parameter :: i_eta = 4

  character(len=name_len) :: sim_name
  character(len=name_len) :: cfg_name, tmp_name, prev_name

  type initcnd_t
     real(dp) :: bg_dens

     real(dp) :: seed_r0(2)
     real(dp) :: seed_r1(2)
     real(dp) :: seed_dens
     real(dp) :: swidth

     real(dp) :: line_r0(2)
     real(dp) :: line_r1(2)
     real(dp) :: line_dens
     real(dp) :: lwidth
  end type initcnd_t

  type(initcnd_t) :: init_cond

  type(LT_table_t)   :: td_tbl  ! Table with transport data vs fld
  type(CFG_t)        :: sim_cfg ! The configuration for the simulation
  type(a2_t)         :: tree    ! This contains the full grid information
  type(mg2_t)        :: mg      ! Multigrid option struct
  type(RNG_t)        :: sim_rng ! Random number generator
  type(ref_info_t)   :: ref_info

  logical          :: photoi_enabled     ! Whether we use phototionization
  real(dp)         :: photoi_frac_O2     ! Oxygen fraction
  real(dp)         :: photoi_eta         ! Photoionization efficiency
  integer          :: photoi_num_photons ! Number of photons to use
  type(LT_table_t) :: photoi_tbl         ! Table for photoionization

  integer            :: i, n, n_steps
  integer            :: output_cnt
  real(dp)           :: dt, time, end_time
  real(dp)           :: dt_amr, dt_output
  character(len=40)  :: fname
  logical            :: write_out

  ! How many multigrid FMG cycles we perform per time step
  integer, parameter :: n_fmg_cycles = 1

  ! The size of the boxes that we use to construct our mesh
  integer :: box_size

  ! The length of the (square) domain
  real(dp) :: domain_len

  ! The applied electric field
  real(dp) :: applied_fld

  ! Pressure of the gas in bar
  real(dp) :: gas_pressure

  call create_cfg(sim_cfg)

  sim_name = ""
  prev_name = ""
  do n = 1, command_argument_count()
     call get_command_argument(n, cfg_name)
     call CFG_read_file(sim_cfg, trim(cfg_name))

     call CFG_get(sim_cfg, "sim_name", tmp_name)
     if (sim_name == "") then
        sim_name = tmp_name
     else if (tmp_name /= "" .and. tmp_name /= prev_name) then
        sim_name = trim(sim_name) // "_" // trim(tmp_name)
     end if
     prev_name = tmp_name
  end do

  call initialize(sim_cfg)

  call CFG_get(sim_cfg, "end_time", end_time)
  call CFG_get(sim_cfg, "sim_name", sim_name)
  call CFG_get(sim_cfg, "box_size", box_size)
  call CFG_get(sim_cfg, "domain_len", domain_len)
  call CFG_get(sim_cfg, "applied_fld", applied_fld)
  call CFG_get(sim_cfg, "dt_output", dt_output)
  call CFG_get(sim_cfg, "dt_amr", dt_amr)

  call get_init_cond(sim_cfg, init_cond)

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the multigrid options. First define the variables to use
  mg%i_phi        = i_phi
  mg%i_tmp        = i_fld
  mg%i_rhs        = i_rhs
  mg%i_res        = i_res

  ! The number of cycles at the lowest level
  mg%n_cycle_base = 8

  ! Routines to use for ...
  mg%sides_bc     => sides_bc_pot ! Filling ghost cell on physical boundaries

  ! This routine always needs to be called when using multigrid
  call mg2_init_mg(mg)

  output_cnt = 0          ! Number of output files written
  time       = 0          ! Simulation time (all times are in s)

  ! Set up the initial conditions
  do i = 1, 10
     call a2_loop_box(tree, set_init_cond)
     call a2_restrict_tree(tree, i_rhs)
     call compute_fld(tree, n_fmg_cycles)
     call a2_adjust_refinement(tree, set_ref_flags, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  if (photoi_enabled) &
       call set_photoionization(tree, photoi_eta, photoi_num_photons)

  do
     ! Get a new time step, which is at most dt_amr
     dt      = get_max_dt(tree)
     n_steps = ceiling(dt_amr/dt)
     dt      = dt_amr / n_steps

     if (dt < 1e-14) then
        print *, "dt getting too small, instability?"
        time = end_time + 1.0_dp
     end if

     ! Every dt_output, write output
     if (output_cnt * dt_output <= time) then
        write_out = .true.
        output_cnt = output_cnt + 1
        write(fname, "(A,I0,A)") trim(sim_name) // "_", &
             output_cnt, ".silo"
     else
        write_out = .false.
     end if

     if (write_out) call a2_write_silo(tree, trim(fname), &
          cc_names, output_cnt, time)

     if (time > end_time) exit

     ! We perform n_steps between mesh-refinements
     do n = 1, n_steps
        time = time + dt

        ! Copy previous solution
        call a2_tree_copy_cc(tree, i_elec, i_elec_old)
        call a2_tree_copy_cc(tree, i_pion, i_pion_old)

        ! Two forward Euler steps over dt
        do i = 1, 2
           ! First calculate fluxes
           call a2_loop_boxes(tree, fluxes_koren, .true.)

           call compute_fld(tree, n_fmg_cycles)

           ! Update the solution
           call a2_loop_box_arg(tree, update_solution, [dt], .true.)

           ! Restrict the electron and ion densities to lower levels
           call a2_restrict_tree(tree, i_elec)
           call a2_restrict_tree(tree, i_pion)

           ! Fill ghost cells
           call a2_gc_sides(tree, i_elec, a2_sides_interp, sides_bc_dens)
           call a2_gc_sides(tree, i_pion, a2_sides_interp, sides_bc_dens)
        end do

        ! Take average of phi_old and phi (explicit trapezoidal rule)
        call a2_loop_box(tree, average_dens, .true.)
     end do

     call a2_adjust_refinement(tree, set_ref_flags, ref_info)

     if (ref_info%n_add > 0 .or. ref_info%n_rm > 0) then
        ! For boxes which just have been refined, set data on their children
        call prolong_to_new_children(tree, ref_info)

        ! Compute the field on the
        call compute_fld(tree, n_fmg_cycles)

        ! This will every now-and-then clean up the data in the tree
        call a2_tidy_up(tree, 0.9_dp, 0.5_dp, 5000, .false.)
     end if

     if (photoi_enabled) &
          call set_photoionization(tree, photoi_eta, photoi_num_photons)
  end do

  call a2_destroy(tree)

contains

  ! Initialize the AMR tree
  subroutine init_tree(tree)
    type(a2_t), intent(inout) :: tree

    ! Variables used below to initialize tree
    real(dp)                  :: dr
    integer                   :: id
    integer                   :: ix_list(2, 1) ! Spatial indices of initial boxes
    integer                   :: nb_list(4, 1) ! Neighbors of initial boxes
    integer                   :: n_boxes_init = 10*1000

    dr = domain_len / box_size

    ! Initialize tree
    call a2_init(tree, box_size, n_var_cell, n_var_face, dr, &
         coarsen_to=4, n_boxes = n_boxes_init)

    ! Set up geometry
    id             = 1          ! One box ...
    ix_list(:, id) = [1,1]      ! With index 1,1 ...
    nb_list(:, id) = -1         ! And neighbors -1 (physical boundary)

    ! Create the base mesh
    call a2_set_base(tree, ix_list, nb_list)

  end subroutine init_tree

  ! Refinement function
  subroutine set_ref_flags(boxes, id, ref_flags)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flags(:)
    integer                  :: nc
    real(dp)                 :: crv_phi, dr2, max_edens, max_fld

    nc        = boxes(id)%n_cell
    dr2       = boxes(id)%dr**2
    crv_phi   = dr2 * maxval(abs(boxes(id)%cc(1:nc, 1:nc, i_rhs)))
    max_edens = maxval(boxes(id)%cc(1:nc, 1:nc, i_elec))
    max_fld   = maxval(boxes(id)%cc(1:nc, 1:nc, i_fld))

    if (crv_phi < 4.0_dp) &
         ref_flags(id) = a5_rm_ref

    if (max_edens > 100 * init_cond%bg_dens .and. &
         boxes(id)%dr > 0.2e-3_dp) &
         ref_flags(id) = a5_do_ref

    if (max_edens > 100 * init_cond%bg_dens .and. &
         boxes(id)%dr > 0.1e-3_dp) &
         ref_flags(id) = max(ref_flags(id), a5_kp_ref)

    if (boxes(id)%dr > 0.75e-3_dp .or. (boxes(id)%dr > 1e-4_dp .and. &
         a2_r_inside(boxes(id), init_cond%seed_r1, 1.0e-3_dp))) &
         ref_flags(id) = a5_do_ref

    if (crv_phi > 2.0e1_dp .and. max_fld > 2e6_dp &
         .and. boxes(id)%dr > 5e-6_dp) &
         ref_flags(id) = a5_do_ref

  end subroutine set_ref_flags

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)
    real(dp)                    :: dens

    nc = box%n_cell

    do j = 0, nc+1
       do i = 0, nc+1
          xy   = a2_r_cc(box, [i,j])
          dens = init_cond%seed_dens * rod_dens(xy, init_cond%seed_r0, &
               init_cond%seed_r1, init_cond%swidth, 3)
          dens = dens + init_cond%line_dens * rod_dens(xy, init_cond%line_r0, &
               init_cond%line_r1, init_cond%lwidth, 4)

          box%cc(i, j, i_elec) = init_cond%bg_dens + dens
          box%cc(i, j, i_pion) = init_cond%bg_dens + dens
       end do
    end do

    box%cc(:, :, i_phi) = 0     ! Inital potential set to zero
  end subroutine set_init_cond

  ! Get maximum time step based on e.g. CFL criteria
  real(dp) function get_max_dt(tree)
    type(a2_t), intent(in) :: tree
    real(dp), parameter    :: UC_eps0        = 8.8541878176d-12
    real(dp), parameter    :: UC_elem_charge = 1.6022d-19
    real(dp)               :: max_fld, min_fld, max_dns, dr_min
    real(dp)               :: mobility, diff_coeff, alpha, max_mobility
    real(dp)               :: dt_cfl, dt_dif, dt_drt, dt_alpha

    call a2_tree_max_cc(tree, i_fld, max_fld)
    call a2_tree_min_cc(tree, i_fld, min_fld)
    call a2_tree_max_cc(tree, i_elec, max_dns)

    dr_min       = a2_min_dr(tree)
    mobility     = LT_get_col(td_tbl, i_mobility, max_fld)
    max_mobility = LT_get_col(td_tbl, i_mobility, min_fld)
    diff_coeff   = LT_get_col(td_tbl, i_diffusion, max_fld)
    alpha        = LT_get_col(td_tbl, i_alpha, max_fld)

    ! CFL condition
    dt_cfl = dr_min / (mobility * max_fld)

    ! Diffusion condition
    dt_dif = dr_min**2 / diff_coeff

    ! Dielectric relaxation time
    dt_drt = UC_eps0 / (UC_elem_charge * max_mobility * max_dns)

    ! Ionization limit
    dt_alpha =  1 / (mobility * max_fld * alpha)

    get_max_dt = 0.8_dp * min(1/(1/dt_cfl + 1/dt_dif), dt_drt, dt_alpha)
  end function get_max_dt

  ! Used to create initial electron/ion seed
  real(dp) function rod_dens(xy, xy0, xy1, sigma, falloff_type)
    real(dp), intent(in) :: xy(2), xy0(2), xy1(2), sigma
    integer, intent(in)  :: falloff_type
    real(dp)             :: distance, temp

    distance = dist_line(xy, xy0, xy1)

    select case (falloff_type)
    case (1)                    ! Sigmoid
       rod_dens    = 2 / (1 + exp(distance / sigma))
    case (2)                    ! Gaussian
       rod_dens    = exp(-(distance/sigma)**2)
    case (3)                    ! Smooth-step
       if (distance < sigma) then
          rod_dens = 1
       else if (distance < 2 * sigma) then
          temp = distance/sigma - 1
          rod_dens = (1- (3 * temp**2 - 2 * temp**3))
       else
          rod_dens = 0.0_dp
       end if
    case (4)                    ! Hard boundary
       if (distance < sigma) then
          rod_dens = 1
       else
          rod_dens = 0
       end if
    case default
       rod_dens = 0.0_dp
    end select
  end function rod_dens

  real(dp) function dist_line(xy, xy0, xy1)
    real(dp), intent(in) :: xy(2), xy0(2), xy1(2)
    real(dp) :: line_len2, temp
    real(dp) :: projection(2)

    line_len2 = sum((xy1 - xy0)**2)
    temp = sum((xy - xy0) * (xy1 - xy0)) / line_len2

    if (temp < 0.0_dp) then
       dist_line = sqrt(sum((xy-xy0)**2))
    else if (temp > 1.0_dp) then
       dist_line = sqrt(sum((xy-xy1)**2))
    else
       projection = xy0 + temp * (xy1 - xy0)
       dist_line = sqrt(sum((xy-projection)**2))
    end if
  end function dist_line

  ! Compute electric field on the tree. First perform multigrid to get electric
  ! potential, then take numerical gradient to geld field.
  subroutine compute_fld(tree, n_fmg)
    use m_units_constants
    type(a2_t), intent(inout) :: tree
    integer, intent(in) :: n_fmg
    real(dp), parameter :: fac = UC_elem_charge / UC_eps0
    integer :: lvl, i, id, nc

    nc = tree%n_cell

    ! Set the source term (rhs)
    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          tree%boxes(id)%cc(:, :, i_rhs) = fac * (&
               tree%boxes(id)%cc(:, :, i_elec) - &
               tree%boxes(id)%cc(:, :, i_pion))
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    ! Restrict the rhs
    call a2_restrict_tree(tree, i_rhs)

    ! Perform n_fmg full-multigrid cycles
    do i = 1, n_fmg
       call mg2_fas_fmg(tree, mg)
    end do

    ! Compute field from potential
    call a2_loop_box(tree, fld_from_pot)

    ! Set the field norm also in ghost cells
    call a2_gc_sides(tree, i_fld, a2_sides_interp, sides_bc_dens)
  end subroutine compute_fld

  ! Compute electric field from electrical potential
  subroutine fld_from_pot(box)
    type(box2_t), intent(inout) :: box
    integer                     :: nc
    real(dp)                    :: inv_dr

    nc     = box%n_cell
    inv_dr = 1 / box%dr

    box%fx(:, :, f_fld) = inv_dr * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi))
    box%fy(:, :, f_fld) = inv_dr * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi))

    box%cc(1:nc, 1:nc, i_fld) = sqrt(&
         0.25_dp * (box%fx(1:nc, 1:nc, f_fld) + box%fx(2:nc+1, 1:nc, f_fld))**2 + &
         0.25_dp * (box%fy(1:nc, 1:nc, f_fld) + box%fy(1:nc, 2:nc+1, f_fld))**2)
  end subroutine fld_from_pot

  ! Koren limiter for use in advection scheme
  elemental function limiter_koren(theta)
    real(dp), intent(in) :: theta
    real(dp)             :: limiter_koren
    real(dp), parameter  :: one_sixth = 1.0_dp / 6.0_dp
    limiter_koren = max(0.0d0, min(1.0_dp, theta, &
         (1.0_dp + 2.0_dp * theta) * one_sixth))
  end function limiter_koren

  ! Safe division for use in Koren limiter
  elemental function ratio(numerator, denominator)
    real(dp), intent(in) :: numerator, denominator
    real(dp)             :: ratio
    if (denominator /= 0.0_dp) then
       ratio = numerator / denominator
    else
       ratio = 1 / epsilon(1.0_dp)
    end if
  end function ratio

  ! Compute the electron fluxes due to drift and diffusion
  subroutine fluxes_koren(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp)                    :: inv_dr, theta
    real(dp)                    :: gradp, gradc, gradn
    real(dp) :: fld_avg, mobility, diff_coeff, v_drift
    integer                     :: i, j, nc, nb_id
    type(LT_loc_t) :: loc

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    ! x-fluxes interior, advective part with flux limiter
    do j = 1, nc
       do i = 1, nc+1
          fld_avg   = 0.5_dp * (boxes(id)%cc(i, j, i_fld) + &
               boxes(id)%cc(i-1, j, i_fld))
          loc        = LT_get_loc(td_tbl, fld_avg)
          mobility   = LT_get_col_at_loc(td_tbl, i_mobility, loc)
          diff_coeff = LT_get_col_at_loc(td_tbl, i_diffusion, loc)
          v_drift    = -mobility * boxes(id)%fx(i, j, f_fld)

          gradc = boxes(id)%cc(i, j, i_elec) - boxes(id)%cc(i-1, j, i_elec)
          if (v_drift < 0.0_dp) then

             if (i == nc+1) then
                nb_id = boxes(id)%neighbors(a2_nb_hx)
                if (nb_id > a5_no_box) then
                   gradn = boxes(nb_id)%cc(2, j, i_elec) - boxes(id)%cc(i, j, i_elec)
                else
                   gradn = 0
                end if
             else
                gradn = boxes(id)%cc(i+1, j, i_elec) - boxes(id)%cc(i, j, i_elec)
             end if

             theta = ratio(gradc, gradn)
             boxes(id)%fx(i, j, f_elec) = v_drift * &
                  (boxes(id)%cc(i, j, i_elec) - limiter_koren(theta) * gradn)
          else                  ! v_drift > 0

             if (i == 1) then
                nb_id = boxes(id)%neighbors(a2_nb_lx)
                if (nb_id > a5_no_box) then
                   gradp = boxes(id)%cc(i-1, j, i_elec) - boxes(nb_id)%cc(nc-1, j, i_elec)
                else
                   gradp = 0
                end if
             else
                gradp = boxes(id)%cc(i-1, j, i_elec) - boxes(id)%cc(i-2, j, i_elec)
             end if

             theta = ratio(gradc, gradp)
             boxes(id)%fx(i, j, f_elec) = v_drift * &
                  (boxes(id)%cc(i-1, j, i_elec) + limiter_koren(theta) * gradp)
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
          boxes(id)%fx(i, j, f_elec) = boxes(id)%fx(i, j, f_elec) - &
               diff_coeff * gradc * inv_dr
       end do
    end do

    ! y-fluxes interior, advective part with flux limiter
    do j = 1, nc+1
       do i = 1, nc
          fld_avg   = 0.5_dp * (boxes(id)%cc(i, j, i_fld) + &
               boxes(id)%cc(i, j-1, i_fld))
          loc        = LT_get_loc(td_tbl, fld_avg)
          mobility   = LT_get_col_at_loc(td_tbl, i_mobility, loc)
          diff_coeff = LT_get_col_at_loc(td_tbl, i_diffusion, loc)
          v_drift    = -mobility * boxes(id)%fy(i, j, f_fld)

          gradc = boxes(id)%cc(i, j, i_elec) - boxes(id)%cc(i, j-1, i_elec)

          if (v_drift < 0.0_dp) then
             if (j == nc+1) then
                nb_id = boxes(id)%neighbors(a2_nb_hy)
                if (nb_id > a5_no_box) then
                   gradn = boxes(nb_id)%cc(i, 2, i_elec) - boxes(id)%cc(i, j, i_elec)
                else
                   gradn = 0
                end if
             else
                gradn = boxes(id)%cc(i, j+1, i_elec) - boxes(id)%cc(i, j, i_elec)
             end if

             theta = ratio(gradc, gradn)
             boxes(id)%fy(i, j, f_elec) = v_drift * &
                  (boxes(id)%cc(i, j, i_elec) - limiter_koren(theta) * gradn)
          else                  ! v_drift > 0

             if (j == 1) then
                nb_id = boxes(id)%neighbors(a2_nb_ly)
                if (nb_id > a5_no_box) then
                   gradp = boxes(id)%cc(i, j-1, i_elec) - boxes(nb_id)%cc(i, nc-1, i_elec)
                else
                   gradp = 0
                end if
             else
                gradp = boxes(id)%cc(i, j-1, i_elec) - boxes(id)%cc(i, j-2, i_elec)
             end if

             theta = ratio(gradc, gradp)
             boxes(id)%fy(i, j, f_elec) = v_drift * &
                  (boxes(id)%cc(i, j-1, i_elec) + limiter_koren(theta) * gradp)
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
          boxes(id)%fy(i, j, f_elec) = boxes(id)%fy(i, j, f_elec) - &
               diff_coeff * gradc * inv_dr
       end do
    end do

  end subroutine fluxes_koren

  ! Take average of new and old electron/ion density for explicit trapezoidal rule
  subroutine average_dens(box)
    type(box2_t), intent(inout) :: box
    box%cc(:, :, i_elec) = 0.5_dp * (box%cc(:, :, i_elec) + box%cc(:, :, i_elec_old))
    box%cc(:, :, i_pion) = 0.5_dp * (box%cc(:, :, i_pion) + box%cc(:, :, i_pion_old))
  end subroutine average_dens

  ! Advance solution over dt based on the fluxes / source term, using forward Euler
  subroutine update_solution(box, dt)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr, src, fld
    real(dp)                    :: alpha, eta, mobility
    integer                     :: i, j, nc
    type(LT_loc_t) :: loc

    nc                    = box%n_cell
    inv_dr                = 1/box%dr
    do j = 1, nc
       do i = 1, nc
          fld      = box%cc(i,j, i_fld)
          loc      = LT_get_loc(td_tbl, fld)
          alpha    = LT_get_col_at_loc(td_tbl, i_alpha, loc)
          eta      = LT_get_col_at_loc(td_tbl, i_eta, loc)
          mobility = LT_get_col_at_loc(td_tbl, i_mobility, loc)

          src = abs(mobility * fld) * (alpha-eta) * &
               dt(1) * box%cc(i, j, i_elec)
          if (photoi_enabled) &
               src = src + dt(1) * box%cc(i,j, i_pho)

          box%cc(i, j, i_elec) = box%cc(i, j, i_elec) + src
          box%cc(i, j, i_pion) = box%cc(i, j, i_pion) + src
       end do
    end do

    box%cc(1:nc, 1:nc, i_elec) = box%cc(1:nc, 1:nc, i_elec) + dt(1) * ( &
         (box%fx(1:nc, :, f_elec) - box%fx(2:nc+1, :, f_elec)) * inv_dr + &
         (box%fy(:, 1:nc, f_elec) - box%fy(:, 2:nc+1, f_elec)) * inv_dr)

  end subroutine update_solution

  subroutine set_photoionization(tree, eta, num_photons)
    use m_photons
    use m_units_constants

    type(a2_t), intent(inout) :: tree
    real(dp), intent(in)      :: eta
    integer, intent(in)       :: num_photons
    real(dp), parameter       :: p_quench = 30.0D0 * UC_torr_to_bar
    real(dp)                  :: quench_fac

    ! Compute quench factor, because some excited species will be quenched by
    ! collisions, preventing the emission of a UV photon
    quench_fac = p_quench / (gas_pressure + p_quench)

    ! Set photon production rate per cell, which is proportional to the
    ! ionization rate.
    call a2_loop_box_arg(tree, set_photoi_rate, [eta * quench_fac], .true.)

    call PH_set_src(tree, photoi_tbl, sim_rng, num_photons, i_pho, i_pho)

  end subroutine set_photoionization

  subroutine set_photoi_rate(box, coeff)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: coeff(:)
    integer                     :: i, j, nc
    real(dp)                    :: fld, alpha, mobility, dr
    type(LT_loc_t)              :: loc

    nc = box%n_cell

    do j = 1, nc
       do i = 1, nc
          dr       = box%dr
          fld      = box%cc(i, j, i_fld)
          loc      = LT_get_loc(td_tbl, fld)
          alpha    = LT_get_col_at_loc(td_tbl, i_alpha, loc)
          mobility = LT_get_col_at_loc(td_tbl, i_mobility, loc)

          box%cc(i, j, i_pho) = fld * mobility * alpha * &
               box%cc(i, j, i_elec) * coeff(1) * dr**2
       end do
    end do
  end subroutine set_photoi_rate

  ! For each box that gets refined, set data on its children using this routine
  subroutine prolong_to_new_children(tree, ref_info)
    type(a2_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id

    do lvl = 1, tree%max_lvl
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a2_prolong1_to(tree%boxes, id, i_elec)
          call a2_prolong1_to(tree%boxes, id, i_pion)
          call a2_prolong1_to(tree%boxes, id, i_phi)
       end do

       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a2_gc_box_sides(tree%boxes, id, i_elec, &
               a2_sides_interp, sides_bc_dens)
          call a2_gc_box_sides(tree%boxes, id, i_pion, &
               a2_sides_interp, sides_bc_dens)
          call a2_gc_box_sides(tree%boxes, id, i_phi, &
               a2_sides_extrap, sides_bc_pot)
       end do
    end do
  end subroutine prolong_to_new_children

  ! This fills ghost cells near physical boundaries for the potential
  subroutine sides_bc_pot(boxes, id, nb, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    integer                     :: nc
    real(dp)                    :: v0

    nc = boxes(id)%n_cell
    v0 = applied_fld * domain_len

    select case (nb)
    case (a2_nb_lx)
       boxes(id)%cc(0, 1:nc, iv) = boxes(id)%cc(1, 1:nc, iv)
    case (a2_nb_hx)
       boxes(id)%cc(nc+1, 1:nc, iv) = boxes(id)%cc(nc, 1:nc, iv)
    case (a2_nb_ly)
       boxes(id)%cc(1:nc, 0, iv) = -2 * v0 - boxes(id)%cc(1:nc, 1, iv)
    case (a2_nb_hy)
       boxes(id)%cc(1:nc, nc+1, iv) = -boxes(id)%cc(1:nc, nc, iv)
    end select
  end subroutine sides_bc_pot

  ! This fills ghost cells near physical boundaries for the electron density
  subroutine sides_bc_dens(boxes, id, nb, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    integer                     :: nc

    nc = boxes(id)%n_cell

    select case (nb)
    case (a2_nb_lx)
       ! Neumann zero
       boxes(id)%cc(0, 1:nc, iv) = boxes(id)%cc(1, 1:nc, iv)
    case (a2_nb_hx)
       ! Neumann zero
       boxes(id)%cc(nc+1, 1:nc, iv) = boxes(id)%cc(nc, 1:nc, iv)
    case (a2_nb_ly)
       ! Neumann zero
       boxes(id)%cc(1:nc, 0, iv) = boxes(id)%cc(1:nc, 1, iv)
    case (a2_nb_hy)
       ! Neumann zero
       boxes(id)%cc(1:nc, nc+1, iv) = boxes(id)%cc(1:nc, nc, iv)
    end select
  end subroutine sides_bc_dens

  subroutine get_init_cond(cfg, cond)
    type(CFG_t), intent(in)        :: cfg
    type(initcnd_t), intent(inout) :: cond
    real(dp)                       :: dlen

    call CFG_get(cfg, "init_bg_dens", cond%bg_dens)
    call CFG_get(cfg, "domain_len", dlen)

    call CFG_get(cfg, "init_seed_dens", cond%seed_dens)
    call CFG_get(cfg, "init_seed_rel_r0", cond%seed_r0)
    call CFG_get(cfg, "init_seed_rel_r1", cond%seed_r1)
    call CFG_get(cfg, "init_seed_width", cond%swidth)
    cond%seed_r0 = cond%seed_r0 * dlen
    cond%seed_r1 = cond%seed_r1 * dlen

    call CFG_get(cfg, "init_line_dens", cond%line_dens)
    call CFG_get(cfg, "init_line_rel_r0", cond%line_r0)
    call CFG_get(cfg, "init_line_rel_r1", cond%line_r1)
    call CFG_get(cfg, "init_line_width", cond%lwidth)
    cond%line_r0 = cond%line_r0 * dlen
    cond%line_r1 = cond%line_r1 * dlen
  end subroutine get_init_cond

  subroutine create_cfg(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add(cfg, "end_time", 10.0d-9, &
         "The desired endtime in seconds of the simulation")
    call CFG_add(cfg, "sim_name", "sim", &
         "The name of the simulation")
    call CFG_add(cfg, "box_size", 8, &
         "The number of grid cells per coordinate in a box")
    call CFG_add(cfg, "domain_len", 32e-3_dp, &
         "The length of the (square) domain")
    call CFG_add(cfg, "gas_name", "N2", &
         "The name of the gas mixture used")
    call CFG_add(cfg, "gas_pressure", 1.0_dp, &
         "The gas pressure in bar (used for photoionization)")
    call CFG_add(cfg, "applied_fld", 1.0d7, &
         "The applied electric field")

    call CFG_add(cfg, "init_bg_dens", 1.0d12, &
         "The background ion and electron density in 1/m^3")
    call CFG_add(cfg, "init_seed_dens", 5.0d19 , &
         "Initial density of the seed")
    call CFG_add(cfg, "init_seed_rel_r0", [0.5d0, 0.0d0], &
         "The relative start position of the initial seed")
    call CFG_add(cfg, "init_seed_rel_r1", [0.5d0, 0.4d0], &
         "The relative end position of the initial seed")
    call CFG_add(cfg, "init_seed_width", 0.5d-3, &
         "Seed width")
    call CFG_add(cfg, "init_line_dens", 5.0d15 , &
         "Initial density of the line")
    call CFG_add(cfg, "init_line_rel_r0", [0.0d0, 0.7d0], &
         "The relative start position of the initial line")
    call CFG_add(cfg, "init_line_rel_r1", [1.0d0, 0.7d0], &
         "The relative end position of the initial line")
    call CFG_add(cfg, "init_line_width", 0.5d-3, &
         "Line width")

    call CFG_add(cfg, "dt_output", 1.0d-10, &
         "The timestep for writing output")
    call CFG_add(cfg, "dt_amr", 1.0d-11, &
         "The timestep for adaptively refining the grid")

    call CFG_add(cfg, "photoi_enabled", .true., &
         "Whether photoionization is enabled")
    call CFG_add(cfg, "photoi_frac_O2", 0.2_dp, &
         "Fraction of oxygen")
    call CFG_add(cfg, "photoi_eta", 0.05_dp, &
         "Photoionization efficiency factor")
    call CFG_add(cfg, "photoi_num_photons", 10*1000, &
         "Number of discrete photons to use for photoionization")

    call CFG_add(cfg, "input_file", "transport_data_file.txt", &
         "Input file with transport data")
    call CFG_add(cfg, "lkptbl_size", 1000, &
         "The transport data table size in the fluid model")
    call CFG_add(cfg, "lkptbl_max_fld", 3.0d7, &
         "The maximum electric field in the fluid model coefficients")
    call CFG_add(cfg, "td_mobility_name", "efield[V/m]_vs_mu[m2/Vs]", &
         "The name of the mobility coefficient")
    call CFG_add(cfg, "td_diffusion_name", "efield[V/m]_vs_dif[m2/s]", &
         "The name of the diffusion coefficient")
    call CFG_add(cfg, "td_alpha_name", "efield[V/m]_vs_alpha[1/m]", &
         "The name of the eff. ionization coeff.")
    call CFG_add(cfg, "td_eta_name", "efield[V/m]_vs_eta[1/m]", &
         "The name of the eff. attachment coeff.")
  end subroutine create_cfg

  subroutine initialize(cfg)
    use m_transport_data
    use m_photons
    use m_config

    type(CFG_t), intent(in) :: cfg
    character(len=name_len) :: input_file, gas_name
    integer                 :: table_size
    real(dp)                :: max_fld
    real(dp), allocatable   :: x_data(:), y_data(:)
    character(len=name_len) :: data_name

    call CFG_get(cfg, "input_file", input_file)
    call CFG_get(cfg, "gas_name", gas_name)

    call CFG_get(cfg, "lkptbl_size", table_size)
    call CFG_get(cfg, "lkptbl_max_fld", max_fld)

    ! Create a lookup table for the model coefficients
    td_tbl = LT_create(0.0_dp, max_fld, table_size, n_var_td)

    ! Fill table with data
    call CFG_get(cfg, "td_mobility_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    call LT_set_col(td_tbl, i_mobility, x_data, y_data)

    call CFG_get(cfg, "td_diffusion_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    call LT_set_col(td_tbl, i_diffusion, x_data, y_data)

    call CFG_get(cfg, "td_alpha_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    call LT_set_col(td_tbl, i_alpha, x_data, y_data)

    call CFG_get(cfg, "td_eta_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    call LT_set_col(td_tbl, i_eta, x_data, y_data)

    ! Create table for photoionization
    call CFG_get(cfg, "gas_pressure", gas_pressure)
    call CFG_get(cfg, "photoi_enabled", photoi_enabled)
    call CFG_get(cfg, "photoi_frac_O2", photoi_frac_O2)
    call CFG_get(cfg, "photoi_eta", photoi_eta)
    call CFG_get(cfg, "photoi_num_photons", photoi_num_photons)

    if (photoi_enabled) &
       photoi_tbl = PH_get_tbl_air(photoi_frac_O2 * gas_pressure)

  end subroutine initialize

end program streamer_2d
