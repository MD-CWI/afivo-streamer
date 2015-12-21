module m_streamer

  use m_config
  use m_random
  use m_photons
  use m_lookup_table

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  ! Default length of strings
  integer, parameter :: ST_slen = 200

  ! ** Indices of cell-centered variables **
  integer, parameter :: n_var_cell = 9 ! Number of variables
  integer, parameter :: i_elec     = 1 ! Electron density
  integer, parameter :: i_pion     = 2 ! Positive ion density
  integer, parameter :: i_elec_old = 3 ! For time-stepping scheme
  integer, parameter :: i_pion_old = 4 ! For time-stepping scheme
  integer, parameter :: i_phi      = 5 ! Electrical potential
  integer, parameter :: i_fld      = 6 ! Electric field norm
  integer, parameter :: i_rhs      = 7 ! Source term Poisson
  integer, parameter :: i_pho      = 8 ! Phototionization rate
  integer, parameter :: i_eps      = 9 ! Level set function

  ! Names of the cell-centered variables
  character(len=10) :: cc_names(n_var_cell) = &
       [character(len=10) :: "elec", "pion", "elec_old", &
       "pion_old", "phi", "fld", "rhs", "pho", "eps"]

  ! ** Indices of face-centered variables **
  integer, parameter :: n_var_face = 2 ! Number of variables
  integer, parameter :: f_elec     = 1 ! Electron flux
  integer, parameter :: f_fld      = 2 ! Electric field vector

  ! ** Indices of transport data **
  integer, parameter :: n_var_td    = 4 ! Number of transport coefficients
  integer, parameter :: i_mobility  = 1 ! Electron mobility
  integer, parameter :: i_diffusion = 2 ! Electron diffusion constant
  integer, parameter :: i_alpha     = 3 ! Ionization coeff (1/m)
  integer, parameter :: i_eta       = 4 ! Attachment coeff (1/m)

  ! How many multigrid FMG cycles we perform per time step
  integer, parameter :: n_fmg_cycles = 1

  ! Type to store initial conditions in
  type initcnd_t
     real(dp)              :: bg_dens
     integer               :: n_cond
     real(dp), allocatable :: seed_r0(:, :)
     real(dp), allocatable :: seed_r1(:, :)
     real(dp), allocatable :: seed_dens(:)
     real(dp), allocatable :: seed_width(:)
     integer, allocatable  :: seed_falloff(:)
  end type initcnd_t

  ! This will contain the initial conditions
  type(initcnd_t)   :: ST_init_cond

  ! Table with transport data vs fld
  type(LT_table_t)  :: ST_td_tbl

  ! The configuration for the simulation
  type(CFG_t)       :: ST_cfg

  ! Random number generator
  type(RNG_t)       :: ST_rng

  ! Whether we use phototionization
  logical           :: ST_photoi_enabled

  ! Oxygen fraction
  real(dp)          :: ST_photoi_frac_O2

  ! Photoionization efficiency
  real(dp)          :: ST_photoi_eta

  ! Number of photons to use
  integer           :: ST_photoi_num_photons

  ! Table for photoionization
  type(PH_tbl_t)    :: ST_photoi_tbl

  ! Start modifying the background field after this times
  real(dp)          :: ST_fld_mod_t0

  ! Amplitude of sinusoidal modification
  real(dp)          :: ST_fld_sin_amplitude

  ! Frequency (Hz) of sinusoidal modification
  real(dp)          :: ST_fld_sin_freq

  ! Linear derivative of background field
  real(dp)          :: ST_fld_lin_deriv

  ! Name of the simulations
  character(len=ST_slen) :: ST_sim_name

  ! Output directory
  character(len=ST_slen) :: ST_output_dir

  ! Number of steps between updating the mesh
  integer  :: ST_steps_amr

  ! Number of output files written
  integer  :: ST_out_cnt

  ! Current time step
  real(dp) :: ST_dt

  ! Maximum allowed time step
  real(dp) :: ST_dt_max

  ! Time between writing output
  real(dp) :: ST_dt_out

  ! Current time
  real(dp) :: ST_time

  ! End time of the simulation
  real(dp) :: ST_end_time

  ! The size of the boxes that we use to construct our mesh
  integer  :: ST_box_size

  ! The length of the (square) domain
  real(dp) :: ST_domain_len

  ! The applied electric field
  real(dp) :: ST_applied_fld

  ! The applied voltage
  real(dp) :: ST_applied_voltage

  ! Pressure of the gas in bar
  real(dp) :: ST_gas_pressure

  ! Dielectric constant
  real(dp) :: ST_epsilon_diel

contains

  subroutine ST_create_cfg(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add(cfg, "end_time", 10.0d-9, &
         "The desired endtime (s) of the simulation")
    call CFG_add(cfg, "sim_name", "sim", &
         "The name of the simulation")
    call CFG_add(cfg, "output_dir", "", &
         "Directory where the output should be written")
    call CFG_add(cfg, "box_size", 8, &
         "The number of grid cells per coordinate in a box")
    call CFG_add(cfg, "domain_len", 32e-3_dp, &
         "The length of the domain (m)")
    call CFG_add(cfg, "gas_name", "N2", &
         "The name of the gas mixture used")
    call CFG_add(cfg, "gas_pressure", 1.0_dp, &
         "The gas pressure (bar), used for photoionization")
    call CFG_add(cfg, "applied_fld", 1.0d7, &
         "The applied electric field (V/m)")
    call CFG_add(cfg, "epsilon_diel", 1.5_dp, &
         "The relative dielectric constant (if a dielectric is specified)")

    call CFG_add(cfg, "fld_mod_t0", 1.0e99_dp, &
         "Modify electric field after this time (s)")
    call CFG_add(cfg, "fld_sin_amplitude", 0.0_dp, &
         "Amplitude of sinusoidal modification (V/m)")
    call CFG_add(cfg, "fld_sin_freq", 0.2e9_dp, &
         "Frequency of sinusoidal modification (Hz)")
    call CFG_add(cfg, "fld_lin_deriv", 0.0_dp, &
         "Linear derivative of field [V/(ms)]")

    call CFG_add(cfg, "bg_dens", 1.0d12, &
         "The background ion and electron density (1/m3)")
    call CFG_add(cfg, "seed_dens", [5.0d19], &
         "Initial density of the seed (1/m3)", .true.)
    call CFG_add(cfg, "seed_rel_r0", [0.5d0, 0.4d0], &
         "The relative start position of the initial seed", .true.)
    call CFG_add(cfg, "seed_rel_r1", [0.5d0, 0.6d0], &
         "The relative end position of the initial seed", .true.)
    call CFG_add(cfg, "seed_width", [0.5d-3], &
         "Seed width (m)", .true.)
    call CFG_add(cfg, "seed_falloff", [1], &
         "Fallof type for seed, see m_geom.f90", .true.)

    call CFG_add(cfg, "dt_output", 1.0d-10, &
         "The timestep for writing output (s)")
    call CFG_add(cfg, "dt_max", 1.0d-11, &
         "The maximum timestep (s)")
    call CFG_add(cfg, "num_steps_amr", 2, &
         "The number of steps after which the mesh is updated")

    call CFG_add(cfg, "photoi_enabled", .true., &
         "Whether photoionization is enabled")
    call CFG_add(cfg, "photoi_frac_O2", 0.2_dp, &
         "Fraction of oxygen (0-1)")
    call CFG_add(cfg, "photoi_eta", 0.05_dp, &
         "Photoionization efficiency factor, typically around 0.05")
    call CFG_add(cfg, "photoi_num_photons", 50*1000, &
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

    call CFG_add(cfg, "td_alpha_fac", 1.0_dp, &
         "Modify alpha by this factor")
    call CFG_add(cfg, "td_eta_fac", 1.0_dp, &
         "Modify eta by this factor")
    call CFG_add(cfg, "td_mobility_fac", 1.0_dp, &
         "Modify mobility by this factor")
    call CFG_add(cfg, "td_diffusion_fac", 1.0_dp, &
         "Modify diffusion by this factor")
  end subroutine ST_create_cfg

  subroutine ST_get_init_cond(cfg, cond, n_dim)
    type(CFG_t), intent(in)        :: cfg
    type(initcnd_t), intent(inout) :: cond
    integer, intent(in)            :: n_dim
    integer                        :: n_cond, varsize
    real(dp)                       :: dlen
    real(dp), allocatable          :: tmp_vec(:)

    call CFG_get(cfg, "bg_dens", cond%bg_dens)
    call CFG_get(cfg, "domain_len", dlen)

    call CFG_get_size(cfg, "seed_dens", n_cond)

    call CFG_get_size(cfg, "seed_rel_r0", varsize)
    if (varsize /= n_dim * n_cond) &
         stop "seed_... variables have incompatible size"

    call CFG_get_size(cfg, "seed_rel_r1", varsize)
    if (varsize /= n_dim * n_cond) &
         stop "seed_... variables have incompatible size"

    call CFG_get_size(cfg, "seed_width", varsize)
    if (varsize /= n_cond) &
         stop "seed_... variables have incompatible size"

    cond%n_cond = n_cond
    allocate(cond%seed_dens(n_cond))
    allocate(cond%seed_r0(n_dim, n_cond))
    allocate(cond%seed_r1(n_dim, n_cond))
    allocate(cond%seed_width(n_cond))
    allocate(cond%seed_falloff(n_cond))

    allocate(tmp_vec(n_dim * n_cond))
    call CFG_get(cfg, "seed_rel_r0", tmp_vec)
    cond%seed_r0 = dlen * reshape(tmp_vec, [n_dim, n_cond])
    call CFG_get(cfg, "seed_rel_r1", tmp_vec)
    cond%seed_r1 = dlen * reshape(tmp_vec, [n_dim, n_cond])

    call CFG_get(cfg, "seed_dens", cond%seed_dens)
    call CFG_get(cfg, "seed_width", cond%seed_width)
    call CFG_get(cfg, "seed_falloff", cond%seed_falloff)
  end subroutine ST_get_init_cond

  subroutine ST_load_transport_data(cfg, td_tbl, photoi_tbl)
    use m_transport_data
    use m_config

    type(CFG_t), intent(in)         :: cfg
    type(LT_table_t), intent(inout) :: td_tbl
    type(PH_tbl_t), intent(inout)   :: photoi_tbl

    character(len=ST_slen) :: input_file, gas_name
    integer                 :: table_size
    real(dp)                :: max_fld, alpha_fac, eta_fac
    real(dp)                :: mobility_fac, diffusion_fac
    real(dp), allocatable   :: x_data(:), y_data(:)
    character(len=ST_slen) :: data_name

    call CFG_get(cfg, "input_file", input_file)
    call CFG_get(cfg, "gas_name", gas_name)

    call CFG_get(cfg, "lkptbl_size", table_size)
    call CFG_get(cfg, "lkptbl_max_fld", max_fld)

    call CFG_get(cfg, "td_alpha_fac", alpha_fac)
    call CFG_get(cfg, "td_eta_fac", eta_fac)
    call CFG_get(cfg, "td_mobility_fac", mobility_fac)
    call CFG_get(cfg, "td_diffusion_fac", diffusion_fac)

    ! Create a lookup table for the model coefficients
    td_tbl = LT_create(0.0_dp, max_fld, table_size, n_var_td)

    ! Fill table with data
    call CFG_get(cfg, "td_mobility_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * mobility_fac
    call LT_set_col(td_tbl, i_mobility, x_data, y_data)

    call CFG_get(cfg, "td_diffusion_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * diffusion_fac
    call LT_set_col(td_tbl, i_diffusion, x_data, y_data)

    call CFG_get(cfg, "td_alpha_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * alpha_fac
    call LT_set_col(td_tbl, i_alpha, x_data, y_data)

    call CFG_get(cfg, "td_eta_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * eta_fac
    call LT_set_col(td_tbl, i_eta, x_data, y_data)

    ! Create table for photoionization
    if (ST_photoi_enabled) then
       call PH_get_tbl_air(photoi_tbl, ST_photoi_frac_O2 * ST_gas_pressure, &
            2 * ST_domain_len)
    end if

  end subroutine ST_load_transport_data

  !> Modified implementation of Koren limiter, to avoid division and the min/max
  !> functions, which can be problematic / expensive. In most literature, you
  !> have r = a / b (ratio of gradients). Then the limiter phi(r) is multiplied
  !> with b. With this implementation, you get phi(r) * b
  elemental function koren_mlim(a, b) result(bphi)
    real(dp), intent(in) :: a  !< Density gradient (numerator)
    real(dp), intent(in) :: b  !< Density gradient (denominator)
    real(dp), parameter  :: sixth = 1/6.0_dp
    real(dp)             :: bphi, aa, ab

    aa = a * a
    ab = a * b

    if (ab <= 0) then
       ! a and b have different sign or one of them is zero, so r is either 0,
       ! inf or negative (special case a == b == 0 is ignored)
       bphi = 0
    else if (aa <= 0.25_dp * ab) then
       ! 0 < a/b <= 1/4, limiter has value a/b
       bphi = a
    else if (aa <= 2.5_dp * ab) then
       ! 1/4 < a/b <= 2.5, limiter has value (1+2*a/b)/6
       bphi = sixth * (b + 2*a)
    else
       ! (1+2*a/b)/6 >= 1, limiter has value 1
       bphi = b
    end if
  end function koren_mlim

  subroutine ST_read_cfg_files(cfg)
    type(CFG_t), intent(inout) :: cfg
    character(len=ST_slen)     :: tmp_name, prev_name, cfg_name
    integer                    :: n

    ST_sim_name = ""
    prev_name = ""
    do n = 1, command_argument_count()
       call get_command_argument(n, cfg_name)
       call CFG_read_file(cfg, trim(cfg_name))

       call CFG_get(cfg, "sim_name", tmp_name)
       if (ST_sim_name == "") then
          ST_sim_name = tmp_name
       else if (tmp_name /= "" .and. tmp_name /= prev_name) then
          ST_sim_name = trim(ST_sim_name) // "_" // trim(tmp_name)
       end if
       prev_name = tmp_name
    end do

  end subroutine ST_read_cfg_files

  subroutine ST_load_cfg(cfg)
    type(CFG_t), intent(inout) :: cfg
    character(len=ST_slen)    :: tmp_name

    call CFG_get(cfg, "end_time", ST_end_time)
    call CFG_get(cfg, "box_size", ST_box_size)
    call CFG_get(cfg, "output_dir", ST_output_dir)
    call CFG_get(cfg, "domain_len", ST_domain_len)
    call CFG_get(cfg, "applied_fld", ST_applied_fld)
    call CFG_get(cfg, "dt_output", ST_dt_out)
    call CFG_get(cfg, "num_steps_amr", ST_steps_amr)
    call CFG_get(cfg, "dt_max", ST_dt_max)
    call CFG_get(cfg, "epsilon_diel", ST_epsilon_diel)

    call CFG_get(cfg, "gas_pressure", ST_gas_pressure)
    call CFG_get(cfg, "photoi_enabled", ST_photoi_enabled)
    call CFG_get(cfg, "photoi_frac_O2", ST_photoi_frac_O2)
    call CFG_get(cfg, "photoi_eta", ST_photoi_eta)
    call CFG_get(cfg, "photoi_num_photons", ST_photoi_num_photons)

    call CFG_get(cfg, "fld_mod_t0", ST_fld_mod_t0)
    call CFG_get(cfg, "fld_sin_amplitude", ST_fld_sin_amplitude)
    call CFG_get(cfg, "fld_sin_freq", ST_fld_sin_freq)
    call CFG_get(cfg, "fld_lin_deriv", ST_fld_lin_deriv)

    ST_applied_voltage = -ST_domain_len * ST_applied_fld

    tmp_name = trim(ST_output_dir) // "/" // trim(ST_sim_name) // "_config.txt"
    print *, "Settings written to ", trim(tmp_name)
    call CFG_write(cfg, trim(tmp_name))

  end subroutine ST_load_cfg

  function ST_get_fld(time) result(fld)
    use m_units_constants

    real(dp), intent(in) :: time
    real(dp)             :: fld

    if (time > ST_fld_mod_t0) then
       fld = ST_applied_fld + (time - ST_fld_mod_t0) * &
            ST_fld_lin_deriv + ST_fld_sin_amplitude * &
            sin((time - ST_fld_mod_t0) * 2 * UC_pi * ST_fld_sin_freq)
    else
       fld = ST_applied_fld
    end if
  end function ST_get_fld

end module m_streamer
