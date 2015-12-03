module m_streamer

  use m_config
  use m_random
  use m_photons
  use m_lookup_table

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)
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
  integer, parameter :: i_pho      = 8 ! Phototionization rate
  integer, parameter :: i_eps      = 9 ! Level set function
  character(len=10)  :: cc_names(n_var_cell) = &
       [character(len=10) :: "elec", "pion", "elec_old", &
       "pion_old", "phi", "fld", "rhs", "pho", "eps"]

  ! Indices of face-centered variables
  integer, parameter :: n_var_face = 2
  integer, parameter :: f_elec     = 1 ! Electron flux
  integer, parameter :: f_fld      = 2 ! Electric field vector

  ! Indices of transport data
  integer, parameter :: n_var_td    = 4
  integer, parameter :: i_mobility  = 1
  integer, parameter :: i_diffusion = 2
  integer, parameter :: i_alpha     = 3
  integer, parameter :: i_eta       = 4

  ! How many multigrid FMG cycles we perform per time step
  integer, parameter :: n_fmg_cycles = 1

    type initcnd_t
     real(dp)              :: bg_dens
     integer               :: n_cond
     real(dp), allocatable :: seed_r0(:, :)
     real(dp), allocatable :: seed_r1(:, :)
     real(dp), allocatable :: seed_dens(:)
     real(dp), allocatable :: seed_width(:)
     integer, allocatable  :: seed_falloff(:)
  end type initcnd_t

  type(initcnd_t)   :: init_cond

  type(LT_table_t)  :: td_tbl             ! Table with transport data vs fld
  type(CFG_t)       :: sim_cfg            ! The configuration for the simulation
  type(RNG_t)       :: sim_rng            ! Random number generator

  logical           :: photoi_enabled     ! Whether we use phototionization
  real(dp)          :: photoi_frac_O2     ! Oxygen fraction
  real(dp)          :: photoi_eta         ! Photoionization efficiency
  integer           :: photoi_num_photons ! Number of photons to use
  type(PH_tbl_t)    :: photoi_tbl         ! Table for photoionization

  integer            :: i, n, n_steps_amr
  integer            :: output_cnt
  real(dp)           :: dt, time, end_time
  real(dp)           :: dt_output, dt_max
  character(len=40)  :: fname
  logical            :: write_out

  ! The size of the boxes that we use to construct our mesh
  integer :: box_size

  ! The length of the (square) domain
  real(dp) :: domain_len

  ! The applied electric field
  real(dp) :: applied_fld

  ! The applied voltage
  real(dp) :: applied_voltage

  ! Pressure of the gas in bar
  real(dp) :: gas_pressure

  ! Dielectric constant
  real(dp) :: epsilon_diel

contains

  subroutine create_cfg(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add(cfg, "end_time", 10.0d-9, &
         "The desired endtime in seconds of the simulation")
    call CFG_add(cfg, "sim_name", "sim", &
         "The name of the simulation")
    call CFG_add(cfg, "output_dir", "", &
         "Directory where the output should be written")
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
    call CFG_add(cfg, "epsilon_diel", 1.5_dp, &
         "The dielectric constant of the dielectric")

    call CFG_add(cfg, "fld_mod_t0", 1.0e99_dp, &
         "Modify electric field after this time")
    call CFG_add(cfg, "fld_sin_amplitude", 0.0_dp, &
         "Amplitude of sinusoidal modification")
    call CFG_add(cfg, "fld_sin_freq", 0.2e9_dp, &
         "Frequency of sinusoidal modification")
    call CFG_add(cfg, "fld_lin_deriv", 0.0_dp, &
         "Linear derivative of field")

    call CFG_add(cfg, "bg_dens", 1.0d12, &
         "The background ion and electron density in 1/m^3")
    call CFG_add(cfg, "seed_dens", [5.0d19], &
         "Initial density of the seed", .true.)
    call CFG_add(cfg, "seed_rel_r0", [0.5d0, 0.4d0], &
         "The relative start position of the initial seed", .true.)
    call CFG_add(cfg, "seed_rel_r1", [0.5d0, 0.6d0], &
         "The relative end position of the initial seed", .true.)
    call CFG_add(cfg, "seed_width", [0.5d-3], &
         "Seed width", .true.)
    call CFG_add(cfg, "seed_falloff", [1], &
         "Fallof type for seed, see m_geom.f90", .true.)

    call CFG_add(cfg, "dt_output", 1.0d-10, &
         "The timestep for writing output")
    call CFG_add(cfg, "dt_max", 1.0d-11, &
         "The maximum timestep")
    call CFG_add(cfg, "num_steps_amr", 2, &
         "The number of steps after which the mesh is updated")

    call CFG_add(cfg, "photoi_enabled", .true., &
         "Whether photoionization is enabled")
    call CFG_add(cfg, "photoi_frac_O2", 0.2_dp, &
         "Fraction of oxygen")
    call CFG_add(cfg, "photoi_eta", 0.05_dp, &
         "Photoionization efficiency factor")
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
  end subroutine create_cfg

  subroutine get_init_cond(cfg, cond, n_dim)
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
  end subroutine get_init_cond

  subroutine init_transport_coeff(cfg)
    use m_transport_data
    use m_config

    type(CFG_t), intent(in) :: cfg
    character(len=name_len) :: input_file, gas_name
    integer                 :: table_size
    real(dp)                :: max_fld, alpha_fac, eta_fac
    real(dp)                :: mobility_fac, diffusion_fac
    real(dp)                :: domain_len
    real(dp), allocatable   :: x_data(:), y_data(:)
    character(len=name_len) :: data_name

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
    call CFG_get(cfg, "gas_pressure", gas_pressure)
    call CFG_get(cfg, "photoi_enabled", photoi_enabled)
    call CFG_get(cfg, "photoi_frac_O2", photoi_frac_O2)
    call CFG_get(cfg, "photoi_eta", photoi_eta)
    call CFG_get(cfg, "photoi_num_photons", photoi_num_photons)

    if (photoi_enabled) then
       call CFG_get(cfg, "domain_len", domain_len)
       call PH_get_tbl_air(photoi_tbl, photoi_frac_O2 * gas_pressure, &
            2 * domain_len)
    end if

  end subroutine init_transport_coeff

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

end module m_streamer
