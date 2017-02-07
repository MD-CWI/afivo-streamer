!> Module  that contains subroutines to handle with configuration files
!!
!! * Creating configuration files
!! * Reading configuration files
!! * Loading configuration files
!! * Setting the initial conditions from the configuration
!!
!! Further it contains a subroutine for
!! * loading the transport data
!! * Setting  the voltage on fixed time
!!
!! Moreover, it contains several pre-defined variables like 
!!
!! * Indices of cell-centered variables
!! * Names of the cell-centered variables
!! * Indices of face-centered variables
!! * Indices of transport data
!! * How many multigrid FMG cycles we perform per time step
!!
!! The subroutines and variables are used by the examples:
!! * streamer_2d
!! * streamer_3d
!! * streamer_cyl
!!
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
  integer, parameter :: n_var_cell     = 8 ! Number of variables
  integer, parameter :: i_electron     = 1 ! Electron density
  integer, parameter :: i_pos_ion      = 2 ! Positive ion density
  integer, parameter :: i_electron_old = 3 ! For time-stepping scheme
  integer, parameter :: i_pos_ion_old  = 4 ! For time-stepping scheme
  integer, parameter :: i_phi          = 5 ! Electrical potential
  integer, parameter :: i_electric_fld = 6 ! Electric field norm
  integer, parameter :: i_rhs          = 7 ! Source term Poisson
  integer, parameter :: i_photo        = 8 ! Phototionization rate

  ! Names of the cell-centered variables
  character(len=12) :: ST_cc_names(n_var_cell) = &
       [character(len=12) :: "electron", "pos_ion", "electron_old", &
       "pos_ion_old", "phi", "electric_fld", "rhs", "pho"]

  ! ** Indices of face-centered variables **
  integer, parameter :: n_var_face   = 2 ! Number of variables
  integer, parameter :: flux_elec    = 1 ! Electron flux
  integer, parameter :: electric_fld = 2 ! Electric field vector

  ! ** Indices of transport data **
  integer, parameter :: n_var_td    = 4 ! Number of transport coefficients
  integer, parameter :: i_mobility  = 1 ! Electron mobility
  integer, parameter :: i_diffusion = 2 ! Electron diffusion constant
  integer, parameter :: i_alpha     = 3 ! Ionization coeff (1/m)
  integer, parameter :: i_eta       = 4 ! Attachment coeff (1/m)

  ! Whether cylindrical coordinates are used
  logical :: ST_cylindrical = .false.

  ! Table with transport data vs electric field
  type(lookup_table_t), protected :: ST_td_tbl

  ! Random number generator
  type(RNG_t) :: ST_rng

  ! Whether we use phototionization
  logical, protected :: ST_photoi_enabled

  ! Oxygen fraction
  real(dp), protected :: ST_photoi_frac_O2

  ! Photoionization efficiency
  real(dp), protected :: ST_photoi_eta

  ! Number of photons to use
  integer, protected :: ST_photoi_num_photons

  ! Table for photoionization
  type(photoi_tbl_t), protected :: ST_photoi_tbl

  ! Name of the simulations
  character(len=ST_slen), protected :: ST_simulation_name

  ! Output directory
  character(len=ST_slen), protected :: ST_output_dir

  ! The number of steps after which the mesh is updated
  integer, protected :: ST_refine_per_steps

  ! The grid spacing will always be larger than this value
  real(dp), protected :: ST_refine_min_dx

  ! The grid spacing will always be smaller than this value
  real(dp), protected :: ST_refine_max_dx

  ! Refine if alpha*dx is larger than this value
  real(dp), protected :: ST_refine_adx

  ! Refine if the curvature in phi is larger than this value
  real(dp), protected :: ST_refine_cphi

  ! Only derefine if grid spacing if smaller than this value
  real(dp), protected :: ST_derefine_dx

  ! Refine around initial conditions up to this time
  real(dp), protected :: ST_refine_init_time

  ! Refine until dx is smaller than this factor times the seed width
  real(dp), protected :: ST_refine_init_fac

  ! Number of output files written
  integer :: ST_out_cnt

  ! Current time step
  real(dp) :: ST_dt

  ! Number of time step restrictions
  integer, parameter :: ST_dt_num_cond = 3

  ! Array of time step restrictions
  real(dp) :: ST_dt_vec(ST_dt_num_cond)

  ! Index of CFL condition
  integer, parameter :: ST_ix_cfl = 1

  ! Index of diffusion time step condition
  integer, parameter :: ST_ix_diff = 2

  ! Index of dielectric relaxation time condition
  integer, parameter :: ST_ix_drt = 3

  ! Maximum allowed time step
  real(dp), protected :: ST_dt_max

  ! Time between writing output
  real(dp), protected :: ST_dt_out

  ! Current time
  real(dp)  :: ST_time

  ! End time of the simulation
  real(dp), protected :: ST_end_time

  ! The size of the boxes that we use to construct our mesh
  integer, protected :: ST_box_size

  ! The length of the (square) domain
  real(dp), protected :: ST_domain_len

  ! Pressure of the gas in bar
  real(dp), protected :: ST_gas_pressure

  ! Dielectric constant
  real(dp), protected :: ST_epsilon_diel

contains

  !> Create the configuration file with default values
  subroutine ST_initialize(cfg)
    type(CFG_t), intent(inout) :: cfg !< The configuration for the simulation

    call CFG_add_get(cfg, "cylindrical", ST_cylindrical, &
         "Whether cylindrical coordinates are used (only in 2D)")
    call CFG_add(cfg, "end_time", ST_end_time, &
         "The desired endtime (s) of the simulation")
    call CFG_add(cfg, "simulation_name", ST_simulation_name, &
         "The name of the simulation")
    call CFG_add(cfg, "output_dir", ST_output_dir, &
         "Directory where the output should be written")
    call CFG_add(cfg, "box_size", ST_box_size, &
         "The number of grid cells per coordinate in a box")
    call CFG_add(cfg, "domain_len", ST_domain_len-3_dp, &
         "The length of the domain (m)")
    call CFG_add(cfg, "gas_name", ST_gas_name, &
         "The name of the gas mixture used")
    call CFG_add(cfg, "gas_pressure", ST_gas_pressure, &
         "The gas pressure (bar), used for photoionization")

    call CFG_add(cfg, "dt_output", 1.0d-10, &
         "The timestep for writing output (s)")
    call CFG_add(cfg, "dt_max", 1.0d-11, &
         "The maximum timestep (s)")

    call CFG_add(cfg, "refine_per_steps", 2, &
         "The number of steps after which the mesh is updated")
    call CFG_add(cfg, "refine_min_dx", 1.0e-6_dp, &
         "The grid spacing will always be larger than this value")
    call CFG_add(cfg, "refine_max_dx", 1.0e-3_dp, &
         "The grid spacing will always be smaller than this value")

    call CFG_add(cfg, "refine_adx", 1.0_dp, &
         "Refine if alpha*dx is larger than this value")
    call CFG_add(cfg, "refine_cphi", 1e99_dp, &
         "Refine if the curvature in phi is larger than this value")
    call CFG_add(cfg, "derefine_dx", 1e-4_dp, &
         "Only derefine if grid spacing if smaller than this value")
    call CFG_add(cfg, "refine_init_time", 10.0e-9_dp, &
         "Refine around initial conditions up to this time")
    call CFG_add(cfg, "refine_init_fac", 0.25_dp, &
         "Refine until dx is smaller than this factor times the seed width")

    call CFG_add(cfg, "photoi_enabled", .true., &
         "Whether photoionization is enabled")
    call CFG_add(cfg, "photoi_frac_O2", 0.2_dp, &
         "Fraction of oxygen (0-1)")
    call CFG_add(cfg, "photoi_eta", 0.05_dp, &
         "Photoionization efficiency factor, typically around 0.05")
    call CFG_add(cfg, "photoi_num_photons", 50*1000, &
         "Number of discrete photons to use for photoionization")
    call CFG_add(cfg, "photoi_rng_seed", [8123, 91234, 12399, 293434], &
         "Seed for the photoionization random number generator")

    call CFG_add(cfg, "input_file", "transport_data_file.txt", &
         "Input file with transport data")
    call CFG_add(cfg, "lookup_table_size", 1000, &
         "The transport data table size in the fluid model")
    call CFG_add(cfg, "lookup_table_max_electric_fld", 3.0d7, &
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
  end subroutine ST_create_config

  ! Set the initial conditions from the configuration
  subroutine ST_get_init_cond(n_dim)
    integer, intent(in)            :: n_dim
    integer                        :: n_cond, varsize
    real(dp)                       :: dlen
    real(dp), allocatable          :: tmp_vec(:)
    type(initcnd_t) :: ic

    call CFG_get(cfg, "background_density", ic%background_density)
    call CFG_get(cfg, "domain_len", dlen)

    call CFG_get_size(cfg, "seed_density", n_cond)

    call CFG_get_size(cfg, "seed_rel_r0", varsize)
    if (varsize /= n_dim * n_cond) &
         stop "seed_rel_r0 variable has incompatible size"

    call CFG_get_size(cfg, "seed_rel_r1", varsize)
    if (varsize /= n_dim * n_cond) &
         stop "seed_rel_r1 variable has incompatible size"

    call CFG_get_size(cfg, "seed_charge_type", varsize)
    if (varsize /= n_cond) &
         stop "seed_charge_type variable has incompatible size"

    call CFG_get_size(cfg, "seed_width", varsize)
    if (varsize /= n_cond) &
         stop "seed_width variable has incompatible size"

    ic%n_cond = n_cond
    allocate(ic%seed_density(n_cond))
    allocate(ic%seed_charge_type(n_cond))
    allocate(ic%seed_r0(n_dim, n_cond))
    allocate(ic%seed_r1(n_dim, n_cond))
    allocate(ic%seed_width(n_cond))
    allocate(ic%seed_falloff(n_cond))

    allocate(tmp_vec(n_dim * n_cond))
    call CFG_get(cfg, "seed_rel_r0", tmp_vec)
    ic%seed_r0 = dlen * reshape(tmp_vec, [n_dim, n_cond])
    call CFG_get(cfg, "seed_rel_r1", tmp_vec)
    ic%seed_r1 = dlen * reshape(tmp_vec, [n_dim, n_cond])

    call CFG_get(cfg, "seed_density", ic%seed_density)
    call CFG_get(cfg, "seed_charge_type", ic%seed_charge_type)
    call CFG_get(cfg, "seed_width", ic%seed_width)
    call CFG_get(cfg, "seed_falloff", ic%seed_falloff)
    ST_init_cond = ic
  end subroutine ST_get_init_cond

  !> Initialize the transport coefficients
  subroutine ST_load_transport_data()
    use iso_fortran_env, only: int64
    use m_transport_data
    use m_config

    character(len=ST_slen) :: input_file, gas_name
    integer                :: table_size, rng_int4_seed(4)
    integer(int64)         :: rng_int8_seed(2)
    real(dp)               :: max_electric_fld, alpha_fac, eta_fac
    real(dp)               :: mobility_fac, diffusion_fac
    real(dp), allocatable  :: x_data(:), y_data(:)
    character(len=ST_slen) :: data_name

    call CFG_get(cfg, "input_file", input_file)
    call CFG_get(cfg, "gas_name", gas_name)

    call CFG_get(cfg, "lookup_table_size", table_size)
    call CFG_get(cfg, "lookup_table_max_electric_fld", max_electric_fld)

    call CFG_get(cfg, "td_alpha_fac", alpha_fac)
    call CFG_get(cfg, "td_eta_fac", eta_fac)
    call CFG_get(cfg, "td_mobility_fac", mobility_fac)
    call CFG_get(cfg, "td_diffusion_fac", diffusion_fac)

    ! Create a lookup table for the model coefficients
    ST_td_tbl = LT_create(0.0_dp, max_electric_fld, table_size, n_var_td)

    ! Fill table with data
    call CFG_get(cfg, "td_mobility_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * mobility_fac
    call LT_set_col(ST_td_tbl, i_mobility, x_data, y_data)

    call CFG_get(cfg, "td_diffusion_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * diffusion_fac
    call LT_set_col(ST_td_tbl, i_diffusion, x_data, y_data)

    call CFG_get(cfg, "td_alpha_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * alpha_fac
    call LT_set_col(ST_td_tbl, i_alpha, x_data, y_data)

    call CFG_get(cfg, "td_eta_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * eta_fac
    call LT_set_col(ST_td_tbl, i_eta, x_data, y_data)

    ! Create table for photoionization
    if (ST_photoi_enabled) then
       call CFG_get(cfg, "photoi_rng_seed", rng_int4_seed)
       rng_int8_seed = transfer(rng_int4_seed, rng_int8_seed)
       call ST_rng%set_seed(rng_int8_seed)

       call photoi_get_table_air(ST_photoi_tbl, ST_photoi_frac_O2 * &
            ST_gas_pressure, 2 * ST_domain_len)
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

  !> Update the configuration file with the values read in 
  ! previous call of ST_read_config_files()
  subroutine ST_load_config()
    character(len=ST_slen)    :: tmp_name

    call CFG_get(cfg, "cylindrical", ST_cylindrical)
    call CFG_get(cfg, "end_time", ST_end_time)
    call CFG_get(cfg, "box_size", ST_box_size)
    call CFG_get(cfg, "output_dir", ST_output_dir)
    call CFG_get(cfg, "domain_len", ST_domain_len)
    call CFG_get(cfg, "applied_electric_fld_y", ST_applied_electric_fld_y)
    call CFG_get(cfg, "applied_electric_fld_x", ST_applied_electric_fld_x)
    call CFG_get(cfg, "dt_output", ST_dt_out)

    call CFG_get(cfg, "refine_per_steps", ST_refine_per_steps)
    call CFG_get(cfg, "refine_min_dx", ST_refine_min_dx)
    call CFG_get(cfg, "refine_max_dx", ST_refine_max_dx)
    call CFG_get(cfg, "refine_adx", ST_refine_adx)
    call CFG_get(cfg, "refine_cphi", ST_refine_cphi)
    call CFG_get(cfg, "derefine_dx", ST_derefine_dx)

    call CFG_get(cfg, "refine_init_time", ST_refine_init_time)
    call CFG_get(cfg, "refine_init_fac", ST_refine_init_fac)

    call CFG_get(cfg, "dt_max", ST_dt_max)
    call CFG_get(cfg, "epsilon_diel", ST_epsilon_diel)

    call CFG_get(cfg, "gas_pressure", ST_gas_pressure)
    call CFG_get(cfg, "photoi_enabled", ST_photoi_enabled)
    call CFG_get(cfg, "photoi_frac_O2", ST_photoi_frac_O2)
    call CFG_get(cfg, "photoi_eta", ST_photoi_eta)
    call CFG_get(cfg, "photoi_num_photons", ST_photoi_num_photons)

  end subroutine ST_load_config

  !> Compute the electric field at a given time
  function ST_get_electric_field(time) result(electric_fld)
    use m_units_constants

    real(dp), intent(in) :: time
    real(dp)             :: electric_fld, t_rel

    t_rel = time - ST_electric_fld_y_mod_t0
    if (t_rel > 0) then
       electric_fld = ST_applied_electric_fld_y * exp(-t_rel/ST_electric_fld_y_decay) + &
                      t_rel * ST_electric_fld_y_lin_deriv + &
                      ST_electric_fld_y_sin_amplitude * &
                      sin(t_rel * ST_electric_fld_y_sin_freq * 2 * UC_pi)
    else
       electric_fld = ST_applied_electric_fld_y
    end if
  end function ST_get_electric_field

  function ST_get_min(a, b, n) result(min_vec)
    integer, intent(in)  :: n
    real(dp), intent(in) :: a(n), b(n)
    real(dp)             :: min_vec(n)

    min_vec = min(a, b)
  end function ST_get_min

  !> Compute the voltage at a given time
  subroutine ST_set_voltage(time)
    use m_units_constants

    real(dp), intent(in) :: time
    real(dp)             :: electric_fld, t_rel

    t_rel = time - ST_electric_fld_y_mod_t0
    if (t_rel > 0) then
       electric_fld = ST_applied_electric_fld_y * exp(-t_rel/ST_electric_fld_y_decay) + &
                      t_rel * ST_electric_fld_y_lin_deriv + &
                      ST_electric_fld_y_sin_amplitude * &
                      sin(t_rel * ST_electric_fld_y_sin_freq * 2 * UC_pi)
    else
       electric_fld = ST_applied_electric_fld_y
    end if
    ST_applied_voltage = -ST_domain_len * electric_fld

  end subroutine ST_set_voltage

end module m_streamer
