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
  logical, protected :: ST_photoi_enabled = .false.

  ! Oxygen fraction
  real(dp), protected :: ST_photoi_frac_O2 = 0.2_dp

  ! Photoionization efficiency
  real(dp), protected :: ST_photoi_eta = 0.05_dp

  ! Number of photons to use
  integer, protected :: ST_photoi_num_photons = 100*1000

  ! Table for photoionization
  type(photoi_tbl_t), protected :: ST_photoi_tbl

  ! Name of the simulations
  character(len=ST_slen), protected :: ST_simulation_name = "sim"

  ! Output directory
  character(len=ST_slen), protected :: ST_output_dir = "output"

  ! The number of steps after which the mesh is updated
  integer, protected :: ST_refine_per_steps = 2

  ! The grid spacing will always be larger than this value
  real(dp), protected :: ST_refine_min_dx = 1.0e-6_dp

  ! The grid spacing will always be smaller than this value
  real(dp), protected :: ST_refine_max_dx = 1.0e-3_dp

  ! Refine if alpha*dx is larger than this value
  real(dp), protected :: ST_refine_adx = 1.0_dp

  ! Refine if the curvature in phi is larger than this value
  real(dp), protected :: ST_refine_cphi = 1e99_dp

  ! Only derefine if grid spacing if smaller than this value
  real(dp), protected :: ST_derefine_dx = 1e-4_dp

  ! Refine around initial conditions up to this time
  real(dp), protected :: ST_refine_init_time = 10e-9_dp

  ! Refine until dx is smaller than this factor times the seed width
  real(dp), protected :: ST_refine_init_fac = 0.25_dp

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
  real(dp), protected :: ST_dt_max = 1.0e-11_dp

  ! Time between writing output
  real(dp), protected :: ST_dt_output = 1.0e-10_dp

  ! Write output along a line
  logical, protected :: ST_lineout_write = .false.

  ! Use this many points for lineout data
  integer, protected :: ST_lineout_npoints = 500

  ! Relative position of line minimum coordinate
  real(dp) :: ST_lineout_rmin(3) = 0.0_dp

  ! Relative position of line maximum coordinate
  real(dp) :: ST_lineout_rmax(3) = 1.0_dp

  ! Current time
  real(dp)  :: ST_time

  ! End time of the simulation
  real(dp), protected :: ST_end_time = 10e-9_dp

  ! The size of the boxes that we use to construct our mesh
  integer, protected :: ST_box_size = 8

  ! The length of the (square) domain
  real(dp), protected :: ST_domain_len = 16e-3_dp

  ! Pressure of the gas in bar
  real(dp), protected :: ST_gas_pressure = 1.0_dp

contains

  !> Create the configuration file with default values
  subroutine ST_initialize(cfg)
    use m_config
    type(CFG_t), intent(inout) :: cfg !< The configuration for the simulation

    call CFG_add_get(cfg, "cylindrical", ST_cylindrical, &
         "Whether cylindrical coordinates are used (only in 2D)")
    call CFG_add_get(cfg, "end_time", ST_end_time, &
         "The desired endtime (s) of the simulation")
    call CFG_add_get(cfg, "simulation_name", ST_simulation_name, &
         "The name of the simulation")
    call CFG_add_get(cfg, "output_dir", ST_output_dir, &
         "Directory where the output should be written")
    call CFG_add_get(cfg, "box_size", ST_box_size, &
         "The number of grid cells per coordinate in a box")
    call CFG_add_get(cfg, "domain_len", ST_domain_len, &
         "The length of the domain (m)")
    call CFG_add_get(cfg, "gas_pressure", ST_gas_pressure, &
         "The gas pressure (bar), used for photoionization")

    call CFG_add_get(cfg, "dt_output", ST_dt_output, &
         "The timestep for writing output (s)")
    call CFG_add_get(cfg, "dt_max", ST_dt_max, &
         "The maximum timestep (s)")

    call CFG_add_get(cfg, "lineout_write", ST_lineout_write, &
         "Write output along a line")
    call CFG_add_get(cfg, "lineout_npoints", ST_lineout_npoints, &
         "Use this many points for lineout data")
    call CFG_add_get(cfg, "lineout_rmin", ST_lineout_rmin, &
         "Relative position of line minimum coordinate")
    call CFG_add_get(cfg, "lineout_rmax", ST_lineout_rmax, &
         "Relative position of line maximum coordinate")

    call CFG_add_get(cfg, "refine_per_steps", ST_refine_per_steps, &
         "The number of steps after which the mesh is updated")
    call CFG_add_get(cfg, "refine_min_dx", ST_refine_min_dx, &
         "The grid spacing will always be larger than this value")
    call CFG_add_get(cfg, "refine_max_dx", ST_refine_max_dx, &
         "The grid spacing will always be smaller than this value")

    call CFG_add_get(cfg, "refine_adx", ST_refine_adx, &
         "Refine if alpha*dx is larger than this value")
    call CFG_add_get(cfg, "refine_cphi", ST_refine_cphi, &
         "Refine if the curvature in phi is larger than this value")
    call CFG_add_get(cfg, "derefine_dx", ST_derefine_dx, &
         "Only derefine if grid spacing if smaller than this value")
    call CFG_add_get(cfg, "refine_init_time", ST_refine_init_time, &
         "Refine around initial conditions up to this time")
    call CFG_add_get(cfg, "refine_init_fac", ST_refine_init_fac, &
         "Refine until dx is smaller than this factor times the seed width")

    call CFG_add_get(cfg, "photoi_enabled", ST_photoi_enabled, &
         "Whether photoionization is enabled")
    call CFG_add_get(cfg, "photoi_frac_O2", ST_photoi_frac_O2, &
         "Fraction of oxygen (0-1)")
    call CFG_add_get(cfg, "photoi_eta", ST_photoi_eta, &
         "Photoionization efficiency factor, typically around 0.05")
    call CFG_add_get(cfg, "photoi_num_photons", ST_photoi_num_photons, &
         "Number of discrete photons to use for photoionization")

  end subroutine ST_initialize

  !> Initialize the transport coefficients
  subroutine ST_load_transport_data(cfg)
    use iso_fortran_env, only: int64
    use m_transport_data
    use m_config

    type(CFG_t), intent(inout) :: cfg

    character(len=ST_slen)     :: td_file = "td_input_file.txt"
    character(len=ST_slen)     :: gas_name         = "AIR"
    integer                    :: table_size       = 1000
    integer                    :: rng_int4_seed(4) = &
         [8123, 91234, 12399, 293434]
    integer(int64)             :: rng_int8_seed(2)
    real(dp)                   :: max_electric_fld = 3e7_dp
    real(dp)                   :: alpha_fac        = 1.0_dp
    real(dp)                   :: eta_fac          = 1.0_dp
    real(dp)                   :: mobility_fac     = 1.0_dp
    real(dp)                   :: diffusion_fac    = 1.0_dp
    real(dp), allocatable      :: x_data(:), y_data(:)
    character(len=ST_slen)     :: data_name

    call CFG_add_get(cfg, "transport_data_file", td_file, &
         "Input file with transport data")
    call CFG_add_get(cfg, "gas_name", gas_name, &
         "The name of the gas mixture used")

    call CFG_add_get(cfg, "lookup_table_size", table_size, &
         "The transport data table size in the fluid model")
    call CFG_add_get(cfg, "lookup_table_max_electric_fld", max_electric_fld, &
         "The maximum electric field in the fluid model coefficients")

    call CFG_add_get(cfg, "td_alpha_fac", alpha_fac, &
         "Modify alpha by this factor")
    call CFG_add_get(cfg, "td_eta_fac", eta_fac, &
         "Modify eta by this factor")
    call CFG_add_get(cfg, "td_mobility_fac", mobility_fac, &
         "Modify mobility by this factor")
    call CFG_add_get(cfg, "td_diffusion_fac", diffusion_fac, &
         "Modify diffusion by this factor")

    ! Create a lookup table for the model coefficients
    ST_td_tbl = LT_create(0.0_dp, max_electric_fld, table_size, n_var_td)

    ! Fill table with data
    data_name = "efield[V/m]_vs_mu[m2/Vs]"
    call CFG_add_get(cfg, "td_mobility_name", data_name, &
         "The name of the mobility coefficient")
    call TD_get_td_from_file(td_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * mobility_fac
    call LT_set_col(ST_td_tbl, i_mobility, x_data, y_data)

    data_name = "efield[V/m]_vs_dif[m2/s]"
    call CFG_add_get(cfg, "td_diffusion_name", data_name, &
         "The name of the diffusion coefficient")
    call TD_get_td_from_file(td_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * diffusion_fac
    call LT_set_col(ST_td_tbl, i_diffusion, x_data, y_data)

    data_name = "efield[V/m]_vs_alpha[1/m]"
    call CFG_add_get(cfg, "td_alpha_name", data_name, &
         "The name of the eff. ionization coeff.")
    call TD_get_td_from_file(td_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * alpha_fac
    call LT_set_col(ST_td_tbl, i_alpha, x_data, y_data)

    data_name = "efield[V/m]_vs_eta[1/m]"
    call CFG_add_get(cfg, "td_eta_name", data_name, &
         "The name of the eff. attachment coeff.")
    call TD_get_td_from_file(td_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * eta_fac
    call LT_set_col(ST_td_tbl, i_eta, x_data, y_data)

    ! Create table for photoionization
    if (ST_photoi_enabled) then
       call CFG_add_get(cfg, "photoi_rng_seed", rng_int4_seed, &
            "Seed for the photoionization random number generator")
       rng_int8_seed = transfer(rng_int4_seed, rng_int8_seed)
       call ST_rng%set_seed(rng_int8_seed)

       call photoi_get_table_air(ST_photoi_tbl, ST_photoi_frac_O2 * &
            ST_gas_pressure, 2 * ST_domain_len)
    end if

  end subroutine ST_load_transport_data

end module m_streamer
