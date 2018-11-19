!> This module contains several pre-defined variables like:
!! * Indices of cell-centered variables
!! * Names of the cell-centered variables
!! * Indices of face-centered variables
!! * Indices of transport data
module m_streamer
  use m_types
  use m_af_all
  use m_config
  use m_random
  use m_lookup_table

  implicit none
  public

  ! ** Indices of cell-centered variables **
  integer, protected :: n_var_cell     = 0  ! Number of variables
  integer, protected :: i_phi          = -1 ! Electrical potential
  integer, protected :: i_electron     = -1 ! Electron density
  integer, protected :: i_1pos_ion     = -1 ! First positive ion species
  integer, protected :: i_electric_fld = -1 ! Electric field norm
  integer, protected :: i_rhs          = -1 ! Source term Poisson
  integer, protected :: i_tmp          = -1 ! Temporary variable

  ! Optional variable (when using photoionization)
  integer :: i_photo = -1 ! Photoionization rate

  ! Optional variable (to show ionization source term)
  integer :: i_src = -1 ! Source term

  integer, parameter :: name_len = 12

  ! If true, only include n_e, n_i and |E| in output files
  logical :: ST_small_output = .false.

  ! ** Indices of face-centered variables **
  integer, protected :: n_var_face   = 0 ! Number of variables
  integer, protected :: flux_elec    = -1 ! Electron flux
  integer, protected :: flux_ion     = -1 ! Positive ion flux
  integer, protected :: electric_fld = -1 ! Electric field vector

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

  !> Maximal electric field for the lookup table
  real(dp), protected :: ST_max_field = 3e7_dp

  !> The diffusion coefficient for positive ions (m2/s)
  real(dp) :: ST_ion_diffusion = 0.0_dp

  !> The mobility of positive ions (m2/Vs)
  real(dp) :: ST_ion_mobility = 0.0_dp

  !> Whether to update ions (depends on ion diffusion/mobility)
  logical :: ST_update_ions = .false.

  ! Random number generator
  type(rng_t) :: ST_rng

  ! Parallel random number generator
  type(prng_t) :: ST_prng

  ! Name of the simulations
  character(len=ST_slen), protected :: ST_simulation_name = "sim"

  ! Output directory
  character(len=ST_slen), protected :: ST_output_dir = "output"

  ! Print status every this many seconds
  real(dp), protected :: ST_print_status_sec = 60.0_dp

  ! The refinement buffer width in cells (around flagged cells)
  integer, protected :: ST_refine_buffer_width = 4

  ! The number of steps after which the mesh is updated
  integer, protected :: ST_refine_per_steps = 2

  ! The grid spacing will always be larger than this value
  real(dp), protected :: ST_refine_min_dx = 1.0e-7_dp

  ! The grid spacing will always be smaller than this value
  real(dp), protected :: ST_refine_max_dx = 1.0e-3_dp

  ! Refine if alpha*dx is larger than this value
  real(dp), protected :: ST_refine_adx = 1.0_dp

  ! For refinement, use alpha(f * E)/f, where f is this factor
  real(dp), protected :: ST_refine_adx_fac = 1.0_dp

  ! Refine if the curvature in phi is larger than this value
  real(dp), protected :: ST_refine_cphi = 1e99_dp

  ! Allow derefinement if the curvature in phi is smaller than this value
  real(dp), protected :: ST_derefine_cphi = 1e99_dp

  ! Only derefine if grid spacing if smaller than this value
  real(dp), protected :: ST_derefine_dx = 1e-4_dp

  ! Refine around initial conditions up to this time
  real(dp), protected :: ST_refine_init_time = 10e-9_dp

  ! Refine until dx is smaller than this factor times the seed width
  real(dp), protected :: ST_refine_init_fac = 0.25_dp

  ! Minimum electron density for adding grid refinement
  real(dp), protected :: ST_refine_min_dens = -1.0e99_dp

  ! Refine a region up to this grid spacing
  real(dp), protected, allocatable :: ST_refine_regions_dr(:)

  ! Refine regions up to this simulation time
  real(dp), protected, allocatable :: ST_refine_regions_tstop(:)

  ! Minimum coordinate of the refinement regions
  real(dp), protected, allocatable :: ST_refine_regions_rmin(:,:)

  ! Maximum coordinate of the refinement regions
  real(dp), protected, allocatable :: ST_refine_regions_rmax(:,:)

  ! Current time step
  real(dp) :: ST_dt

  logical, protected :: ST_drt_limit_flux = .false.

  real(dp), protected :: ST_drt_limit_flux_Emax = 1e10_dp

  real(dp), protected :: ST_elec_density_threshold = -1.0e100_dp
  real(dp), protected :: ST_src_min_density = -1.0e100_dp
  real(dp), protected :: ST_src_max_density = 1.0e100_dp
  real(dp), protected :: ST_src_min_count = -1.0_dp
  real(dp), protected :: ST_src_max_field = -1.0_dp
  real(dp), protected :: ST_max_src = -1.0_dp
  real(dp), protected :: ST_velocity_max_field = -1.0_dp
  real(dp), protected :: ST_max_velocity = -1.0_dp
  real(dp), protected :: ST_diffusion_field_limit = 1.0e100_dp

  ! Number of time step restrictions
  integer, parameter :: ST_dt_num_cond = 3

  ! Array of time step restrictions per thread
  real(dp), allocatable :: ST_dt_matrix(:, :)

  ! Index of CFL condition
  integer, parameter :: ST_ix_cfl = 1

  ! Index of diffusion time step condition
  integer, parameter :: ST_ix_diff = 2

  ! Index of dielectric relaxation time condition
  integer, parameter :: ST_ix_drt = 3

  ! Safety factor for the time step
  real(dp) :: ST_dt_safety_factor = 0.9_dp

  ! Maximum allowed time step
  real(dp), protected :: ST_dt_max = 1.0e-11_dp

  ! Minimum allowed time step
  real(dp), protected :: ST_dt_min = 1.0e-14_dp

  ! Time between writing output
  real(dp), protected :: ST_dt_output = 1.0e-10_dp

  ! Include ionization source term in output
  logical, protected :: ST_output_src_term = .false.

  ! If positive: decay rate for source term (1/s) for time-averaged values
  real(dp), protected :: ST_output_src_decay_rate = -1.0_dp

  ! Write output along a line
  logical, protected :: ST_lineout_write = .false.

  ! Write Silo output
  logical, protected :: ST_silo_write = .true.

  ! Write binary output files (to resume later)
  logical, protected :: ST_datfile_write = .false.

  ! Use this many points for lineout data
  integer, protected :: ST_lineout_npoints = 500

  ! Relative position of line minimum coordinate
  real(dp), protected :: ST_lineout_rmin(3) = 0.0_dp

  ! Relative position of line maximum coordinate
  real(dp), protected :: ST_lineout_rmax(3) = 1.0_dp

  ! Write uniform output in a plane
  logical, protected :: ST_plane_write = .false.

  ! Which variable to include in plane
  integer, protected :: ST_plane_ivar

  ! Use this many points for plane data
  integer, protected :: ST_plane_npixels(2) = [64, 64]

  ! Relative position of plane minimum coordinate
  real(dp), protected :: ST_plane_rmin(3) = 0.0_dp

  ! Relative position of plane maximum coordinate
  real(dp), protected :: ST_plane_rmax(3) = 1.0_dp

  ! Current time
  real(dp)  :: ST_time

  ! End time of the simulation
  real(dp), protected :: ST_end_time = 10e-9_dp

  ! The size of the boxes that we use to construct our mesh
  integer, protected :: ST_box_size = 8

  ! The length of the (square) domain
  real(dp), protected :: ST_domain_len = 16e-3_dp

  ! Whether the domain is periodic in the x (and y)
  logical, protected :: ST_periodic = .false.

  ! Pressure of the gas in bar
  real(dp), protected :: ST_gas_pressure = 1.0_dp

  ! Fraction of O2
  real(dp), protected :: ST_gas_frac_O2 = 0.2_dp

  ! Number of V-cycles to perform per time step
  integer, protected :: ST_multigrid_num_vcycles = 2

  ! Used for the extrapolation of the ionization coefficient
  real(dp), protected :: ST_extrap_src_dydx
  real(dp), protected :: ST_extrap_src_y0

contains

  !> Create the configuration file with default values
  subroutine ST_initialize(tree, cfg, ndim)
    use iso_fortran_env, only: int64
    use m_config
    use omp_lib
    use m_chemistry
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg  !< The configuration for the simulation
    integer, intent(in)        :: ndim !< Number of dimensions
    integer                    :: n, n_threads
    real(dp)                   :: vec(ndim), tmp
    real(dp), allocatable      :: dbuffer(:)
    character(len=name_len)    :: varname
    integer                    :: rng_int4_seed(4) = &
         [8123, 91234, 12399, 293434]
    integer(int64)             :: rng_int8_seed(2)

    call CFG_add_get(cfg, "small_output", ST_small_output, &
         "If true, only include n_e, n_i and |E| in output files")

    ! Set index of electrons
    i_electron = af_find_cc_variable(tree, "e")

    ! Set index of first positive ion species
    do n = 1, n_species
       if (species_charge(n) == 1) then
          i_1pos_ion = n
          exit
       end if
    end do

    if (i_1pos_ion == -1) error stop "No positive ion species (1+) found"

    call af_add_cc_variable(tree, "phi", write_out=(.not. ST_small_output), ix=i_phi)
    call af_add_cc_variable(tree, "electric_fld", ix=i_electric_fld)
    call af_add_cc_variable(tree, "rhs", write_out=.false., ix=i_rhs)
    call af_add_cc_variable(tree, "tmp", write_out=.false., ix=i_tmp)

    call af_add_fc_variable(tree, "flux_elec", ix=flux_elec)
    call af_add_fc_variable(tree, "field", ix=electric_fld)

    call CFG_add_get(cfg, "ion_mobility", ST_ion_mobility, &
         "The mobility of positive ions (m2/Vs)")
    call CFG_add_get(cfg, "ion_diffusion", ST_ion_diffusion, &
         "The diffusion coefficient for positive ions (m2/s)")

    ST_update_ions = (abs(ST_ion_mobility) > 0 .or. abs(ST_ion_diffusion) > 0)
    if (ST_update_ions) then
       call af_add_fc_variable(tree, "flux_ion", ix=flux_ion)
    end if

    n_threads = af_get_max_threads()
    ! Prevent cache invalidation issues by enlarging the array
    allocate(ST_dt_matrix(ST_dt_num_cond+32, n_threads))

    call CFG_add_get(cfg, "cylindrical", ST_cylindrical, &
         "Whether cylindrical coordinates are used (only in 2D)")
    call CFG_add_get(cfg, "end_time", ST_end_time, &
         "The desired endtime (s) of the simulation")
    call CFG_add_get(cfg, "simulation_name", ST_simulation_name, &
         "The name of the simulation")
    call CFG_add_get(cfg, "output_dir", ST_output_dir, &
         "Directory where the output should be written")
    call CFG_add_get(cfg, "print_status_sec", ST_print_status_sec, &
         "Print status every this many seconds")
    call CFG_add_get(cfg, "box_size", ST_box_size, &
         "The number of grid cells per coordinate in a box")
    call CFG_add_get(cfg, "domain_len", ST_domain_len, &
         "The length of the domain (m)")
    call CFG_add_get(cfg, "periodic", ST_periodic, &
         "Whether the domain is periodic")
    call CFG_add_get(cfg, "gas_pressure", ST_gas_pressure, &
         "The gas pressure (bar), used for photoionization")
    call CFG_add_get(cfg, "gas_frac_O2", ST_gas_frac_O2, &
         "Fraction of O2, used for photoionization")

    call CFG_add_get(cfg, "multigrid_num_vcycles", ST_multigrid_num_vcycles, &
         "Number of V-cycles to perform per time step")

    call CFG_add_get(cfg, "fixes%drt_limit_flux", ST_drt_limit_flux, &
         "Avoid dielectric relaxation time step constraint by limiting flux")
    call CFG_add_get(cfg, "fixes%drt_limit_flux_Emax", ST_drt_limit_flux_Emax, &
         "Compute electron impact ionization source based on fluxes")
    call CFG_add_get(cfg, "fixes%velocity_max_field", ST_velocity_max_field, &
         "If positive, limit velocities to the value in this field (V/m)")
    call CFG_add_get(cfg, "fixes%src_max_field", ST_src_max_field, &
         "If positive, limit ionization to the value in this field (V/m)")
    call CFG_add_get(cfg, "fixes%src_min_density", ST_src_min_density, &
         "Disable impact ionization source below this density")
    call CFG_add_get(cfg, "fixes%src_max_density", ST_src_max_density, &
         "Disable impact ionization source above this density")
    call CFG_add_get(cfg, "fixes%src_min_count", ST_src_min_count, &
         "Disable impact ionization source below this electron count")
    call CFG_add_get(cfg, "fixes%elec_density_threshold", ST_elec_density_threshold, &
         "Set densities below this value to zero (1/m3)")
    call CFG_add_get(cfg, "fixes%diffusion_field_limit", ST_diffusion_field_limit, &
         "Disable diffusion parallel to fields above this threshold (V/m)")

    call CFG_add_get(cfg, "dt_output", ST_dt_output, &
         "The timestep for writing output (s)")
    call CFG_add_get(cfg, "dt_max", ST_dt_max, &
         "The maximum timestep (s)")
    call CFG_add_get(cfg, "dt_min", ST_dt_min, &
         "The minimum timestep (s)")
    call CFG_add_get(cfg, "dt_safety_factor", ST_dt_safety_factor, &
         "Safety factor for the time step")

    call CFG_add_get(cfg, "output_src_term", ST_output_src_term, &
         "Include ionization source term in output")

    tmp = 1/ST_output_src_decay_rate
    call CFG_add_get(cfg, "output_src_decay_time", tmp, &
         "If positive: decay time for source term (s) for time-averaged values")
    ST_output_src_decay_rate = 1/tmp

    if (ST_output_src_term) then
       call af_add_cc_variable(tree, "src", ix=i_src)
    end if

    call CFG_add_get(cfg, "lineout%write", ST_lineout_write, &
         "Write output along a line")
    call CFG_add_get(cfg, "silo_write", ST_silo_write, &
         "Write silo output")
    call CFG_add_get(cfg, "datfile_write", ST_datfile_write, &
         "Write binary output files (to resume later)")
    call CFG_add_get(cfg, "lineout%npoints", ST_lineout_npoints, &
         "Use this many points for lineout data")
    call CFG_add_get(cfg, "lineout%rmin", ST_lineout_rmin(1:ndim), &
         "Relative position of line minimum coordinate")
    call CFG_add_get(cfg, "lineout%rmax", ST_lineout_rmax(1:ndim), &
         "Relative position of line maximum coordinate")

    call CFG_add_get(cfg, "plane%write", ST_plane_write, &
         "Write uniform output in a plane")
    varname = "electron"
    call CFG_add_get(cfg, "plane%varname", varname, &
         "Names of variable to write in a plane")
    ST_plane_ivar = af_find_cc_variable(tree, trim(varname))
    call CFG_add_get(cfg, "plane%npixels", ST_plane_npixels, &
         "Use this many pixels for plane data")
    call CFG_add_get(cfg, "plane%rmin", ST_plane_rmin(1:ndim), &
         "Relative position of plane minimum coordinate")
    call CFG_add_get(cfg, "plane%rmax", ST_plane_rmax(1:ndim), &
         "Relative position of plane maximum coordinate")

    call CFG_add_get(cfg, "refine_buffer_width", ST_refine_buffer_width, &
         "The refinement buffer width in cells (around flagged cells)")
    call CFG_add_get(cfg, "refine_per_steps", ST_refine_per_steps, &
         "The number of steps after which the mesh is updated")
    call CFG_add_get(cfg, "refine_min_dx", ST_refine_min_dx, &
         "The grid spacing will always be larger than this value")
    call CFG_add_get(cfg, "refine_max_dx", ST_refine_max_dx, &
         "The grid spacing will always be smaller than this value")

    if (ST_refine_min_dx > ST_refine_max_dx) &
         error stop "Cannot have refine_min_dx < refine_max_dx"

    call CFG_add_get(cfg, "refine_adx", ST_refine_adx, &
         "Refine if alpha*dx is larger than this value")
    call CFG_add_get(cfg, "refine_adx_fac", ST_refine_adx_fac, &
         "For refinement, use alpha(f * E)/f, where f is this factor")
    call CFG_add_get(cfg, "refine_cphi", ST_refine_cphi, &
         "Refine if the curvature in phi is larger than this value")
    call CFG_add_get(cfg, "derefine_cphi", ST_derefine_cphi, &
         "Allow derefinement if the curvature in phi is smaller than this value")
    call CFG_add_get(cfg, "derefine_dx", ST_derefine_dx, &
         "Only derefine if grid spacing if smaller than this value")
    call CFG_add_get(cfg, "refine_init_time", ST_refine_init_time, &
         "Refine around initial conditions up to this time")
    call CFG_add_get(cfg, "refine_init_fac", ST_refine_init_fac, &
         "Refine until dx is smaller than this factor times the seed width")
    call CFG_add_get(cfg, "refine_min_dens", ST_refine_min_dens, &
         "Minimum electron density for adding grid refinement")

    call CFG_add(cfg, "refine_regions_dr", [1.0e99_dp], &
         "Refine regions up to this grid spacing", .true.)
    call CFG_add(cfg, "refine_regions_tstop", [-1.0e99_dp], &
         "Refine regions up to this simulation time", .true.)
    vec = 0.0_dp
    call CFG_add(cfg, "refine_regions_rmin", vec, &
         "Minimum coordinate of the refinement regions", .true.)
    call CFG_add(cfg, "refine_regions_rmax", vec, &
         "Maximum coordinate of the refinement regions", .true.)

    call CFG_get_size(cfg, "refine_regions_dr", n)
    allocate(ST_refine_regions_dr(n))
    allocate(ST_refine_regions_tstop(n))
    allocate(ST_refine_regions_rmin(ndim, n))
    allocate(ST_refine_regions_rmax(ndim, n))
    allocate(dbuffer(ndim * n))

    call CFG_get(cfg, "refine_regions_dr", ST_refine_regions_dr)
    call CFG_get(cfg, "refine_regions_tstop", ST_refine_regions_tstop)
    call CFG_get(cfg, "refine_regions_rmin", dbuffer)
    ST_refine_regions_rmin = reshape(dbuffer, [ndim, n])
    call CFG_get(cfg, "refine_regions_rmax", dbuffer)
    ST_refine_regions_rmax = reshape(dbuffer, [ndim, n])

    call CFG_add_get(cfg, "rng_seed", rng_int4_seed, &
         "Seed for random numbers. If all zero, generate from clock.")

    if (all(rng_int4_seed == 0)) then
       rng_int4_seed = get_random_seed()
       print *, "RNG seed: ", rng_int4_seed
    end if

    rng_int8_seed = transfer(rng_int4_seed, rng_int8_seed)
    call ST_rng%set_seed(rng_int8_seed)
    call ST_prng%init_parallel(n_threads, ST_rng)

  end subroutine ST_initialize

  !> Get a random seed based on the current time
  function get_random_seed() result(seed)
    integer :: seed(4)
    integer :: time, i

    call system_clock(time)
    do i = 1, 4
       seed(i) = ishftc(time, i*8)
    end do
  end function get_random_seed

  !> Initialize the transport coefficients
  subroutine ST_load_transport_data(cfg)
    use m_transport_data
    use m_config

    type(CFG_t), intent(inout) :: cfg
    character(len=ST_slen)     :: td_file       = "td_input_file.txt"
    character(len=ST_slen)     :: gas_name      = "AIR"
    integer                    :: table_size    = 1000
    real(dp)                   :: alpha_fac     = 1.0_dp
    real(dp)                   :: eta_fac       = 1.0_dp
    real(dp)                   :: mobility_fac  = 1.0_dp
    real(dp)                   :: diffusion_fac = 1.0_dp
    real(dp)                   :: x(2), y(2)
    real(dp), allocatable      :: x_data(:), y_data(:)
    character(len=ST_slen)     :: data_name

    call CFG_add_get(cfg, "transport_data_file", td_file, &
         "Input file with transport data")
    call CFG_add_get(cfg, "gas_name", gas_name, &
         "The name of the gas mixture used")

    call CFG_add_get(cfg, "lookup_table_size", table_size, &
         "The transport data table size in the fluid model")
    call CFG_add_get(cfg, "lookup_table_max_electric_fld", ST_max_field, &
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
    ST_td_tbl = LT_create(0.0_dp, ST_max_field, table_size, n_var_td)

    ! Fill table with data
    data_name = "efield[V/m]_vs_mu[m2/Vs]"
    call CFG_add_get(cfg, "td_mobility_name", data_name, &
         "The name of the mobility coefficient")
    call TD_get_from_file(td_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * mobility_fac
    call LT_set_col(ST_td_tbl, i_mobility, x_data, y_data)

    data_name = "efield[V/m]_vs_dif[m2/s]"
    call CFG_add_get(cfg, "td_diffusion_name", data_name, &
         "The name of the diffusion coefficient")
    call TD_get_from_file(td_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * diffusion_fac
    call LT_set_col(ST_td_tbl, i_diffusion, x_data, y_data)

    data_name = "efield[V/m]_vs_alpha[1/m]"
    call CFG_add_get(cfg, "td_alpha_name", data_name, &
         "The name of the eff. ionization coeff.")
    call TD_get_from_file(td_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * alpha_fac
    call LT_set_col(ST_td_tbl, i_alpha, x_data, y_data)

    data_name = "efield[V/m]_vs_eta[1/m]"
    call CFG_add_get(cfg, "td_eta_name", data_name, &
         "The name of the eff. attachment coeff.")
    call TD_get_from_file(td_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * eta_fac
    call LT_set_col(ST_td_tbl, i_eta, x_data, y_data)

    ! Determine extrapolation coefficients for transport data
    x = [ST_max_field, 0.8_dp * ST_max_field]

    ! Linear extrapolation for ionization coefficient
    y = LT_get_col(ST_td_tbl, i_alpha, x)

    ST_extrap_src_dydx = (y(2) - y(1)) / (x(2) - x(1))
    ST_extrap_src_y0   = y(2)

    ! Set maximal velocity if velocity_max_field is given
    if (ST_velocity_max_field > 0) then
       ST_max_velocity = ST_velocity_max_field * &
            LT_get_col(ST_td_tbl, i_mobility, ST_velocity_max_field)
    end if

    if (ST_src_max_field > 0) then
       ST_max_src = ST_src_max_field * &
            LT_get_col(ST_td_tbl, i_mobility, ST_src_max_field) * &
            LT_get_col(ST_td_tbl, i_alpha, ST_src_max_field)
    end if

  end subroutine ST_load_transport_data

end module m_streamer
