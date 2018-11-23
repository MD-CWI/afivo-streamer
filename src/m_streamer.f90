!> This module contains several pre-defined variables like:
!! * Indices of cell-centered variables
!! * Names of the cell-centered variables
!! * Indices of face-centered variables
!! * Indices of transport data
module m_streamer
  use m_types
  use m_af_all
  use m_random
  use m_lookup_table

  implicit none
  private

  ! ** Indices of cell-centered variables **
  integer, public, protected :: n_var_cell     = 0  ! Number of variables
  integer, public, protected :: i_phi          = -1 ! Electrical potential
  integer, public, protected :: i_electron     = -1 ! Electron density
  integer, public, protected :: i_1pos_ion     = -1 ! First positive ion species
  integer, public, protected :: i_electric_fld = -1 ! Electric field norm
  integer, public, protected :: i_rhs          = -1 ! Source term Poisson
  integer, public, protected :: i_tmp          = -1 ! Temporary variable

  ! Optional variable (to show ionization source term)
  integer, public, protected :: i_src = -1 ! Source term

  ! If true, only include n_e, n_i and |E| in output files
  logical, public, protected :: ST_small_output = .false.

  ! ** Indices of face-centered variables **
  integer, public, protected :: n_var_face   = 0 ! Number of variables
  integer, public, protected :: flux_elec    = -1 ! Electron flux
  integer, public, protected :: electric_fld = -1 ! Electric field vector
  integer, public, protected :: flux_species(1) = -1

  ! ** Indices of transport data **
  integer, parameter, public :: n_var_td    = 3 ! Number of transport coefficients
  integer, parameter, public :: i_mobility  = 1 ! Electron mobility
  integer, parameter, public :: i_diffusion = 2 ! Electron diffusion constant
  integer, parameter, public :: i_alpha_eff = 3 ! Net ionization coefficient

  ! Whether cylindrical coordinates are used
  logical, public, protected :: ST_cylindrical = .false.

  ! Table with transport data vs electric field
  type(LT_t), public, protected :: ST_td_tbl

  !> Maximal electric field for the lookup table
  real(dp), public, protected :: ST_max_field = 3e7_dp

  !> The diffusion coefficient for positive ions (m2/s)
  real(dp), public, protected :: ST_ion_diffusion = 0.0_dp

  !> The mobility of positive ions (m2/Vs)
  real(dp), public, protected :: ST_ion_mobility = 0.0_dp

  !> Whether to update ions (depends on ion diffusion/mobility)
  logical, public, protected :: ST_update_ions = .false.

  ! Random number generator
  type(rng_t), public :: ST_rng

  ! Parallel random number generator
  type(prng_t), public :: ST_prng

  ! Name of the simulations
  character(len=string_len), public, protected :: ST_simulation_name = "sim"

  ! Output directory
  character(len=string_len), public, protected :: ST_output_dir = "output"

  ! Print status every this many seconds
  real(dp), public, protected :: ST_print_status_sec = 60.0_dp

  logical, public, protected :: ST_drt_limit_flux = .false.

  real(dp), public, protected :: ST_src_max_density = 1.0e100_dp
  real(dp), public, protected :: ST_max_velocity = -1.0_dp
  real(dp), public, protected :: ST_diffusion_field_limit = 1.0e100_dp

  ! Time between writing output
  real(dp), public, protected :: ST_dt_output = 1.0e-10_dp

  ! Include ionization source term in output
  logical, public, protected :: ST_output_src_term = .false.

  ! If positive: decay rate for source term (1/s) for time-averaged values
  real(dp), public, protected :: ST_output_src_decay_rate = -1.0_dp

  ! Write output along a line
  logical, public, protected :: ST_lineout_write = .false.

  ! Write Silo output
  logical, public, protected :: ST_silo_write = .true.

  ! Write binary output files (to resume later)
  logical, public, protected :: ST_datfile_write = .false.

  ! Use this many points for lineout data
  integer, public, protected :: ST_lineout_npoints = 500

  ! Relative position of line minimum coordinate
  real(dp), public, protected :: ST_lineout_rmin(3) = 0.0_dp

  ! Relative position of line maximum coordinate
  real(dp), public, protected :: ST_lineout_rmax(3) = 1.0_dp

  ! Write uniform output in a plane
  logical, public, protected :: ST_plane_write = .false.

  ! Which variable to include in plane
  integer, public, protected :: ST_plane_ivar

  ! Use this many points for plane data
  integer, public, protected :: ST_plane_npixels(2) = [64, 64]

  ! Relative position of plane minimum coordinate
  real(dp), public, protected :: ST_plane_rmin(3) = 0.0_dp

  ! Relative position of plane maximum coordinate
  real(dp), public, protected :: ST_plane_rmax(3) = 1.0_dp

  ! End time of the simulation
  real(dp), public, protected :: ST_end_time = 10e-9_dp

  ! The size of the boxes that we use to construct our mesh
  integer, public, protected :: ST_box_size = 8

  ! The length of the (square) domain
  real(dp), public, protected :: ST_domain_len = 16e-3_dp

  ! Whether the domain is periodic in the x (and y)
  logical, public, protected :: ST_periodic = .false.

  ! Number of V-cycles to perform per time step
  integer, public, protected :: ST_multigrid_num_vcycles = 2

  ! Global time
  real(dp), public :: global_time = 0.0_dp

  public :: ST_initialize
  public :: ST_load_transport_data

contains

  !> Create the configuration file with default values
  subroutine ST_initialize(tree, cfg, ndim)
    use iso_fortran_env, only: int64
    use m_config
    use omp_lib
    use m_chemistry
    use m_units_constants
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg  !< The configuration for the simulation
    integer, intent(in)        :: ndim !< Number of dimensions
    integer                    :: n, n_threads
    real(dp)                   :: tmp
    character(len=name_len)    :: varname
    integer                    :: rng_int4_seed(4) = &
         [8123, 91234, 12399, 293434]
    integer(int64)             :: rng_int8_seed(2)

    call CFG_add_get(cfg, "small_output", ST_small_output, &
         "If true, only include n_e, n_i and |E| in output files")

    ! Set index of electrons
    i_electron = af_find_cc_variable(tree, "e")

    flux_species(1) = i_electron

    ! Set index of first positive ion species
    do n = 1, n_species
       if (species_charge(n) == 1) then
          i_1pos_ion = species_ix(n)
          exit
       end if
    end do

    if (i_1pos_ion == -1) error stop "No positive ion species (1+) found"

    call af_add_cc_variable(tree, "phi", write_out=(.not. ST_small_output), ix=i_phi)
    call af_add_cc_variable(tree, "electric_fld", ix=i_electric_fld)
    call af_add_cc_variable(tree, "rhs", write_out=(.not. ST_small_output), ix=i_rhs)
    call af_add_cc_variable(tree, "tmp", write_out=.false., ix=i_tmp)

    call af_add_fc_variable(tree, "flux_elec", ix=flux_elec)
    call af_add_fc_variable(tree, "field", ix=electric_fld)

    ! call CFG_add_get(cfg, "ion_mobility", ST_ion_mobility, &
    !      "The mobility of positive ions (m2/Vs)")
    ! call CFG_add_get(cfg, "ion_diffusion", ST_ion_diffusion, &
    !      "The diffusion coefficient for positive ions (m2/s)")

    ! ST_update_ions = (abs(ST_ion_mobility) > 0 .or. abs(ST_ion_diffusion) > 0)
    ! if (ST_update_ions) then
    !    call af_add_fc_variable(tree, "flux_ion", ix=flux_ion)
    ! end if

    call CFG_add_get(cfg, "cylindrical", ST_cylindrical, &
         "Whether cylindrical coordinates are used (only in 2D)")
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

    call CFG_add_get(cfg, "multigrid_num_vcycles", ST_multigrid_num_vcycles, &
         "Number of V-cycles to perform per time step")

    call CFG_add_get(cfg, "fixes%drt_limit_flux", ST_drt_limit_flux, &
         "Avoid dielectric relaxation time step constraint by limiting flux")
    call CFG_add_get(cfg, "fixes%max_velocity", ST_max_velocity, &
         "Limit velocities to this value (m/s)")
    call CFG_add_get(cfg, "fixes%src_max_density", ST_src_max_density, &
         "Disable impact ionization source above this density")
    call CFG_add_get(cfg, "fixes%diffusion_field_limit", ST_diffusion_field_limit, &
         "Disable diffusion parallel to fields above this threshold (V/m)")

    call CFG_add_get(cfg, "dt_output", ST_dt_output, &
         "The timestep for writing output (s)")
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
    varname = "e"
    call CFG_add_get(cfg, "plane%varname", varname, &
         "Names of variable to write in a plane")
    ST_plane_ivar = af_find_cc_variable(tree, trim(varname))
    call CFG_add_get(cfg, "plane%npixels", ST_plane_npixels, &
         "Use this many pixels for plane data")
    call CFG_add_get(cfg, "plane%rmin", ST_plane_rmin(1:ndim), &
         "Relative position of plane minimum coordinate")
    call CFG_add_get(cfg, "plane%rmax", ST_plane_rmax(1:ndim), &
         "Relative position of plane maximum coordinate")

    call CFG_add_get(cfg, "rng_seed", rng_int4_seed, &
         "Seed for random numbers. If all zero, generate from clock.")

    if (all(rng_int4_seed == 0)) then
       rng_int4_seed = get_random_seed()
       print *, "RNG seed: ", rng_int4_seed
    end if

    rng_int8_seed = transfer(rng_int4_seed, rng_int8_seed)
    call ST_rng%set_seed(rng_int8_seed)
    n_threads = af_get_max_threads()
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
    use m_gas
    use m_chemistry

    type(CFG_t), intent(inout) :: cfg
    character(len=string_len)  :: td_file = "td_input_file.txt"
    real(dp), allocatable      :: x_data(:), y_data(:)

    ! Create a lookup table for the model coefficients
    ST_td_tbl = LT_create(rate_min_townsend, rate_max_townsend, &
         rate_table_size, n_var_td)

    call CFG_add_get(cfg, "transport_data_file", td_file, &
         "Input file with transport data")

    ! Fill table with data
    call TD_get_from_file(td_file, "Mobility *N (1/m/V/s)", x_data, y_data)
    call LT_set_col(ST_td_tbl, i_mobility, x_data, y_data)

    call TD_get_from_file(td_file, "Diffusion coefficient *N (1/m/s)", &
         x_data, y_data)
    call LT_set_col(ST_td_tbl, i_diffusion, x_data, y_data)

    call TD_get_from_file(td_file, "Townsend ioniz. coef. alpha/N (m2)", &
         x_data, y_data)
    call LT_set_col(ST_td_tbl, i_alpha_eff, x_data, y_data)

  end subroutine ST_load_transport_data

end module m_streamer
