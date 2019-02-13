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
  integer, public, protected :: ix_electron    = -1 ! Electron density (in species list)
  integer, public, protected :: i_1pos_ion     = -1 ! First positive ion species
  integer, public, protected :: ix_1pos_ion    = -1 ! First positive ion (in species list)
  integer, public, protected :: i_electric_fld = -1 ! Electric field norm
  integer, public, protected :: i_rhs          = -1 ! Source term Poisson
  integer, public, protected :: i_tmp          = -1 ! Temporary variable
  integer, public, protected :: i_eps          = -1 ! Can be set to include a dielectric
  integer, public, protected :: i_gas_dens     = -1 ! Set for a variable gas density
  integer, public, protected, allocatable :: i_gas_comp(:) ! Set for a variable gas density

  ! ** Indices of face-centered variables **
  integer, public, protected :: n_var_face   = 0 ! Number of variables
  integer, public, protected :: flux_elec    = -1 ! Electron flux
  integer, public, protected :: electric_fld = -1 ! Electric field vector
  integer, public, protected :: flux_species(1) = -1

  ! Whether cylindrical coordinates are used
  logical, public, protected :: ST_cylindrical = .false.

  !> Whether a dielectric is used
  logical, public, protected :: ST_use_dielectric = .false.

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

  logical, public, protected :: ST_drt_limit_flux = .false.

  real(dp), public, protected :: ST_src_max_density = 1.0e100_dp
  real(dp), public, protected :: ST_max_velocity = -1.0_dp
  real(dp), public, protected :: ST_diffusion_field_limit = 1.0e100_dp

  ! End time of the simulation
  real(dp), public, protected :: ST_end_time = 10e-9_dp

  ! The size of the boxes that we use to construct our mesh
  integer, public, protected :: ST_box_size = 8

  ! Size of the coarse grid
  integer, public, protected :: ST_coarse_grid_size(NDIM) = 8

  ! Domain length per dimension
  real(dp), public, protected :: ST_domain_len(NDIM) = 16e-3_dp

  ! Whether the domain is periodic (per dimension)
  logical, public, protected :: ST_periodic(NDIM) = .false.

  ! Number of V-cycles to perform per time step
  integer, public, protected :: ST_multigrid_num_vcycles = 2

  ! Global time
  real(dp), public :: global_time = 0.0_dp

  ! Method used to prolong (interpolate) densities
  procedure(af_subr_prolong), pointer, public, protected :: &
       ST_prolongation_method => null()

  public :: ST_initialize

contains

  !> Create the configuration file with default values
  subroutine ST_initialize(tree, cfg, ndim)
    use iso_fortran_env, only: int64
    use m_config
    use omp_lib
    use m_chemistry
    use m_units_constants
    use m_gas
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg  !< The configuration for the simulation
    integer, intent(in)        :: ndim !< Number of dimensions
    integer                    :: n, n_threads
    character(len=name_len)    :: prolong_method
    integer                    :: rng_int4_seed(4) = &
         [8123, 91234, 12399, 293434]
    integer(int64)             :: rng_int8_seed(2)

    ! Set index of electrons
    i_electron = af_find_cc_variable(tree, "e")
    ix_electron = species_index("e")

    flux_species(1) = i_electron

    ! Set index of first positive ion species
    do n = 1, n_species
       if (species_charge(n) == 1) then
          i_1pos_ion = species_itree(n)
          ix_1pos_ion = n
          exit
       end if
    end do

    if (i_1pos_ion == -1) error stop "No positive ion species (1+) found"

    if (.not. gas_constant_density) then
       i_gas_dens = af_find_cc_variable(tree, "M")

       allocate(i_gas_comp(size(gas_components)))
       do n = 1, size(gas_components)
          i_gas_comp(n) = af_find_cc_variable(tree, trim(gas_components(n)))
       end do
    end if

    call af_add_cc_variable(tree, "phi", ix=i_phi)
    call af_add_cc_variable(tree, "electric_fld", ix=i_electric_fld)
    call af_add_cc_variable(tree, "rhs", ix=i_rhs)
    call af_add_cc_variable(tree, "tmp", write_out=.false., ix=i_tmp)

    call af_add_fc_variable(tree, "flux_elec", ix=flux_elec)
    call af_add_fc_variable(tree, "field", ix=electric_fld)

    call CFG_add_get(cfg, "cylindrical", ST_cylindrical, &
         "Whether cylindrical coordinates are used (only in 2D)")

    call CFG_add_get(cfg, "use_dielectric", ST_use_dielectric, &
         "Whether a dielectric is used (experimental)")
    if (ST_use_dielectric) then
       call af_add_cc_variable(tree, "eps", ix=i_eps)
       call af_set_cc_methods(tree, i_eps, af_bc_neumann_zero, &
            af_gc_prolong_copy, af_prolong_zeroth)
    end if

    call CFG_add_get(cfg, "end_time", ST_end_time, &
       "The desired endtime (s) of the simulation")
    call CFG_add_get(cfg, "box_size", ST_box_size, &
         "The number of grid cells per coordinate in a box")
    call CFG_add_get(cfg, "coarse_grid_size", ST_coarse_grid_size, &
         "The size of the coarse grid")
    call CFG_add_get(cfg, "domain_len", ST_domain_len, &
         "The length of the domain (m)")
    call CFG_add_get(cfg, "periodic", ST_periodic, &
         "Whether the domain is periodic (per dimension)")

    call CFG_add_get(cfg, "multigrid_num_vcycles", ST_multigrid_num_vcycles, &
         "Number of V-cycles to perform per time step")

    prolong_method = "limit"
    call CFG_add_get(cfg, "prolong_density", prolong_method, &
         "Density prolongation method (limit, linear, linear_cons, sparse)")
    select case (prolong_method)
    case ("limit")
       ST_prolongation_method => af_prolong_limit
    case ("linear")
       ST_prolongation_method => af_prolong_linear
    case ("linear_cons")
       ST_prolongation_method => af_prolong_linear_cons
    case ("sparse")
       ST_prolongation_method => af_prolong_sparse
    case ("zeroth")
       ST_prolongation_method => af_prolong_zeroth
    case default
       error stop "Unknown prolong_density method"
    end select

    call CFG_add_get(cfg, "fixes%drt_limit_flux", ST_drt_limit_flux, &
         "Avoid dielectric relaxation time step constraint by limiting flux")
    call CFG_add_get(cfg, "fixes%max_velocity", ST_max_velocity, &
         "Limit velocities to this value (m/s)")
    call CFG_add_get(cfg, "fixes%src_max_density", ST_src_max_density, &
         "Disable impact ionization source above this density")
    call CFG_add_get(cfg, "fixes%diffusion_field_limit", ST_diffusion_field_limit, &
         "Disable diffusion parallel to fields above this threshold (V/m)")

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

end module m_streamer
