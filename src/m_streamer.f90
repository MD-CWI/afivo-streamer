#include "../afivo/src/cpp_macros.h"
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

  !> Index of electrical potential
  integer, public, protected :: i_phi          = -1
  !> Index of electron density
  integer, public, protected :: i_electron     = -1
  !> Index of electron density (in species list)
  integer, public, protected :: ix_electron    = -1
  !> Index of first positive ion species
  integer, public, protected :: i_1pos_ion     = -1
  !> Index of first positive ion (in species list)
  integer, public, protected :: ix_1pos_ion    = -1
  !> Index of electric field norm
  integer, public, protected :: i_electric_fld = -1
  !> Index of source term Poisson
  integer, public, protected :: i_rhs          = -1
  !> Index of temporary variable
  integer, public, protected :: i_tmp          = -1
  !> Index of can be set to include a dielectric
  integer, public, protected :: i_eps          = -1

  !> Index of electron flux
  integer, public, protected :: flux_elec    = -1
  !> Index of electric field vector
  integer, public, protected :: electric_fld = -1
  !> List of all flux variables (face-centered index)
  integer, public, protected, allocatable :: flux_variables(:)
  !> List of all flux species (cell-centered index)
  integer, public, protected, allocatable :: flux_species(:)
  !> List of the charges of the flux species
  integer, public, protected, allocatable :: flux_species_charge(:)
  !> List of positive ion fluxes (useful for secondary emission)
  integer, public, protected, allocatable :: flux_pos_ion(:)

  !> Whether cylindrical coordinates are used
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

  !> Random number generator
  type(rng_t), public :: ST_rng

  !> Parallel random number generator
  type(prng_t), public :: ST_prng

  !> Avoid dielectric relaxation time step constraint by limiting flux
  logical, public, protected :: ST_drt_limit_flux = .false.

  !> Limit velocities to this value (m/s)
  real(dp), public, protected :: ST_max_velocity = -1.0_dp

  !> Disable diffusion parallel to fields above this threshold (V/m)
  real(dp), public, protected :: ST_diffusion_field_limit = 1.0e100_dp

  !> End time of the simulation
  real(dp), public, protected :: ST_end_time = 10e-9_dp

  !> If we are using ST_end_time
  logical, public, protected :: ST_use_end_time = .true.

  !> Whether streamer length is used as a simulation stopping
  logical, public, protected :: ST_use_end_streamer_length = .false.

  !> Wait n steps before initializing streamer begin position
  integer, public, protected :: ST_initial_streamer_pos_steps_wait = 5

  !> Streamer length at which the simulation will stop
  real(dp), public, protected :: ST_end_streamer_length = 15e-3


  !> The size of the boxes that we use to construct our mesh
  integer, public, protected :: ST_box_size = 8

  !> Size of the coarse grid
  integer, public, protected :: ST_coarse_grid_size(NDIM) = 8

  !> Domain length per dimension
  real(dp), public, protected :: ST_domain_len(NDIM) = 16e-3_dp

  !> Origin of domain
  real(dp), public, protected :: ST_domain_origin(NDIM) = 0.0_dp

  !> Whether the domain is periodic (per dimension)
  logical, public, protected :: ST_periodic(NDIM) = .false.

  !> Number of V-cycles to perform per time step
  integer, public, protected :: ST_multigrid_num_vcycles = 2

  ! Stop multigrid when residual is smaller than this factor times max(|rhs|)
  real(dp), public, protected :: ST_multigrid_relative_residual = 1e-4_dp

  ! Global time
  real(dp), public :: global_time = 0.0_dp

  !> Method used to prolong (interpolate) densities
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
    use m_transport_data
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg  !< The configuration for the simulation
    integer, intent(in)        :: ndim !< Number of dimensions
    integer                    :: n, k, n_threads, ix_chemistry
    character(len=name_len)    :: prolong_method
    integer                    :: rng_int4_seed(4) = &
         [8123, 91234, 12399, 293434]
    integer(int64)             :: rng_int8_seed(2)

    ! Set index of electrons
    i_electron = af_find_cc_variable(tree, "e")
    ix_electron = species_index("e")

    ! Set index of first positive ion species
    do n = n_gas_species+1, n_species
       if (species_charge(n) == 1) then
          i_1pos_ion = species_itree(n)
          ix_1pos_ion = n
          exit
       end if
    end do

    if (i_1pos_ion == -1) error stop "No positive ion species (1+) found"

    ! Set flux species
    call af_add_fc_variable(tree, "flux_elec", ix=flux_elec, &
         write_binary=.false.)
    call af_add_fc_variable(tree, "field", ix=electric_fld)

    allocate(flux_species(1+transport_data_ions%n_mobile_ions))
    allocate(flux_species_charge(1+transport_data_ions%n_mobile_ions))
    allocate(flux_variables(1+transport_data_ions%n_mobile_ions))
    flux_species(1)        = i_electron
    flux_species_charge(1) = -1
    flux_variables(1)      = flux_elec

    do n = 1, transport_data_ions%n_mobile_ions
       flux_species(1+n) = af_find_cc_variable(tree, &
            trim(transport_data_ions%names(n)))

       ! Get index in chemistry list and determine charge
       ix_chemistry = species_index(trim(transport_data_ions%names(n)))
       flux_species_charge(1+n) = species_charge(ix_chemistry)

       call af_add_fc_variable(tree, trim(transport_data_ions%names(n)), &
            ix=flux_variables(1+n), write_binary=.false.)
    end do

    ! Create a list of positive ion fluxes for secondary emission
    n = count(flux_species_charge > 0)
    allocate(flux_pos_ion(n))

    k = 0
    do n = 1, size(flux_species_charge)
       if (flux_species_charge(n) > 0) then
          k = k + 1
          flux_pos_ion(k) = flux_variables(n)
       end if
    end do

    call af_add_cc_variable(tree, "phi", ix=i_phi)
    call af_add_cc_variable(tree, "electric_fld", ix=i_electric_fld)
    call af_add_cc_variable(tree, "rhs", ix=i_rhs)
    call af_add_cc_variable(tree, "tmp", write_out=.false., &
         write_binary=.false., ix=i_tmp)

    call CFG_add_get(cfg, "cylindrical", ST_cylindrical, &
         "Whether cylindrical coordinates are used (only in 2D)")

    call CFG_add_get(cfg, "use_dielectric", ST_use_dielectric, &
         "Whether a dielectric is used (experimental)")
    if (ST_use_dielectric) then
       call af_add_cc_variable(tree, "eps", ix=i_eps)
       call af_set_cc_methods(tree, i_eps, af_bc_neumann_zero, &
            af_gc_prolong_copy, af_prolong_zeroth)
    end if

    call CFG_add_get(cfg, "use_end_time", ST_use_end_time, &
	       "Whether end_time is used to end the simulation")
    call CFG_add_get(cfg, "use_end_streamer_length", ST_use_end_streamer_length, &
         "Whether the length of the streamer is used to end the simulation")
    call CFG_add_get(cfg, "end_streamer_length", ST_end_streamer_length, &
         "Streamer length at which the simulation will end.")
    call CFG_add_get(cfg, "initial_streamer_pos_steps_wait", ST_initial_streamer_pos_steps_wait, &
          "Number of simulation steps to wait before initializing the starting position of the streamer")

    call CFG_add_get(cfg, "end_time", ST_end_time, &
       "The desired endtime (s) of the simulation")
    call CFG_add_get(cfg, "box_size", ST_box_size, &
         "The number of grid cells per coordinate in a box")
    call CFG_add_get(cfg, "coarse_grid_size", ST_coarse_grid_size, &
         "The size of the coarse grid")
    call CFG_add_get(cfg, "domain_len", ST_domain_len, &
         "The length of the domain (m)")
    call CFG_add_get(cfg, "domain_origin", ST_domain_origin, &
         "The origin of the domain (m)")
    call CFG_add_get(cfg, "periodic", ST_periodic, &
         "Whether the domain is periodic (per dimension)")

    call CFG_add_get(cfg, "multigrid_num_vcycles", ST_multigrid_num_vcycles, &
         "Number of V-cycles to perform per time step")
    call CFG_add_get(cfg, "multigrid_relative_residual", &
         ST_multigrid_relative_residual, &
         "Stop multigrid when residual is smaller than this factor times max(|rhs|)")

    prolong_method = "limit_pos"
    call CFG_add_get(cfg, "prolong_density", prolong_method, &
         "Density prolongation method (limit, linear, linear_cons, sparse)")
    select case (prolong_method)
    case ("limit_pos")
       ST_prolongation_method => af_prolong_limit_pos
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
