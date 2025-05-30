!> This module contains several pre-defined variables like:
!! * Indices of cell-centered variables
!! * Names of the cell-centered variables
!! * Indices of face-centered variables
!! * Indices of transport data
module m_streamer
#include "../afivo/src/cpp_macros.h"
  use m_types
  use m_af_all
  use m_random
  use m_lookup_table
  use m_config

  implicit none
  private

  !> Index of electrical potential
  integer, public, protected :: i_phi          = -1
  !> Index of electron density
  integer, public, protected :: i_electron     = -1
  !> Index of electron density (in species list)
  integer, public, protected :: ix_electron    = -1
  !> Index of electron energy density
  integer, public, protected :: i_electron_energy = -1
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
  !> Index can be set to include a dielectric
  integer, public, protected :: i_eps          = -1
  !> Index can be set to include an electrode
  integer, public, protected :: i_lsf          = -1

  !> Index of all densities that evolve in time
  integer, public, protected, allocatable :: all_densities(:)

  !> Include deposited power density in output
  logical, public, protected :: compute_power_density = .false.
  !> Index of deposited power density
  integer, public, protected :: i_power_density = -1

  !> Index of correction factor for source terms
  integer, public, protected :: i_srcfac = -1

  !> Index of electron flux
  integer, public, protected :: flux_elec    = -1
  !> Index of electron energy flux
  integer, public, protected :: flux_energy  = -1
  !> Index of electric field vector
  integer, public, protected :: electric_fld = -1

  !> Number of flux variables
  integer, public, protected :: flux_num_species = -1
  !> Number of electron flux variables
  integer, public, protected :: flux_num_electron_vars = -1
  !> List of all flux variables (face-centered index)
  integer, public, protected, allocatable :: flux_variables(:)
  !> List of all flux species (cell-centered index)
  integer, public, protected, allocatable :: flux_species(:)
  !> List of the charges of the flux species
  integer, public, protected, allocatable :: flux_species_charge(:)
  !> List of the signs of the charges of the flux species (+- 1)
  integer, public, protected, allocatable :: flux_species_charge_sign(:)
  !> List of positive ion fluxes (useful for secondary emission)
  integer, public, protected, allocatable :: flux_pos_ion(:)

  !> Whether cylindrical coordinates are used
  logical, public, protected :: ST_cylindrical = .false.

  !> Whether a dielectric is used
  logical, public, protected :: ST_use_dielectric = .false.

  !> Whether to include an electrode
  logical, public, protected :: ST_use_electrode = .false.

  !> Boundary condition for the plasma species
  procedure(af_subr_bc), public, protected, pointer :: &
       bc_species => null()

  !> Multigrid option structure
  type(mg_t), public :: mg

  !> Random number generator
  type(rng_t), public :: ST_rng

  !> Parallel random number generator
  type(prng_t), public :: ST_prng

  !> Avoid dielectric relaxation time step constraint by limiting flux
  logical, public, protected :: ST_drt_limit_flux = .false.

  !> Ensure that flux limiting does not lead to fields higher than this
  real(dp), public, protected :: ST_drt_max_field = 1.0e100_dp

  !> Use source factor to prevent unphysical effects due to diffusion
  integer, public, protected :: ST_source_factor

  integer, public, parameter :: source_factor_none = 0
  integer, public, parameter :: source_factor_flux = 1
  integer, public, parameter :: source_factor_original_flux = 2

  !> Minimum number of electrons per cell to include source terms
  real(dp), public, protected :: ST_source_min_electrons_per_cell = -1e100_dp

  !> End time of the simulation
  real(dp), public, protected :: ST_end_time = 10e-9_dp

  !> Whether streamer length is used as a simulation stopping
  logical, public, protected :: ST_use_end_streamer_length = .false.

  !> Wait n steps before initializing streamer begin position
  integer, public, protected :: ST_initial_streamer_pos_steps_wait = 5

  !> Streamer length at which the simulation will stop
  real(dp), public, protected :: ST_end_streamer_length = 15e-3

  !> Abort axisymmetric simulations if there is branching
  logical, public, protected :: ST_abort_axisymmetric_if_branching = .false.

  !> The size of the boxes that we use to construct our mesh
  integer, public, protected :: ST_box_size = 8

  !> Size of the coarse grid
  integer, public, protected :: ST_coarse_grid_size(NDIM) = -1

  !> Domain length per dimension
  real(dp), public, protected :: ST_domain_len(NDIM) = 16e-3_dp

  !> Origin of domain
  real(dp), public, protected :: ST_domain_origin(NDIM) = 0.0_dp

  !> Whether the domain is periodic (per dimension)
  logical, public, protected :: ST_periodic(NDIM) = .false.

  !> Whether to limit plasma reactions to a certain region
  logical, public, protected :: ST_plasma_region_enabled = .false.

  !> Limit plasma reactions to coordinates between rmin and rmax
  real(dp), public, protected :: ST_plasma_region_rmin(NDIM) = -1e100_dp

  !> Limit plasma reactions to coordinates between rmin and rmax
  real(dp), public, protected :: ST_plasma_region_rmax(NDIM) = 1e100_dp

  !> Number of V-cycles to perform per time step
  integer, public, protected :: ST_multigrid_num_vcycles = 2

  ! Stop multigrid when residual is smaller than this factor times max(|rhs|)
  real(dp), public, protected :: ST_multigrid_max_rel_residual = 1e-4_dp

  !> Global time
  real(dp), public :: global_time = 0.0_dp

  !> Global time step
  real(dp), public :: global_dt = 0.0_dp

  !> Current sum of reaction rates per thread
  real(dp), public, allocatable :: ST_current_rates(:, :)

  !> Global sum of reaction rates
  real(dp), public, allocatable :: ST_global_rates(:)

  !> Current sum of J.E per thread
  real(dp), public, allocatable :: ST_current_JdotE(:, :)

  !> Per how many iterations the electric current is computed
  integer, public, protected :: current_update_per_steps = 1000*1000

  !> Electric current through electrodes due to J.E
  real(dp), public :: ST_global_JdotE_current

  !> Electric current through electrodes due to displacement current
  real(dp), public :: ST_global_displ_current

  !> Global sum of J.E
  real(dp), public :: ST_global_JdotE

  ! To keep track of the computational cost of different parts
  real(dp), public :: wc_time_flux = 0.0_dp
  real(dp), public :: wc_time_source = 0.0_dp
  real(dp), public :: wc_time_copy_state = 0.0_dp
  real(dp), public :: wc_time_field = 0.0_dp
  real(dp), public :: wc_time_output = 0.0_dp
  real(dp), public :: wc_time_refine = 0.0_dp
  real(dp), public :: wc_time_photoi = 0.0_dp

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
    use m_dt
    use m_model
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg  !< The configuration for the simulation
    integer, intent(in)        :: ndim !< Number of dimensions
    integer                    :: n, k, n_threads, ix_chemistry
    character(len=name_len)    :: prolong_method, bc_method
    character(len=name_len)    :: source_factor = "none"
    character(len=string_len)  :: tmp_str
    integer                    :: rng_int4_seed(4) = &
         [8123, 91234, 12399, 293434]
    integer(int64)             :: rng_int8_seed(2)
    real(dp)                   :: tmp
    integer                    :: flux_ix
    logical                    :: write_source_factor = .false.

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

    all_densities = species_itree(n_gas_species+1:n_species)

    if (model_has_energy_equation) then
       i_electron_energy = af_find_cc_variable(tree, "e_energy")
       call af_add_fc_variable(tree, "flux_energy", ix=flux_energy, &
            write_binary=.false.)
       flux_num_electron_vars = 2
    else
       flux_num_electron_vars = 1
    end if

    flux_num_species = flux_num_electron_vars+transport_data_ions%n_mobile_ions
    allocate(flux_species(flux_num_species))
    allocate(flux_species_charge(flux_num_species))
    allocate(flux_species_charge_sign(flux_num_species))
    allocate(flux_variables(flux_num_species))

    flux_species(1)             = i_electron
    flux_species_charge(1)      = -1
    flux_species_charge_sign(1) = -1
    flux_variables(1)           = flux_elec

    if (model_has_energy_equation) then
       flux_species(2) = i_electron_energy
       flux_species_charge(2) = 0
       flux_species_charge_sign(2) = -1 ! Used to determine upwind direction
       flux_variables(2) = flux_energy
    end if

    do n = 1, transport_data_ions%n_mobile_ions
       flux_ix = flux_num_electron_vars + n
       flux_species(flux_ix) = af_find_cc_variable(tree, &
            trim(transport_data_ions%names(n)))

       ! Get index in chemistry list and determine charge
       ix_chemistry = species_index(trim(transport_data_ions%names(n)))
       flux_species_charge(flux_ix) = species_charge(ix_chemistry)
       flux_species_charge_sign(flux_ix) = sign(1, species_charge(ix_chemistry))

       call af_add_fc_variable(tree, trim(transport_data_ions%names(n)), &
            ix=flux_variables(flux_ix), write_binary=.false.)
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

    ! Add one copy so that the old value can be restored
    call af_add_cc_variable(tree, "phi", ix=i_phi, n_copies=2)
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

    call CFG_add_get(cfg, "use_electrode", ST_use_electrode, &
         "Whether to include an electrode")
    if (ST_use_electrode) then
       call af_add_cc_variable(tree, "lsf", ix=i_lsf)
    end if

    bc_method = "neumann_zero"
    call CFG_add_get(cfg, "species_boundary_condition", &
         bc_method, &
         "Boundary condition for the plasma species")
    select case (bc_method)
    case ("neumann_zero")
       bc_species => af_bc_neumann_zero
    case ("dirichlet_zero")
       bc_species => bc_species_dirichlet_zero
    case default
       print *, "Unknown boundary condition: ", trim(bc_method)
       print *, "Try neumann_zero or dirichlet_zero"
       error stop
    end select

    call CFG_add_get(cfg, "compute_power_density", compute_power_density, &
         "Whether to compute the deposited power density")

    if (compute_power_density) then
       call af_add_cc_variable(tree, "power_density", ix = i_power_density)
    end if

    call CFG_add_get(cfg, "use_end_streamer_length", ST_use_end_streamer_length, &
         "Whether the length of the streamer is used to end the simulation")
    call CFG_add_get(cfg, "end_streamer_length", ST_end_streamer_length, &
         "Streamer length at which the simulation will end.")
    call CFG_add_get(cfg, "initial_streamer_pos_steps_wait", &
         ST_initial_streamer_pos_steps_wait, &
         "Number of simulation steps to wait before initializing "&
         "the starting position of the streamer")

    call CFG_add_get(cfg, "abort_axisymmetric_if_branching", &
         ST_abort_axisymmetric_if_branching, &
         "Abort axisymmetric simulations if there is branching")

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

    call CFG_add_get(cfg, "plasma_region_enabled", ST_plasma_region_enabled, &
         "Whether to limit plasma reactions to a certain region")
    call CFG_add_get(cfg, "plasma_region_rmin", ST_plasma_region_rmin, &
         "Limit plasma reactions to coordinates between rmin and rmax")
    call CFG_add_get(cfg, "plasma_region_rmax", ST_plasma_region_rmax, &
         "Limit plasma reactions to coordinates between rmin and rmax")

    if (all(ST_coarse_grid_size == -1)) then
       ! Not set, automatically determine size
       ST_coarse_grid_size = ST_box_size * &
            nint(ST_domain_len / minval(ST_domain_len))
    end if

    tmp = maxval(ST_domain_len/ST_coarse_grid_size) / &
         minval(ST_domain_len/ST_coarse_grid_size)
    if (tmp > 1.001_dp) then
       print *, "!!! Warning: using non-square grid cells"
       write(*, "(A,F12.4)") " !!! Maximal aspect ratio:", tmp
    end if

    call CFG_add_get(cfg, "multigrid_num_vcycles", ST_multigrid_num_vcycles, &
         "Number of V-cycles to perform per time step")
    call CFG_add_get(cfg, "multigrid_max_rel_residual", &
         ST_multigrid_max_rel_residual, &
         "Stop multigrid when residual is smaller than this factor times max(|rhs|)")

    call CFG_add_get(cfg, "current_update_per_steps", &
         current_update_per_steps, &
         "Per how many iterations the electric current is computed")

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

    call CFG_add_get(cfg, "fixes%drt_max_field", ST_drt_max_field, &
         "Enable flux limiting, but prevent field from exceeding this value")
    if (ST_drt_max_field < 1e100_dp) then
       error stop "fixes%drt_max_field not yet implemented"
       ST_drt_limit_flux = .true.
    end if

    call CFG_add_get(cfg, "fixes%source_factor", source_factor, &
         "Use source factor to prevent unphysical effects due to diffusion")
    call CFG_add_get(cfg, "fixes%write_source_factor", write_source_factor, &
         "Whether to write the source factor to the output")
    call CFG_add_get(cfg, "fixes%source_min_electrons_per_cell", &
         ST_source_min_electrons_per_cell, &
         "Minimum number of electrons per cell to include source terms")

    select case (source_factor)
    case ("none")
       ST_source_factor = source_factor_none
    case ("flux")
       if (model_has_energy_equation) &
            error stop "source_factor flux incompatible with energy eq."
       ST_source_factor = source_factor_flux
    case default
       print *, "Options fixes%source_factor: none, flux"
       error stop "Unknown fixes%source_factor"
    end select

    if (ST_source_factor > 0 .and. write_source_factor) then
       call af_add_cc_variable(tree, "srcfac", ix=i_srcfac)
    end if

    call CFG_add_get(cfg, "rng_seed", rng_int4_seed, &
         "Seed for random numbers; if all zero, generate randomly")

    if (all(rng_int4_seed == 0)) then
       rng_int4_seed = get_random_seed()
       print *, "RNG seed: ", rng_int4_seed

       ! Store the updated seed in the configuration
       write(tmp_str, *) "rng_seed = ", rng_int4_seed
       call CFG_update_from_line(cfg, tmp_str)
    end if

    rng_int8_seed = transfer(rng_int4_seed, rng_int8_seed)
    call ST_rng%set_seed(rng_int8_seed)
    n_threads = af_get_max_threads()
    call ST_prng%init_parallel(n_threads, ST_rng)

    ! Initialize global storage of reaction rates, +32 to avoid threads writing
    ! to nearby memory
    allocate(ST_current_rates(n_reactions+32, n_threads))
    allocate(ST_global_rates(n_reactions))
    allocate(ST_current_JdotE(1+32, n_threads))
    ST_global_rates = 0
    ST_global_JdotE = 0

  end subroutine ST_initialize

  !> Get a random seed based
  function get_random_seed() result(seed)
    use iso_fortran_env, only: int64
    integer  :: seed(4)
    integer  :: i
    real(dp) :: rr
    integer(int64) :: time

    ! Set a random seed (this does not work on all systems)
    call random_seed()

    ! Get some count of the time
    call system_clock(time)

    do i = 1, 4
       call random_number(rr)
       seed(i) = ieor(int(time), int(huge(1) * rr))
    end do
  end function get_random_seed

  !> Impose a Dirichlet zero boundary condition for plasma species in the last
  !> dimension, which is supposed to have the electrodes. We use Neumann
  !> conditions in the other dimensions. Note that this version avoids
  !> extrapolation (in contrast to the regular Dirichlet b.c.), which is more
  !> suitable for conserved species densities.
  subroutine bc_species_dirichlet_zero(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    if (af_neighb_dim(nb) == NDIM) then
       bc_type = af_bc_dirichlet_copy
       bc_val  = 0.0_dp
    else
       bc_type = af_bc_neumann
       bc_val  = 0.0_dp
    end if
  end subroutine bc_species_dirichlet_zero

end module m_streamer
