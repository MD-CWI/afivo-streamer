!> Module that provides routines for reading in arbritrary transport data
module m_transport_data
  use m_lookup_table
  use m_types

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! ** Indices of transport data **
  integer, parameter, public :: td_mobility  = 1 !< Electron mobility
  integer, parameter, public :: td_diffusion = 2 !< Electron diffusion constant
  integer, parameter, public :: td_alpha     = 3 !< Ionization coefficient
  integer, parameter, public :: td_eta       = 4 !< Attachment coefficient
  !> Electron energy in eV (used with chemistry)
  integer, protected, public :: td_energy_eV = -1

  ! Table with transport data vs electric field
  type(LT_t), public, protected :: td_tbl

  !> Whether old style transport data is used (alpha, eta, mu, D vs V/m)
  logical, public, protected :: td_old_style = .false.

  ! TODO: move this to separate ion module
  type ion_transport_t
     integer                        :: n_mobile_ions ! Number of mobile ions
     real(dp), allocatable          :: mobilities(:) ! Mobility of the ions
     character(len=10), allocatable :: names(:)      ! Names of the ions
  end type ion_transport_t

  type(ion_transport_t), public :: transport_data_ions

  !> Secondary electron emission yield for positive ions
  real(dp), public, protected   :: ion_se_yield = 0.0_dp

  public :: transport_data_initialize

contains

  !> Initialize the transport coefficients
  subroutine transport_data_initialize(cfg)
    use m_config
    use m_table_data
    use m_gas

    type(CFG_t), intent(inout) :: cfg
    character(len=string_len)  :: td_file = undefined_str
    real(dp), allocatable      :: x_data(:), y_data(:)
    real(dp)                   :: dummy_real(0)
    character(len=10)          :: dummy_string(0)
    integer                    :: n

    call CFG_add_get(cfg, "input_data%file", td_file, &
         "Input file with transport (and reaction) data")
    if (td_file == undefined_str) error stop "input_data%file undefined"

    call CFG_add_get(cfg, "input_data%old_style", td_old_style, &
         "Use old style transport data (alpha, eta, mu, D vs V/m)")

    ! Fill table with data
    if (td_old_style) then
       if (.not. gas_constant_density) then
          error stop "Old style transport used with varying gas density"
       end if

       ! Create a lookup table for the model coefficients
       td_tbl = LT_create(table_min_townsend, table_max_townsend, table_size, 4)

       call table_from_file(td_file, "efield[V/m]_vs_mu[m2/Vs]", x_data, y_data)
       x_data = x_data * SI_to_Townsend / gas_number_density
       y_data = y_data * gas_number_density
       call LT_set_col(td_tbl, td_mobility, x_data, y_data)

       call table_from_file(td_file, "efield[V/m]_vs_dif[m2/s]", &
            x_data, y_data)
       x_data = x_data * SI_to_Townsend / gas_number_density
       y_data = y_data * gas_number_density
       call LT_set_col(td_tbl, td_diffusion, x_data, y_data)

       call table_from_file(td_file, "efield[V/m]_vs_alpha[1/m]", &
            x_data, y_data)
       x_data = x_data * SI_to_Townsend / gas_number_density
       y_data = y_data / gas_number_density
       call LT_set_col(td_tbl, td_alpha, x_data, y_data)

       call table_from_file(td_file, "efield[V/m]_vs_eta[1/m]", &
            x_data, y_data)
       x_data = x_data * SI_to_Townsend / gas_number_density
       y_data = y_data / gas_number_density
       call LT_set_col(td_tbl, td_eta, x_data, y_data)
    else
       ! Create a lookup table for the model coefficients
       td_tbl = LT_create(table_min_townsend, table_max_townsend, table_size, 5)

       call table_from_file(td_file, "Mobility *N (1/m/V/s)", x_data, y_data)
       call LT_set_col(td_tbl, td_mobility, x_data, y_data)

       call table_from_file(td_file, "Diffusion coefficient *N (1/m/s)", &
            x_data, y_data)
       call LT_set_col(td_tbl, td_diffusion, x_data, y_data)

       call table_from_file(td_file, "Townsend ioniz. coef. alpha/N (m2)", &
            x_data, y_data)
       call LT_set_col(td_tbl, td_alpha, x_data, y_data)

       call table_from_file(td_file, "Townsend attach. coef. eta/N (m2)", &
            x_data, y_data)
       call LT_set_col(td_tbl, td_eta, x_data, y_data)

       td_energy_eV = 5
       call table_from_file(td_file, "Mean energy (eV)", &
            x_data, y_data)
       call LT_set_col(td_tbl, td_energy_eV, x_data, y_data)
    end if

    call CFG_add(cfg, "input_data%mobile_ions", dummy_string, &
         "List of ions that are considered mobile", .true.)
    call CFG_add(cfg, "input_data%ion_mobilities", dummy_real, &
         "List of ion mobilities (m^2/Vs)", .true.)

    call CFG_get_size(cfg, "input_data%mobile_ions", n)

    transport_data_ions%n_mobile_ions = n
    allocate(transport_data_ions%names(n))
    allocate(transport_data_ions%mobilities(n))

    call CFG_get(cfg, "input_data%mobile_ions", transport_data_ions%names)
    call CFG_get(cfg, "input_data%ion_mobilities", &
         transport_data_ions%mobilities)

    call CFG_add_get(cfg, "input_data%ion_se_yield", ion_se_yield, &
         "Secondary electron emission yield for positive ions")

  end subroutine transport_data_initialize

end module m_transport_data
