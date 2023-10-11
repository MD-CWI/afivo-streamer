!> Module that provides routines for reading in arbritrary transport data
module m_transport_data
  use m_lookup_table
  use m_types
  use m_spline_interp

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

  ! Table with transport data vs electron energy
  type(LT_t), public, protected :: td_ee_tbl

  !> Electron mobility as a function of energy
  integer, protected, public :: td_ee_mobility = 1

  !> Electron diffusion coefficient as a function of energy
  integer, protected, public :: td_ee_diffusion = 2

  !> Electron energy loss
  integer, protected, public :: td_ee_loss = 3

  !> Field as a function of energy
  integer, protected, public :: td_ee_field = 4

  !> Whether old style transport data is used (alpha, eta, mu, D vs V/m)
  logical, public, protected :: td_old_style = .false.

  !> Maximal energy (eV) in input data (automatically updated)
  real(dp), public, protected :: td_max_eV = 20.0_dp

  ! @todo move this to separate ion module
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
    use m_units_constants
    use m_model
    type(CFG_t), intent(inout) :: cfg
    character(len=string_len)  :: td_file = undefined_str
    real(dp), allocatable      :: xx(:), yy(:)
    real(dp), allocatable      :: energy_eV(:), field_Td(:)
    real(dp)                   :: dummy_real(0), max_Td, max_eV
    character(len=10)          :: dummy_string(0)
    integer                    :: n

    call CFG_add_get(cfg, "input_data%file", td_file, &
         "Input file with transport (and reaction) data")
    if (td_file == undefined_str) error stop "input_data%file undefined"

    call CFG_add_get(cfg, "input_data%old_style", td_old_style, &
         "Use old style transport data (alpha, eta, mu, D vs V/m)")

    ! Fill table with data
    if (td_old_style) then
       if (.not. gas_constant_density) &
            error stop "Old style transport used with varying gas density"
       if (model_has_energy_equation) &
            error stop "Old style transport used with energy equation"

       call table_from_file(td_file, "efield[V/m]_vs_mu[m2/Vs]", xx, yy)
       xx = xx * SI_to_Townsend / gas_number_density
       yy = yy * gas_number_density

       ! Create a lookup table for the model coefficients
       if (table_max_townsend < 0) then
          max_Td = xx(size(xx))
       else
          max_Td = table_max_townsend
       end if

       td_tbl = LT_create(table_min_townsend, max_Td, table_size, &
            4, table_xspacing)

       call table_set_column(td_tbl, td_mobility, xx, yy)

       call table_from_file(td_file, "efield[V/m]_vs_dif[m2/s]", xx, yy)
       xx = xx * SI_to_Townsend / gas_number_density
       yy = yy * gas_number_density
       call table_set_column(td_tbl, td_diffusion, xx, yy)

       call table_from_file(td_file, "efield[V/m]_vs_alpha[1/m]", &
            xx, yy)
       xx = xx * SI_to_Townsend / gas_number_density
       yy = yy / gas_number_density
       call table_set_column(td_tbl, td_alpha, xx, yy)

       call table_from_file(td_file, "efield[V/m]_vs_eta[1/m]", &
            xx, yy)
       xx = xx * SI_to_Townsend / gas_number_density
       yy = yy / gas_number_density
       call table_set_column(td_tbl, td_eta, xx, yy)
    else
       call table_from_file(td_file, "Mobility *N (1/m/V/s)", xx, yy)

       ! Create a lookup table for the model coefficients
       if (table_max_townsend < 0) then
          max_Td = xx(size(xx))
       else
          max_Td = table_max_townsend
       end if

       td_tbl = LT_create(table_min_townsend, max_Td, &
            table_size, 5, table_xspacing)

       call table_set_column(td_tbl, td_mobility, xx, yy)

       call table_from_file(td_file, "Diffusion coefficient *N (1/m/s)", &
            xx, yy)
       call table_set_column(td_tbl, td_diffusion, xx, yy)

       call table_from_file(td_file, "Townsend ioniz. coef. alpha/N (m2)", &
            xx, yy)
       call table_set_column(td_tbl, td_alpha, xx, yy)

       call table_from_file(td_file, "Townsend attach. coef. eta/N (m2)", &
            xx, yy)
       call table_set_column(td_tbl, td_eta, xx, yy)

       td_energy_eV = 5
       call table_from_file(td_file, "Mean energy (eV)", &
            xx, yy)
       call table_set_column(td_tbl, td_energy_eV, xx, yy)
       td_max_eV = yy(size(yy))
    end if

    if (model_has_energy_equation) then
       call table_from_file(td_file, "Mean energy (eV)", field_Td, energy_eV)
       max_eV = energy_eV(size(energy_eV))
       td_ee_tbl = LT_create(0.0_dp, max_eV, table_size, 4, table_xspacing)

       call table_from_file(td_file, "Mobility *N (1/m/V/s)", xx, yy)
       if (.not. same_data(xx, field_Td)) &
            error stop "Same reduced field table required in all input data"

       ! Mobility as a function of energy
       call table_set_column(td_ee_tbl, td_ee_mobility, energy_eV, yy)

       ! Energy loss is mu E^2 as a function of energy. Prepend a zero, since at
       ! zero energy there can be no energy loss.
       yy = yy * xx**2 * Townsend_to_SI**2 * gas_number_density
       call table_set_column(td_ee_tbl, td_ee_loss, &
            [0.0_dp, energy_eV], [0.0_dp, yy])

       call table_from_file(td_file, "Diffusion coefficient *N (1/m/s)", xx, yy)
       if (.not. same_data(xx, field_Td)) &
            error stop "Same reduced field table required in all input data"

       ! Also prepend a zero, since at zero energy there can be no diffusion
       call table_set_column(td_ee_tbl, td_ee_diffusion, &
            [0.0_dp, energy_eV], [0.0_dp, yy])

       call table_set_column(td_ee_tbl, td_ee_field, &
            [0.0_dp, energy_eV], [0.0_dp, xx])
    end if

    call CFG_add(cfg, "input_data%mobile_ions", dummy_string, &
         "List of ions that are considered mobile", .true.)
    call CFG_add(cfg, "input_data%ion_mobilities", dummy_real, &
         "List of ion mobilities (m^2/Vs) at 1 bar, 300 K", .true.)

    call CFG_get_size(cfg, "input_data%mobile_ions", n)

    transport_data_ions%n_mobile_ions = n
    allocate(transport_data_ions%names(n))
    allocate(transport_data_ions%mobilities(n))

    call CFG_get(cfg, "input_data%mobile_ions", transport_data_ions%names)
    call CFG_get(cfg, "input_data%ion_mobilities", &
         transport_data_ions%mobilities)

    if (any(transport_data_ions%mobilities < 0)) &
         error stop "Ion mobilities should be given as positive numbers"

    ! Scale ion mobilities with gas number density at 300 K and 1 bar
    transport_data_ions%mobilities = transport_data_ions%mobilities * &
         (1e5_dp / (UC_boltzmann_const * 300))

    call CFG_add_get(cfg, "input_data%ion_se_yield", ion_se_yield, &
         "Secondary electron emission yield for positive ions")

  end subroutine transport_data_initialize

  !> Check whether data is the same
  pure logical function same_data(x1, x2)
    real(dp), intent(in) :: x1(:), x2(:)

    if (size(x1) == size(x2)) then
       same_data = minval(abs(x1-x2)) < tiny_real
    else
       same_data = .false.
    end if
  end function same_data

end module m_transport_data
