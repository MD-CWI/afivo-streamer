!> Module that provides routines for reading in arbritrary transport data
module m_transport_data
  use m_lookup_table
  use m_types

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! ** Indices of transport data **
  integer, parameter, public :: td_num_var   = 4 ! Number of transport coefficients
  integer, parameter, public :: td_mobility  = 1 ! Electron mobility
  integer, parameter, public :: td_diffusion = 2 ! Electron diffusion constant
  integer, parameter, public :: td_alpha     = 3 ! Ionization coefficient
  integer, parameter, public :: td_eta       = 4 ! Attachment coefficient

  ! Table with transport data vs electric field
  type(LT_t), public, protected :: td_tbl

  !> Whether old style transport data is used (alpha, eta, mu, D vs V/m)
  logical, public, protected :: td_old_style = .false.

  public :: transport_data_initialize

contains

  !> Initialize the transport coefficients
  subroutine transport_data_initialize(cfg)
    use m_config
    use m_table_data
    use m_gas

    type(CFG_t), intent(inout) :: cfg
    character(len=string_len)  :: td_file = "td_input_file.txt"
    real(dp), allocatable      :: x_data(:), y_data(:)

    ! Create a lookup table for the model coefficients
    td_tbl = LT_create(table_min_townsend, table_max_townsend, &
         table_size, td_num_var)

    call CFG_add_get(cfg, "transport_data%file", td_file, &
         "Input file with transport data")
    call CFG_add_get(cfg, "transport_data%old_style", td_old_style, &
         "Use old style transport data (alpha, eta, mu, D vs V/m)")

    ! Fill table with data
    if (td_old_style) then
       if (.not. gas_constant_density) then
          error stop "Old style transport used with varying gas density"
       end if

       call table_from_file(td_file, "efield[V/m]_vs_mu[m2/Vs]", x_data, y_data)
       x_data = x_data * SI_to_Townsend
       y_data = y_data * gas_number_density
       call LT_set_col(td_tbl, td_mobility, x_data, y_data)

       call table_from_file(td_file, "efield[V/m]_vs_dif[m2/s]", &
            x_data, y_data)
       x_data = x_data * SI_to_Townsend
       y_data = y_data * gas_number_density
       call LT_set_col(td_tbl, td_diffusion, x_data, y_data)

       call table_from_file(td_file, "efield[V/m]_vs_alpha[1/m]", &
            x_data, y_data)
       x_data = x_data * SI_to_Townsend
       y_data = y_data / gas_number_density
       call LT_set_col(td_tbl, td_alpha, x_data, y_data)

       call table_from_file(td_file, "efield[V/m]_vs_eta[1/m]", &
            x_data, y_data)
       x_data = x_data * SI_to_Townsend
       y_data = y_data / gas_number_density
       call LT_set_col(td_tbl, td_eta, x_data, y_data)
    else
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
    end if
  end subroutine transport_data_initialize

end module m_transport_data
