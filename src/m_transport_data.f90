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

  public :: transport_data_initialize

contains

  !> Initialize the transport coefficients
  subroutine transport_data_initialize(cfg)
    use m_config
    use m_table_data

    type(CFG_t), intent(inout) :: cfg
    character(len=string_len)  :: td_file = "td_input_file.txt"
    real(dp), allocatable      :: x_data(:), y_data(:)

    ! Create a lookup table for the model coefficients
    td_tbl = LT_create(table_min_townsend, table_max_townsend, &
         table_size, td_num_var)

    call CFG_add_get(cfg, "transport_data%file", td_file, &
         "Input file with transport data")

    ! Fill table with data
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
  end subroutine transport_data_initialize

end module m_transport_data
