!> Template for user code, this one simply uses the default routines
module m_user
  use m_af_all
  use m_config
  use m_user_methods
  use m_table_data
  use m_lookup_table
  use m_streamer

  implicit none
  private

  ! Public methods
  public :: user_initialize
  type(LT_t), public, protected :: potential_table
  integer                       :: upper_electrode = 1
  integer                       :: lower_electrode = 2

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    user_initial_conditions => my_init_cond
    user_potential_bc => potential_bc
    call potential_from_table();

  end subroutine user_initialize

  subroutine my_init_cond(box)
    type(box_t), intent(inout) :: box

    ! print *, box%ix
  end subroutine my_init_cond

  subroutine potential_from_table()
    use m_table_data
    use m_lookup_table
    character(len=string_len)  :: td_file = undefined_str
    real(dp), allocatable      :: x_data(:), y_data(:)

    ! Create a lookup table for the applied potential
    potential_table = LT_create(0.0_dp, 0.16_dp, 1000, 2)
    td_file = "applied_voltage_upper.txt";
    call table_from_file(td_file, "location[m]_vs_potential[V]", x_data, y_data)
    call LT_set_col(potential_table, upper_electrode, x_data, y_data)
    td_file = "applied_voltage_lower.txt";
    call table_from_file(td_file, "location[m]_vs_potential[V]", x_data, y_data)
    call LT_set_col(potential_table, lower_electrode, x_data, y_data)

  end subroutine potential_from_table

  !> Dirichlet boundary conditions for the potential in the last dimension,
  !> Neumann zero boundary conditions in the other directions
  subroutine potential_bc(box, nb, iv, coords, bc_val, bc_type)
    use m_field, only: current_voltage
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type
    integer                 :: n

    if (af_neighb_dim(nb) == NDIM) then
       if (af_neighb_low(nb)) then
          bc_type = af_bc_dirichlet
          do n = 1, box%n_cell**(NDIM-1)
             bc_val(n) = current_voltage * &
                  LT_get_col(potential_table, lower_electrode, coords(1, n));
          end do
       else
          bc_type = af_bc_dirichlet
          do n = 1, box%n_cell**(NDIM-1)
             bc_val(n) = current_voltage * &
                  LT_get_col(potential_table, upper_electrode, coords(1, n));
          end do
       end if
    else
       bc_type = af_bc_neumann
       bc_val = 0.0_dp
    end if
  end subroutine potential_bc

end module m_user
