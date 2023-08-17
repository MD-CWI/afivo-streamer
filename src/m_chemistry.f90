#include "../afivo/src/cpp_macros.h"
!> Module for handling chemical reactions
module m_chemistry
  use m_types
  use m_af_all
  use m_lookup_table
  use m_table_data

  implicit none
  private

  !> Identifier for ionization reactions
  integer, parameter, public :: ionization_reaction = 1
  !> Identifier for attachment reactions
  integer, parameter, public :: attachment_reaction = 2
  !> Identifier for recombination reactions
  integer, parameter, public :: recombination_reaction = 3
  !> Identifier for detachment reactions
  integer, parameter, public :: detachment_reaction = 4
  !> Identifier for general reactions (not of any particular type)
  integer, parameter, public :: general_reaction = 5

  character(len=20), parameter, public :: reaction_names(*) = &
       [character(len=20) :: "ionization", "attachment", "recombination", &
       "detachment", "general"]

  !> Maximum number of coefficients for a reaction rate function
  integer, parameter :: rate_max_num_coeff = 9

  !> Basic chemical reaction type
  type reaction_t
     integer, allocatable  :: ix_in(:)            !< Index of input species
     integer, allocatable  :: ix_out(:)           !< Index of output species
     integer, allocatable  :: multiplicity_out(:) !< Multiplicity of output
     integer               :: n_species_in        !< Number of input species
     integer               :: rate_type           !< Type of reaction rate
     !> Type of reaction (e.g. ionization)
     integer               :: reaction_type = general_reaction
     real(dp)              :: rate_factor         !< Multiply rate by this factor
     integer               :: n_coeff             !< Number of stored coefficients
     !> Data for the reaction rate
     real(dp)              :: rate_data(rate_max_num_coeff) = -huge(1.0_dp)
     integer               :: lookup_table_index  !< Index in lookup table
     real(dp), allocatable :: x_data(:)
     real(dp), allocatable :: y_data(:)
     character(len=50)     :: description
  end type reaction_t

  !> Compact reaction type, for efficiency
  type tiny_react_t
     integer, allocatable  :: ix_in(:)
     integer, allocatable  :: ix_out(:)
     integer, allocatable  :: multiplicity_out(:)
  end type tiny_react_t

  !> Reaction with a field-dependent reaction rate
  integer, parameter :: rate_tabulated_field = 1

  !> Reaction of the form c1
  integer, parameter :: rate_analytic_constant = 2

  !> Reaction of the form c1 * (Td - c2)
  integer, parameter :: rate_analytic_linear = 3

  !> Reaction of the form c1 * exp(-(c2/(c3 + Td))**2)
  integer, parameter :: rate_analytic_exp_v1 = 4

  !> Reaction of the form c1 * exp(-(Td/c2)**2)
  integer, parameter :: rate_analytic_exp_v2 = 5

  !> Reaction of the form c1 * (300 / Te)**c2
  integer, parameter :: rate_analytic_k1 = 6

  !> Reaction of the form (c1 * (kB_eV * Te + c2)**2 - c3) * c4
  integer, parameter :: rate_analytic_k3 = 8

  !> Reaction of the form c1 * (Tg / 300)**c2 * exp(-c3 / Tg)
  integer, parameter :: rate_analytic_k4 = 9

  !> Reaction of the form c1 * exp(-c2 / Tg)
  integer, parameter :: rate_analytic_k5 = 10

  !> Reaction of the form c1 * Tg**c2
  integer, parameter :: rate_analytic_k6 = 11

  !> Reaction of the form c1 * (Tg / c2)**c3
  integer, parameter :: rate_analytic_k7 = 12

  !> Reaction of the form c1 * (300 / Tg)**c2
  integer, parameter :: rate_analytic_k8 = 13

  !> Reaction of the form c1 * exp(-c2 * Tg)
  integer, parameter :: rate_analytic_k9 = 14

  !> Reaction of the form 10**(c1 + c2 * (Tg - 300))
  integer, parameter :: rate_analytic_k10 = 15

  !> Reaction of the form c1 * (300 / Tg)**c2 * exp(-c3 / Tg)
  integer, parameter :: rate_analytic_k11 = 16

  !> Reaction of the form c1 * Tg**c2 * exp(-c3 / Tg)
  integer, parameter :: rate_analytic_k12 = 17

  !> Reaction of the form c1 * exp(-(c2 / (c3 + Td))**c4)
  integer, parameter :: rate_analytic_k13 = 18

  !> Reaction of the form c1 * exp(-(Td / c2)**c3)
  integer, parameter :: rate_analytic_k14 = 19

  !> Reaction of the form c1 * exp(-(c2 /(kb * (Tg + Td/c3)))**c4)
  integer, parameter :: rate_analytic_k15 = 20

  !> Maximum number of species
  integer, parameter :: max_num_species      = 100

  !> Maximum number of reactions
  integer, parameter :: max_num_reactions    = 500

  !> Number of species present
  integer, public, protected :: n_species = 0

  !> Number of gas species present
  integer, public, protected :: n_gas_species = 0

  !> Number of plasma species present
  integer, public, protected :: n_plasma_species = 0

  !> Number of reactions present
  integer, public, protected :: n_reactions = 0

  !> List of the species
  character(len=comp_len), public, protected :: species_list(max_num_species)

  !> Charge of the species
  integer, public, protected                 :: species_charge(max_num_species) = 0

  !> species_itree(n) holds the index of species n in the tree (cell-centered variables)
  integer, public, protected                 :: species_itree(max_num_species)

  !> List of reactions
  type(reaction_t), public, protected        :: reactions(max_num_reactions)

  !> A copy of the list of reactions for performance reasons
  type(tiny_react_t)                         :: tiny_react(max_num_reactions)

  !> Lookup table with reaction rates
  type(LT_t)                                 :: chemtbl

  !> List with indices of charged species
  integer, allocatable, protected :: charged_species_itree(:)

  !> List with charges of charged species
  integer, allocatable, protected :: charged_species_charge(:)

  public :: charged_species_itree
  public :: charged_species_charge

  public :: chemistry_initialize
  public :: chemistry_write_summary
  public :: chemistry_get_breakdown_field
  public :: get_rates
  public :: get_derivatives
  public :: species_index

  public :: read_reactions

contains

  !> Initialize module and load chemical reactions
  subroutine chemistry_initialize(tree, cfg)
    use m_config
    use m_units_constants
    use m_table_data
    use m_transport_data
    use m_gas
    use m_dt
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg
    integer                    :: n, i, i_elec
    character(len=string_len)  :: reaction_file
    character(len=comp_len)    :: tmp_name
    logical                    :: read_success

    call CFG_get(cfg, "input_data%file", reaction_file)

    if (.not. gas_constant_density) then
       ! Make sure the gas components are the first species
       n_gas_species                 = size(gas_components)
       n_species                     = n_gas_species
       species_list(1:n_gas_species) = gas_components
    end if

    call read_reactions(trim(reaction_file), read_success)

    if (.not. read_success) then
       print *, "m_chemistry: no reaction table found, using standard model"

       species_list(1) = "e"
       species_list(2) = "M+"
       species_list(3) = "M-"
       n_species       = 3
       n_reactions     = 2

       ! Ionization reaction
       if (gas_constant_density) then
          reactions(1)%ix_in = [1]
          reactions(1)%ix_out = [1, 2]
          reactions(1)%multiplicity_out = [2, 1]
          reactions(1)%n_species_in = 2
          reactions(1)%rate_type = rate_tabulated_field
          reactions(1)%rate_factor = 1.0_dp
          reactions(1)%x_data = td_tbl%x
          reactions(1)%y_data = td_tbl%rows_cols(:, td_alpha) * &
               td_tbl%rows_cols(:, td_mobility) * reactions(1)%x_data * &
               Townsend_to_SI * gas_number_density
          reactions(1)%description = "e + M > e + e + M+"

          ! Attachment reaction
          reactions(2)%ix_in = [1]
          reactions(2)%ix_out = [3]
          reactions(2)%multiplicity_out = [1]
          reactions(2)%n_species_in = 2
          reactions(2)%rate_type = rate_tabulated_field
          reactions(2)%rate_factor = 1.0_dp
          reactions(2)%x_data = td_tbl%x
          reactions(2)%y_data = td_tbl%rows_cols(:, td_eta) * &
               td_tbl%rows_cols(:, td_mobility) * reactions(2)%x_data * &
               Townsend_to_SI * gas_number_density
          reactions(2)%description = "e + M > M-"
       else
          error stop "Varying gas density not yet supported"
       end if
    end if

    ! Convert names to simple ascii
    do n = 1, n_species
       tmp_name = species_list(n)
       call to_simple_ascii(trim(tmp_name), species_list(n), &
            species_charge(n))
    end do

    ! Store reactions of the tabulated field type
    i = count(reactions(1:n_reactions)%rate_type == rate_tabulated_field)
    chemtbl = LT_create(table_min_townsend, table_max_townsend, &
         table_size, i, table_xspacing)

    i = 0
    do n = 1, n_reactions
       if (reactions(n)%rate_type == rate_tabulated_field) then
          i = i + 1
          reactions(n)%lookup_table_index = i
          if (td_bulk_scale_reactions) then
             call table_set_column(chemtbl, i, reactions(n)%x_data, &
                  reactions(n)%y_data * &
                  LT_get_col(td_tbl, td_bulk_scaling, reactions(n)%x_data))
          else
             call table_set_column(chemtbl, i, reactions(n)%x_data, &
                  reactions(n)%y_data)
          end if
       end if
    end do

    ! Also store in more memory-efficient structure
    do n = 1, n_reactions
       tiny_react(n)%ix_in            = reactions(n)%ix_in
       tiny_react(n)%ix_out           = reactions(n)%ix_out
       tiny_react(n)%multiplicity_out = reactions(n)%multiplicity_out
    end do

    ! Gas species are not stored in the tree for now
    species_itree(1:n_gas_species) = -1
    n_plasma_species = n_species - n_gas_species

    do n = n_gas_species+1, n_species
       call af_add_cc_variable(tree, trim(species_list(n)), &
            n_copies=af_advance_num_steps(time_integrator)+1, &
            ix=species_itree(n))
    end do

    ! Store list with only charged species
    n = count(species_charge(1:n_species) /= 0)
    allocate(charged_species_itree(n))
    allocate(charged_species_charge(n))

    i = 0
    do n = 1, n_species
       if (species_charge(n) /= 0) then
          i = i + 1
          charged_species_itree(i) = species_itree(n)
          charged_species_charge(i) = species_charge(n)
       end if
    end do

    call check_charge_conservation()

    ! Specify reaction types
    i_elec = species_index("e")

    do n = 1, n_reactions
       if (any(reactions(n)%ix_in == i_elec) .and. &
            .not. any(reactions(n)%ix_out == i_elec) .and. &
            .not. any(species_charge(reactions(n)%ix_in) > 0)) then
          ! In: an electron and no positive ions, out: no electrons
          reactions(n)%reaction_type = attachment_reaction
       else if (any(reactions(n)%ix_in == i_elec) .and. &
            any(reactions(n)%ix_out == i_elec .and. &
            reactions(n)%multiplicity_out == 2)) then
          ! An electron in and out (with multiplicity 2)
          reactions(n)%reaction_type = ionization_reaction
       else if (any(species_charge(reactions(n)%ix_in) /= 0) .and. &
            .not. any(species_charge(reactions(n)%ix_out) /= 0)) then
          ! In: charged species, out: no charged species
          reactions(n)%reaction_type = recombination_reaction
       else if (all(reactions(n)%ix_in /= i_elec) .and. &
            any(reactions(n)%ix_out == i_elec)) then
          ! In: no electrons, out: an electron
          reactions(n)%reaction_type = detachment_reaction
       end if
    end do

    print *, "--- List of reactions ---"
    do n = 1, n_reactions
       write(*, "(I4,' (',I0,') ',A15,A)") n, reactions(n)%n_species_in, &
            reaction_names(reactions(n)%reaction_type), &
            reactions(n)%description
    end do
    print *, "-------------------------"
    print *, ""
    print *, "--- List of gas species ---"
    do n = 1, n_gas_species
       write(*, "(I4,A)") n, ": " // species_list(n)
    end do
    print *, "-------------------------"
    print *, ""
    print *, "--- List of plasma species ---"
    do n = n_gas_species+1, n_species
       write(*, "(I4,A)") n, ": " // species_list(n)
    end do
    print *, "-------------------------"

  end subroutine chemistry_initialize

  !> Write a summary of the reactions (TODO) and the ionization and attachment
  !> coefficients (if working at constant pressure)
  subroutine chemistry_write_summary(fname)
    use m_gas
    use m_transport_data
    character(len=*), intent(in) :: fname
    real(dp), allocatable        :: fields(:)
    real(dp), allocatable        :: rates(:, :)
    real(dp), allocatable        :: eta(:), alpha(:), src(:), loss(:)
    real(dp), allocatable        :: v(:), mu(:), diff(:)
    integer                      :: n, n_fields, i_elec
    integer                      :: my_unit

    if (gas_constant_density) then
       i_elec = species_index("e")
       n_fields = td_tbl%n_points

       if (n_fields < 3) error stop "Not enough data for linear extrapolation"

       allocate(fields(n_fields))
       fields = td_tbl%x

       allocate(rates(n_fields, n_reactions))
       allocate(eta(n_fields), alpha(n_fields), src(n_fields), loss(n_fields))
       call get_rates(fields, rates, n_fields)

       loss(:)   = 0.0_dp
       src(:) = 0.0_dp

       do n = 1, n_reactions
          if (reactions(n)%reaction_type == attachment_reaction) then
             loss(:) = loss(:) + rates(:, n)
          else if (reactions(n)%reaction_type == ionization_reaction) then
             src(:) = src(:) + rates(:, n)
          end if
       end do

       allocate(diff(n_fields))
       diff = LT_get_col(td_tbl, td_diffusion, fields)

       allocate(mu(n_fields))
       allocate(v(n_fields))
       mu = LT_get_col(td_tbl, td_mobility, fields)
       v = mu * fields * Townsend_to_SI

       ! v(1) is zero, so extrapolate linearly
       eta(2:) = loss(2:) / v(2:)
       eta(1) = 2 * eta(2) - eta(3)
       alpha(2:) = src(2:) / v(2:)
       alpha(1) = 2 * alpha(2) - alpha(3)

       ! Write to a file
       open(newunit=my_unit, file=trim(fname), action="write")
       write(my_unit, "(A)") &
            "E/N[Td] E[V/m] Electron_mobility[m^2/(Vs)] &
            &Electron_diffusion[m^2/s] Townsend_ioniz._coef._alpha[1/m] &
            &Townsend_attach._coef._eta[1/m] Ionization_rate[1/s] &
            &Attachment_rate[1/s]"
       do n = 1, n_fields
          write(my_unit, *) fields(n), &
               fields(n) * Townsend_to_SI * gas_number_density, &
               mu(n) / gas_number_density, &
               diff(n) / gas_number_density, alpha(n), eta(n), &
               src(n), loss(n)
       end do
       write(my_unit, *) ""
       close(my_unit)
    end if
  end subroutine chemistry_write_summary

  subroutine check_charge_conservation()
    integer :: n, q_in, q_out

    do n = 1, n_reactions
       q_in = sum(species_charge(reactions(n)%ix_in))
       q_out = sum(species_charge(reactions(n)%ix_out) * &
            reactions(n)%multiplicity_out)
       if (q_in /= q_out) then
          print *, trim(reactions(n)%description)
          error stop "Charge is not conserved in the above reaction"
       end if
    end do
  end subroutine check_charge_conservation

  !> Get the breakdown field in Townsend
  subroutine chemistry_get_breakdown_field(field_td, min_growth_rate)
    use m_transport_data
    !> Breakdown field in Townsend
    real(dp), intent(out) :: field_td
    !> Minimal growth rate for breakdown
    real(dp), intent(in)  :: min_growth_rate

    integer               :: n, n_fields
    real(dp), allocatable :: fields(:), rates(:, :), src(:), loss(:)

    n_fields = td_tbl%n_points
    allocate(fields(n_fields))
    fields = td_tbl%x

    allocate(rates(n_fields, n_reactions))
    allocate(src(n_fields), loss(n_fields))
    call get_rates(fields, rates, n_fields)

    loss(:) = 0.0_dp
    src(:)  = 0.0_dp

    do n = 1, n_reactions
       if (reactions(n)%reaction_type == attachment_reaction) then
          loss(:) = loss(:) + rates(:, n)
       else if (reactions(n)%reaction_type == ionization_reaction) then
          src(:) = src(:) + rates(:, n)
       end if
    end do

    do n = n_fields, 1, -1
       if (src(n) - loss(n) < min_growth_rate) exit
    end do

    field_td = 0.0_dp
    if (n > 0) field_td = fields(n)
  end subroutine chemistry_get_breakdown_field

  !> Compute reaction rates
  !>
  !> @todo These reactions do not take into account a variable gas_temperature
  subroutine get_rates(fields, rates, n_cells)
    use m_units_constants
    use m_gas
    use m_transport_data
    integer, intent(in)   :: n_cells                     !< Number of cells
    real(dp), intent(in)  :: fields(n_cells)             !< The field (in Td) in the cells
    real(dp), intent(out) :: rates(n_cells, n_reactions) !< The reaction rates
    integer               :: n, n_coeff
    real(dp)              :: c0, c(rate_max_num_coeff)
    real(dp)              :: Te(n_cells)                 !> Electron Temperature in Kelvin
    logical               :: Te_available

    ! Conversion factor to go from eV to Kelvin
    real(dp), parameter   :: electron_eV_to_K = 2 * UC_elec_volt / &
         (3 * UC_boltzmann_const)

    Te_available = .false.

    do n = 1, n_reactions
       ! A factor that the reaction rate is multiplied with, for example to
       ! account for a constant gas number density and if the length unit is cm instead of m
       c0 = reactions(n)%rate_factor

       ! Coefficients for the reaction rate
       n_coeff = reactions(n)%n_coeff
       c(1:n_coeff) = reactions(n)%rate_data(1:n_coeff)

       select case (reactions(n)%rate_type)
       case (rate_tabulated_field)
          rates(:, n) = c0 * LT_get_col(chemtbl, &
               reactions(n)%lookup_table_index, fields)
       case (rate_analytic_constant)
          rates(:, n) = c0 * c(1)
       case (rate_analytic_linear)
          rates(:, n) = c0 * c(1) * (fields - c(2))
       case (rate_analytic_exp_v1)
          rates(:, n) = c0 * c(1) * exp(-(c(2) / (c(3) + fields))**2)
       case (rate_analytic_exp_v2)
          rates(:, n) = c0 * c(1) * exp(-(fields/c(2))**2)
       case (rate_analytic_k1)
          if (.not. Te_available) then
             Te = electron_eV_to_K * LT_get_col(td_tbl, td_energy_eV, fields)
          end if
          rates(:, n) = c0 * c(1) * (300 / Te)**c(2)
       case (rate_analytic_k3)
          if (.not. Te_available) then
             Te = electron_eV_to_K * LT_get_col(td_tbl, td_energy_eV, fields)
          end if
          ! We convert boltzmann_const from J / K to eV / K
          rates(:, n) = c0 * (c(1) * ((UC_boltzmann_const / UC_elec_volt) * Te + c(2))**2 - c(3)) * c(4)
       case (rate_analytic_k4)
          rates(:, n) = c0 * c(1) * (gas_temperature / 300)**c(2) * exp(-c(3) / gas_temperature)
       case (rate_analytic_k5)
          rates(:, n) = c0 * c(1) * exp(-c(2) / gas_temperature)
       case (rate_analytic_k6)
          rates(:, n) = c0 * c(1) * gas_temperature**c(2)
       case (rate_analytic_k7)
          rates(:, n) = c0 * c(1) * (gas_temperature / c(2))**c(3)
       case (rate_analytic_k8)
          rates(:, n) = c0 * c(1) * (300 / gas_temperature)**c(2)
       case (rate_analytic_k9)
          rates(:, n) = c0 * c(1) * exp(-c(2) * gas_temperature)
       case (rate_analytic_k10)
          rates(:, n) = c0 * 10**(c(1) + c(2) * (gas_temperature - 300))
       case (rate_analytic_k11)
          rates(:, n) = c0 * c(1) * (300 / gas_temperature)**c(2) * exp(-c(3) / gas_temperature)
       case (rate_analytic_k12)
          rates(:, n) = c0 * c(1) * gas_temperature**c(2) * exp(-c(3) / gas_temperature)
       case (rate_analytic_k13)
          rates(:, n) = c0 * c(1) * exp(-(c(2) / (c(3) + fields))**c(4))
       case (rate_analytic_k14)
          rates(:, n) = c0 * c(1) * exp(-(fields / c(2))**c(3))
       case (rate_analytic_k15)
          ! Note that this reaction depends on Ti, ionic temperature, which according to Galimberti(1979),
          ! Ti = T_gas + fields/c(3), with c(3) = 0.18 Td/Kelvin, UC_boltzmann_const is in J/Kelvin,
          ! c(2) is given in Joule in the input file
          rates(:, n) = c0 * c(1) * exp(-(c(2) / (UC_boltzmann_const * (gas_temperature + fields/c(3))))**c(4))
      end select
    end do
  end subroutine get_rates

  !> Compute derivatives due to chemical reactions. Note that the 'rates'
  !> argument is modified.
  subroutine get_derivatives(dens, rates, derivs, n_cells)
    integer, intent(in)     :: n_cells
    !> Species densities
    real(dp), intent(in)    :: dens(n_cells, n_species)
    !> On input, reaction rate coefficients. On output, actual reaction rates.
    real(dp), intent(inout) :: rates(n_cells, n_reactions)
    !> Derivatives of the chemical species
    real(dp), intent(out)   :: derivs(n_cells, n_species)
    integer                 :: n, i, ix

    derivs(:, :) = 0.0_dp

    ! Loop over reactions and add to derivative
    do n = 1, n_reactions
       ! Determine production rate of 'full' reaction
       rates(:, n) = rates(:, n) * &
            product(dens(:, tiny_react(n)%ix_in), dim=2)

       ! Input species are removed
       do i = 1, size(tiny_react(n)%ix_in)
          ix = tiny_react(n)%ix_in(i)
          derivs(:, ix) = derivs(:, ix) - rates(:, n)
       end do

       ! Output species are created
       do i = 1, size(tiny_react(n)%ix_out)
          ix = tiny_react(n)%ix_out(i)
          derivs(:, ix) = derivs(:, ix) + rates(:, n) * &
               tiny_react(n)%multiplicity_out(i)
       end do
    end do
  end subroutine get_derivatives

  !> Try to read a list species whose production will be ignored
  subroutine read_ignored_species(filename, ignored_species)
    character(len=*), intent(in) :: filename
    character(len=comp_len), allocatable, intent(inout) :: &
         ignored_species(:)
    character(len=comp_len)      :: tmp(max_num_species)
    character(len=string_len)    :: line
    integer                      :: my_unit, n_ignored

    n_ignored = 0

    open(newunit=my_unit, file=filename, action="read")

    ! Find list of reactions
    do
       read(my_unit, "(A)", end=998) line
       line = adjustl(line)

       if (line == "ignored_species") then
          ! Read next line starting with at least 5 dashes
          read(my_unit, "(A)") line
          if (line(1:5) /= "-----") &
               error stop "ignored_species not followed by -----"
          exit
       end if
    end do

    ! Read ignored species, one per line
    do
       read(my_unit, "(A)", end=999) line
       line = adjustl(line)

       ! Ignore comments
       if (line(1:1) == "#") cycle

       ! Exit when we read a line of dashes
       if (line(1:5) == "-----") exit

       n_ignored = n_ignored + 1
       read(line, *) tmp(n_ignored)
    end do

998 close(my_unit)
    ignored_species = tmp(1:n_ignored)
    return

    ! Error messages
999 error stop "read_ignored_species: no closing dashes"
  end subroutine read_ignored_species

  !> Read reactions from a file
  subroutine read_reactions(filename, read_success)
    character(len=*), intent(in) :: filename
    logical, intent(out)         :: read_success
    character(len=string_len)    :: line
    integer, parameter           :: field_len    = 50
    character(len=field_len)     :: data_value(max_num_reactions)
    character(len=field_len)     :: reaction(max_num_reactions)
    character(len=field_len)     :: how_to_get(max_num_reactions)
    character(len=10)            :: length_unit(max_num_reactions)
    type(reaction_t)             :: new_reaction
    integer                      :: my_unit
    integer                      :: n_reactions_found
    integer, parameter           :: n_fields_max = 40
    integer                      :: i0(n_fields_max), i1(n_fields_max)
    integer                      :: n, i, k, n_found, lo, hi
    logical                      :: keep_reaction

    character(len=comp_len), allocatable :: ignored_species(:)

    type group
       character(len=name_len)               :: name
       character(len=field_len), allocatable :: members(:)
    end type group

    integer, parameter :: max_groups = 10
    integer            :: i_group, group_size
    type(group)        :: groups(max_groups)

    call read_ignored_species(filename, ignored_species)

    open(newunit=my_unit, file=filename, action="read")

    n_reactions_found = 0
    i_group           = 0
    group_size        = 0
    read_success      = .false.

    ! Find list of reactions
    do
       read(my_unit, "(A)", end=998) line
       line = adjustl(line)

       if (line == "reaction_list") then
          ! Read next line starting with at least 5 dashes
          read(my_unit, "(A)") line
          if (line(1:5) /= "-----") &
               error stop "read_reactions: reaction_list not followed by -----"
          exit
       end if
    end do

    ! Start reading reactions
    do
       read(my_unit, "(A)", end=999) line
       line = adjustl(line)

       ! Ignore comments
       if (line(1:1) == "#") cycle

       ! Exit when we read a line of dashes
       if (line(1:5) == "-----") exit

       ! Group syntax for multiple species
       if (line(1:1) == "@") then
          call get_fields_string(line, "=,", n_fields_max, n_found, i0, i1)
          i_group = i_group + 1

          if (i_group > max_groups) error stop "Too many groups"

          if (i_group == 1) then
             group_size = n_found - 1
          else if (n_found - 1 /= group_size) then
             print *, trim(line)
             error stop "Groups for a reaction should have the same size"
          end if

          groups(i_group)%name = line(i0(1):i1(1))
          allocate(groups(i_group)%members(n_found - 1))

          do n = 2, n_found
             groups(i_group)%members(n-1) = adjustl(line(i0(n):i1(n)))
          end do
          cycle
       end if

       if (i_group > 0) then
          ! Handle groups
          lo                   = n_reactions_found
          hi                   = n_reactions_found+group_size-1
          reaction(lo+1:hi)    = reaction(lo)
          how_to_get(lo+1:hi)  = how_to_get(lo)
          data_value(lo+1:hi)  = data_value(lo)
          length_unit(lo+1:hi) = length_unit(lo)
          n_reactions_found    = hi

          do k = 1, group_size
             do i = 1, i_group
                n = lo + k - 1
                reaction(n)   = string_replace(reaction(n), &
                     groups(i)%name, groups(i)%members(k))
                how_to_get(n) = string_replace(how_to_get(n), &
                     groups(i)%name, groups(i)%members(k))
                data_value(n) = string_replace(data_value(n), &
                     groups(i)%name, groups(i)%members(k))
             end do
          end do

          ! Check if there are duplicate reactions
          do n = lo, hi
             if (count(reaction(lo:hi) == reaction(n)) > 1 .and. &
                  count(data_value(lo:hi) == data_value(n)) > 1) then
                do k = lo, hi
                   print *, trim(reaction(k)), ",", data_value(k)
                end do
                error stop "Groups lead to duplicate reactions"
             end if
          end do

          i_group    = 0
          group_size = 0
       end if

       call get_fields_string(line, ",", n_fields_max, n_found, i0, i1)

       if (n_found < 3 .or. n_found > 4) then
          print *, trim(line)
          error stop "Invalid chemistry syntax"
       end if

       if (n_reactions_found >= max_num_reactions) &
            error stop "Too many reactions, increase max_num_reactions"

       n_reactions_found             = n_reactions_found + 1
       reaction(n_reactions_found)   = line(i0(1):i1(1))
       how_to_get(n_reactions_found) = line(i0(2):i1(2))
       data_value(n_reactions_found) = line(i0(3):i1(3))

       ! Fourth entry can hold a custom length unit, default is meter
       if (n_found > 3) then
          length_unit(n_reactions_found) = line(i0(4):i1(4))
       else
          length_unit(n_reactions_found) = "m"
       end if
    end do

998 continue

    ! Close the file (so that we can re-open it for reading data)
    close(my_unit)

    n_reactions = 0

    ! Now parse the reactions
    do n = 1, n_reactions_found
       call parse_reaction(trim(reaction(n)), new_reaction, &
            ignored_species, keep_reaction)

       if (keep_reaction) then
          n_reactions = n_reactions + 1
       else
          cycle
       end if

       new_reaction%description = trim(reaction(n))

       ! IMPORTANT: If you change the reactions below, do not forget to update
       ! documentation/chemistry.md accordingly!
       select case (how_to_get(n))
       case ("field_table")
          ! Reaction data should be present in the same file
          call read_reaction_table(filename, &
               trim(data_value(n)), new_reaction)
          new_reaction%n_coeff = 0
       case ("c1")
          new_reaction%rate_type = rate_analytic_constant
          new_reaction%n_coeff = 1
          read(data_value(n), *) new_reaction%rate_data(1)
       case ("c1*(Td-c2)")
          new_reaction%rate_type = rate_analytic_linear
          new_reaction%n_coeff = 2
          read(data_value(n), *) new_reaction%rate_data(1:2)
       case ("c1*exp(-(c2/(c3+Td))**2)")
          new_reaction%rate_type = rate_analytic_exp_v1
          new_reaction%n_coeff = 3
          read(data_value(n), *) new_reaction%rate_data(1:3)
       case ("c1*exp(-(Td/c2)**2)")
          new_reaction%rate_type = rate_analytic_exp_v2
          new_reaction%n_coeff = 2
          read(data_value(n), *) new_reaction%rate_data(1:2)
       case ("c1*(300/Te)**c2")
          new_reaction%rate_type = rate_analytic_k1
          new_reaction%n_coeff = 2
          read(data_value(n), *) new_reaction%rate_data(1:2)
       case ("(c1*(kB_eV*Te+c2)**2-c3)*c4")
          new_reaction%rate_type = rate_analytic_k3
          new_reaction%n_coeff = 4
          read(data_value(n), *) new_reaction%rate_data(1:4)
       case ("c1*(Tg/300)**c2*exp(-c3/Tg)")
          new_reaction%rate_type = rate_analytic_k4
          new_reaction%n_coeff = 3
          read(data_value(n), *) new_reaction%rate_data(1:3)
       case ("c1*exp(-c2/Tg)")
          new_reaction%rate_type = rate_analytic_k5
          new_reaction%n_coeff = 2
          read(data_value(n), *) new_reaction%rate_data(1:2)
       case ("c1*Tg**c2")
          new_reaction%rate_type = rate_analytic_k6
          new_reaction%n_coeff = 2
          read(data_value(n), *) new_reaction%rate_data(1:2)
       case ("c1*(Tg/c2)**c3")
          new_reaction%rate_type = rate_analytic_k7
          new_reaction%n_coeff = 3
          read(data_value(n), *) new_reaction%rate_data(1:3)
       case ("c1*(300/Tg)**c2")
          new_reaction%rate_type = rate_analytic_k8
          new_reaction%n_coeff = 2
          read(data_value(n), *) new_reaction%rate_data(1:2)
       case ("c1*exp(-c2*Tg)")
          new_reaction%rate_type = rate_analytic_k9
          new_reaction%n_coeff = 2
          read(data_value(n), *) new_reaction%rate_data(1:2)
       case ("10**(c1+c2*(Tg-300))")
          new_reaction%rate_type = rate_analytic_k10
          new_reaction%n_coeff = 2
          read(data_value(n), *) new_reaction%rate_data(1:2)
       case ("c1*(300/Tg)**c2*exp(-c3/Tg)")
          new_reaction%rate_type = rate_analytic_k11
          new_reaction%n_coeff = 3
          read(data_value(n), *) new_reaction%rate_data(1:3)
       case ("c1*Tg**c2*exp(-c3/Tg)")
          new_reaction%rate_type = rate_analytic_k12
          new_reaction%n_coeff = 3
          read(data_value(n), *) new_reaction%rate_data(1:3)
       case ("c1*exp(-(c2/(c3+Td))**c4)")
          new_reaction%rate_type = rate_analytic_k13
          new_reaction%n_coeff = 4
          read(data_value(n), *) new_reaction%rate_data(1:4)
       case ("c1*exp(-(Td/c2)**c3)")
          new_reaction%rate_type = rate_analytic_k14
          new_reaction%n_coeff = 3
          read(data_value(n), *) new_reaction%rate_data(1:3)
       case ("c1*exp(-(c2/(kb*(Tg+Td/c3)))**c4)")
          new_reaction%rate_type = rate_analytic_k15
          new_reaction%n_coeff = 4
          read(data_value(n), *) new_reaction%rate_data(1:4)
       case default
          print *, "Unknown rate type: ", trim(how_to_get(n))
          print *, "For reaction:      ", trim(reaction(n))
          print *, "In file:           ", trim(filename)
          print *, "See documentation/chemistry.md"
          if (how_to_get(n) /= "field_table" .and. &
               index(how_to_get(n), "c1") == 0) then
             print *, "You probably use the old reaction format"
             print *, "Try to use tools/chemistry_update_reactions.sh"
             print *, "See also the chemistry documentation"
          end if
          error stop "Unknown chemical reaction"
       end select

       ! Correct for length unit in the rate function (e.g. [k] = cm3/s)
       select case (length_unit(n))
       case ("m")
          continue
       case ("cm")
          ! 1 cm^3 = 1e-6 m^3
          new_reaction%rate_factor = new_reaction%rate_factor * &
               1e-6_dp**(new_reaction%n_species_in-1)
       case default
          print *, "For reaction: ", trim(reaction(n))
          print *, "Invalid length unit: ", trim(length_unit(n))
          error stop
       end select

       reactions(n_reactions) = new_reaction
    end do

    if (n_reactions > 0) read_success = .true.
    return

    ! Error messages
999 error stop "read_reactions: no closing dashes for reaction list"
  end subroutine read_reactions

  !> Read a reaction table
  subroutine read_reaction_table(filename, dataname, rdata)
    use m_transport_data
    character(len=*), intent(in)    :: filename
    character(len=*), intent(in)    :: dataname
    type(reaction_t), intent(inout) :: rdata

    rdata%rate_type = rate_tabulated_field
    call table_from_file(filename, dataname, rdata%x_data, rdata%y_data)
  end subroutine read_reaction_table

  !> Parse a reaction and store it
  subroutine parse_reaction(reaction_text, reaction, ignored_species, &
       keep_reaction)
    use m_gas
    character(len=*), intent(in)  :: reaction_text
    type(reaction_t), intent(out) :: reaction
    character(len=comp_len), intent(in) :: ignored_species(:)
    logical, intent(out)          :: keep_reaction
    integer, parameter            :: max_components = 100
    character(len=comp_len)       :: component
    integer                       :: i, ix, n, n_found, multiplicity
    integer                       :: n_in, n_out
    integer                       :: i0(max_components)
    integer                       :: i1(max_components)
    integer                       :: ix_in(max_components)
    integer                       :: ix_out(max_components)
    integer                       :: multiplicity_out(max_components)
    logical                       :: left_side, is_gas_species
    real(dp)                      :: rfactor

    call get_fields_string(reaction_text, " ", max_components, n_found, i0, i1)

    left_side             = .true.
    keep_reaction         = .true.
    n_in                  = 0
    n_out                 = 0
    rfactor               = 1.0_dp
    reaction%n_species_in = 0

    do n = 1, n_found
       component = reaction_text(i0(n):i1(n))

       if (component == "+") cycle

       if (component == "->") then
          left_side = .false.
          cycle
       end if

       ! Assume we have a multiplicity less than 10
       if (lge(component(1:1), '1') .and. lle(component(1:1), '9')) then
          read(component(1:1), *) multiplicity
          component = component(2:)
       else
          multiplicity = 1
       end if

       if (left_side) then
          reaction%n_species_in = reaction%n_species_in + multiplicity
       end if

       ! If the gas density is constant, remove gas species from the reaction
       if (gas_constant_density) then
          ix = gas_index(component)
          if (ix /= -1) then
             ! Multiply reaction by constant density
             if (left_side) rfactor = rfactor * gas_densities(ix)
             cycle
          end if

          if (component == "M") then
             ! Assume this stands for 'any molecule'
             if (left_side) rfactor = rfactor * gas_number_density
             cycle
          end if
       end if

       ! Handle ignored species
       if (findloc(ignored_species, component, dim=1) > 0) then
          is_gas_species = (gas_index(component) > 0 .or. component == "M")

          if (left_side .and. .not. is_gas_species) then
             ! Ignore the whole reaction, since this species will not be
             ! produced (and will thus have zero density)
             keep_reaction = .false.
             return
          else
             ! Ignore the production of this species, but keep the reaction
             cycle
          end if
       end if

       ix = species_index(component)
       if (ix == -1) then
          if (n_species >= max_num_species) &
               error stop "Too many species, increase max_num_species"
          n_species        = n_species + 1
          ix               = n_species
          species_list(ix) = trim(component)
       end if

       if (left_side) then
          do i = 1, multiplicity
             n_in = n_in + 1
             ix_in(n_in) = ix
          end do
       else
          ! Check if species is already present in right-hand side
          do i = 1, n_out
             if (ix == ix_out(i)) exit
          end do

          if (i <= n_out) then
             multiplicity_out(i) = multiplicity_out(i) + multiplicity
          else
             ! If not already present, add the species
             n_out = n_out + 1
             ix_out(n_out) = ix
             multiplicity_out(n_out) = multiplicity
          end if
       end if
    end do

    reaction%ix_in            = ix_in(1:n_in)
    reaction%ix_out           = ix_out(1:n_out)
    reaction%multiplicity_out = multiplicity_out(1:n_out)
    reaction%rate_factor      = rfactor

    if (n_in == 0) then
       print *, "Error in the following reaction:"
       print *, trim(reaction_text)
       error stop "No input species"
    end if
  end subroutine parse_reaction

  !> Find index of a species, return -1 if not found
  elemental integer function species_index(name)
    character(len=*), intent(in) :: name
    do species_index = 1, n_species
       if (species_list(species_index) == name) exit
    end do
    if (species_index == n_species+1) species_index = -1
  end function species_index

  !> Routine to find the indices of entries in a string
  subroutine get_fields_string(line, delims, n_max, n_found, ixs_start, ixs_end)
    !> The line from which we want to read
    character(len=*), intent(in)  :: line
    !> A string with delimiters. For example delims = " ,'"""//char(9)
    character(len=*), intent(in)  :: delims
    !> Maximum number of entries to read in
    integer, intent(in)           :: n_max
    !> Number of entries found
    integer, intent(inout)        :: n_found
    !> On return, ix_start(i) holds the starting point of entry i
    integer, intent(inout)        :: ixs_start(n_max)
    !> On return, ix_end(i) holds the end point of entry i
    integer, intent(inout)        :: ixs_end(n_max)

    integer                       :: ix, ix_prev

    ix_prev = 0
    n_found = 0

    do while (n_found < n_max)

       ! Find the starting point of the next entry (a non-delimiter value)
       ix = verify(line(ix_prev+1:), delims)
       if (ix == 0) exit

       n_found            = n_found + 1
       ixs_start(n_found) = ix_prev + ix ! This is the absolute position in 'line'

       ! Get the end point of the current entry (next delimiter index minus one)
       ix = scan(line(ixs_start(n_found)+1:), delims) - 1

       if (ix == -1) then              ! If there is no last delimiter,
          ixs_end(n_found) = len(line) ! the end of the line is the endpoint
       else
          ixs_end(n_found) = ixs_start(n_found) + ix
       end if

       ix_prev = ixs_end(n_found) ! We continue to search from here
    end do

  end subroutine get_fields_string

  !> Replace text in string, inspired by
  !> http://fortranwiki.org/fortran/show/String_Functions
  pure function string_replace(string, text, replacement) result(outs)
    character(len=*), intent(in)  :: string
    character(len=*), intent(in)  :: text
    character(len=*), intent(in)  :: replacement
    character(len=:), allocatable :: outs
    integer, parameter :: buffer_space = 256
    character(len=(len(string)+buffer_space)) :: buffer
    integer                       :: i, nt, nr, len_outs

    buffer = string
    nt = len_trim(text)
    nr = len_trim(replacement)

    do
       i = index(buffer, text(:nt))
       if (i == 0) exit
       buffer = buffer(:i-1) // replacement(:nr) // trim(buffer(i+nt:))
    end do

    len_outs = len_trim(buffer)
    allocate(character(len=len_outs) :: outs)
    outs = buffer(1:len_outs)
  end function string_replace

  !> An inefficient routine to replace *^+- characters in a string
  subroutine to_simple_ascii(text, simple, charge)
    character(len=*), intent(in)    :: text
    character(len=*), intent(inout) :: simple
    integer, intent(out)            :: charge
    integer                         :: n
    logical                         :: in_brackets = .false.

    charge = 0
    simple = ""

    do n = 1, len_trim(text)
       select case (text(n:n))
       case ('(')
         in_brackets = .true.
         simple = trim(simple) // "_"
       case (')')
         in_brackets = .false.
       case ('*')
          simple = trim(simple) // "_star"
       case ('+')
          if (.not. in_brackets) then
            charge = charge + 1
          end if
          simple = trim(simple) // "_plus"
       case ('-')
          if (.not. in_brackets) then
            charge = charge - 1
          end if
          simple = trim(simple) // "_min"
       case ('^')
          simple = trim(simple) // "_hat"
       case ("'")
          simple = trim(simple) // "p"
       case default
          simple = trim(simple) // text(n:n)
       end select
    end do

    ! Handle some species separately
    if (simple == "e") charge = -1
  end subroutine to_simple_ascii

end module m_chemistry
