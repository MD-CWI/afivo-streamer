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
  integer, parameter :: rate_max_num_coeff = 4

  !> Basic chemical reaction type
  type reaction_t
     integer, allocatable  :: ix_in(:)            !< Index of input species
     integer, allocatable  :: ix_out(:)           !< Index of output species
     integer, allocatable  :: multiplicity_out(:) !< Multiplicity of output
     integer               :: rate_type           !< Type of reaction rate
     !> Type of reaction (e.g. ionization)
     integer               :: reaction_type = general_reaction
     real(dp)              :: rate_factor         !< Multiply rate by this factor
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

  !> Indicates a reaction with a field-dependent reaction rate
  integer, parameter :: rate_tabulated_field = 1

  !> Indicates a reaction with a constant reaction rate
  integer, parameter :: rate_analytic_constant = 2

  !> Indicates a reaction of the form c1 * (Td - c2)
  integer, parameter :: rate_analytic_linear = 3

  !> Indicates a reaction of the form c1 * exp(-(c2/(c3 + Td))**2)
  integer, parameter :: rate_analytic_exp_v1 = 4

  !> Indicates a reaction of the form c1 * exp(-(Td/c2)**2)
  integer, parameter :: rate_analytic_exp_v2 = 5

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
  public :: get_rates
  public :: get_derivatives

  public :: species_index

contains

  !> Initialize module and load chemical reactions
  subroutine chemistry_initialize(tree, cfg)
    use m_config
    use m_units_constants
    use m_table_data
    use m_transport_data
    use m_gas
    use m_advance_base, only: advance_num_states
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
          reactions(1)%rate_type = rate_tabulated_field
          reactions(1)%rate_factor = 1.0_dp
          reactions(1)%x_data = &
               LT_get_xdata(td_tbl%x_min, td_tbl%dx, td_tbl%n_points)
          reactions(1)%y_data = td_tbl%rows_cols(:, td_alpha) * &
               td_tbl%rows_cols(:, td_mobility) * reactions(1)%x_data * &
               Townsend_to_SI * gas_number_density
          reactions(1)%description = "e + M > e + e + M+"

          ! Attachment reaction
          reactions(2)%ix_in = [1]
          reactions(2)%ix_out = [3]
          reactions(2)%multiplicity_out = [1]
          reactions(2)%rate_type = rate_tabulated_field
          reactions(2)%rate_factor = 1.0_dp
          reactions(2)%x_data = &
               LT_get_xdata(td_tbl%x_min, td_tbl%dx, td_tbl%n_points)
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
    chemtbl = LT_create(table_min_townsend, table_max_townsend, table_size, i)

    i = 0
    do n = 1, n_reactions
       if (reactions(n)%rate_type == rate_tabulated_field) then
          i = i + 1
          reactions(n)%lookup_table_index = i
          call LT_set_col(chemtbl, i, reactions(n)%x_data, reactions(n)%y_data)
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
            n_copies=advance_num_states, ix=species_itree(n))
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
       write(*, "(I4,' ',A15,A)") n, reaction_names(reactions(n)%reaction_type), &
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
    real(dp), allocatable        :: eta(:), alpha(:)
    real(dp), allocatable        :: v(:), mu(:), diff(:)
    integer                      :: n, n_fields, i_elec
    integer                      :: my_unit

    if (gas_constant_density) then
       i_elec = species_index("e")
       n_fields = td_tbl%n_points

       allocate(fields(n_fields))
       fields = LT_get_xdata(td_tbl%x_min, td_tbl%dx, n_fields)

       allocate(rates(n_fields, n_reactions))
       allocate(eta(n_fields), alpha(n_fields))
       call get_rates(fields, rates, n_fields)

       eta(:)   = 0.0_dp
       alpha(:) = 0.0_dp

       do n = 1, n_reactions
          if (reactions(n)%reaction_type == attachment_reaction) then
             eta(:) = eta(:) + rates(:, n)
          else if (reactions(n)%reaction_type == ionization_reaction) then
             alpha(:) = alpha(:) + rates(:, n)
          end if
       end do

       allocate(diff(n_fields))
       diff = LT_get_col(td_tbl, td_diffusion, fields)

       allocate(mu(n_fields))
       allocate(v(n_fields))
       mu = LT_get_col(td_tbl, td_mobility, fields)
       v = mu * fields * Townsend_to_SI

       ! v(1) is zero, so extrapolate linearly
       eta(2:) = eta(2:) / v(2:)
       eta(1) = 2 * eta(2) - eta(3)
       alpha(2:) = alpha(2:) / v(2:)
       alpha(1) = 2 * alpha(2) - alpha(3)

       ! Write to a file
       open(newunit=my_unit, file=trim(fname), action="write")
       write(my_unit, "(A)") "# Description of columns"
       write(my_unit, "(A)") "# 1: E/N [Td]"
       write(my_unit, "(A)") "# 2: E [V/m]"
       write(my_unit, "(A)") "# 3: Electron mobility [m^2/(V s)]"
       write(my_unit, "(A)") "# 4: Electron diffusion [m^2/s]"
       write(my_unit, "(A)") "# 5: Townsend ioniz. coef. alpha [1/m]"
       write(my_unit, "(A)") "# 6: Townsend attach. coef. eta [1/m]"
       do n = 1, n_fields
          write(my_unit, *) fields(n), fields(n) * Townsend_to_SI, &
               mu(n) / gas_number_density, &
               diff(n) / gas_number_density, alpha(n), eta(n)
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

  !> Compute reaction rates
  subroutine get_rates(fields, rates, n_cells)
    integer, intent(in)   :: n_cells                     !< Number of cells
    real(dp), intent(in)  :: fields(n_cells)             !< The field (in Td) in the cells
    real(dp), intent(out) :: rates(n_cells, n_reactions) !< The reaction rates
    integer               :: n
    real(dp)              :: c0, c(rate_max_num_coeff)

    do n = 1, n_reactions
       ! A factor that the reaction rate is multiplied with, for example to
       ! account for a constant gas number density
       c0   = reactions(n)%rate_factor

       ! Coefficients for the reaction rate
       c(:) = reactions(n)%rate_data

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
       end select
    end do
  end subroutine get_rates

  !> Compute derivatives due to chemical reactions
  subroutine get_derivatives(dens, rates, derivs, n_cells)
    integer, intent(in)   :: n_cells
    real(dp), intent(in)  :: dens(n_cells, n_species)
    real(dp), intent(in)  :: rates(n_cells, n_reactions)
    real(dp), intent(out) :: derivs(n_cells, n_species)
    real(dp)              :: prate(n_cells, n_reactions)
    integer               :: n, i, ix

    derivs(:, :) = 0.0_dp

    ! Loop over reactions and add to derivative
    do n = 1, n_reactions
       ! Determine production rate of 'full' reaction
       prate(:, n) = rates(:, n) * &
            product(dens(:, tiny_react(n)%ix_in), dim=2)

       ! Input species are removed
       do i = 1, size(tiny_react(n)%ix_in)
          ix = tiny_react(n)%ix_in(i)
          derivs(:, ix) = derivs(:, ix) - prate(:, n)
       end do

       ! Output species are created
       do i = 1, size(tiny_react(n)%ix_out)
          ix = tiny_react(n)%ix_out(i)
          derivs(:, ix) = derivs(:, ix) + prate(:, n) * &
               tiny_react(n)%multiplicity_out(i)
       end do
    end do
  end subroutine get_derivatives

  !> Read reactions from a file
  subroutine read_reactions(filename, read_success)
    character(len=*), intent(in) :: filename
    logical, intent(out)         :: read_success
    character(len=string_len)    :: line
    character(len=50)            :: data_value(max_num_reactions)
    character(len=50)            :: reaction(max_num_reactions)
    character(len=50)            :: how_to_get(max_num_reactions)
    type(reaction_t)             :: new_reaction
    integer                      :: my_unit
    integer, parameter           :: n_fields = 3
    integer                      :: i0(n_fields), i1(n_fields)
    integer                      :: n, n_found

    open(newunit=my_unit, file=filename, action="read")

    n_reactions  = 0
    read_success = .false.

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

       call get_fields_string(line, ",", n_fields, n_found, i0, i1)

       if (n_found /= n_fields) then
          print *, trim(line)
          error stop "Invalid chemistry syntax"
       end if

       if (n_reactions >= max_num_reactions) &
            error stop "Too many reactions, increase max_num_reactions"

       n_reactions             = n_reactions + 1
       reaction(n_reactions)   = line(i0(1):i1(1))
       how_to_get(n_reactions) = line(i0(2):i1(2))
       data_value(n_reactions) = line(i0(3):i1(3))
    end do

    ! Close the file (so that we can re-open it for reading data)
    close(my_unit)

    ! Now parse the reactions
    do n = 1, n_reactions
       call parse_reaction(trim(reaction(n)), new_reaction)
       new_reaction%description = trim(reaction(n))

       select case (how_to_get(n))
       case ("field_table")
          ! Reaction data should be present in the same file
          call read_reaction_table(filename, &
               trim(data_value(n)), new_reaction)
       case ("constant")
          new_reaction%rate_type = rate_analytic_constant
          read(data_value(n), *) new_reaction%rate_data(1)
       case ("linear")
          new_reaction%rate_type = rate_analytic_linear
          read(data_value(n), *) new_reaction%rate_data(1:2)
       case ("exp_v1")
          new_reaction%rate_type = rate_analytic_exp_v1
          read(data_value(n), *) new_reaction%rate_data(1:3)
       case ("exp_v2")
          new_reaction%rate_type = rate_analytic_exp_v2
          read(data_value(n), *) new_reaction%rate_data(1:2)
       case default
          print *, "Unknown rate type: ", trim(how_to_get(n))
          print *, "For reaction:      ", trim(reaction(n))
          print *, "In file:           ", trim(filename)
          error stop
       end select

       reactions(n) = new_reaction
    end do

    if (n_reactions > 0) read_success = .true.
998 return

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
  subroutine parse_reaction(reaction_text, reaction)
    use m_gas
    character(len=*), intent(in)  :: reaction_text
    type(reaction_t), intent(out) :: reaction
    integer, parameter            :: max_components = 100
    character(len=comp_len)       :: component
    integer                       :: i, ix, n, n_found, multiplicity
    integer                       :: n_in, n_out
    integer                       :: i0(max_components)
    integer                       :: i1(max_components)
    integer                       :: ix_in(max_components)
    integer                       :: ix_out(max_components)
    integer                       :: multiplicity_out(max_components)
    logical                       :: left_side
    real(dp)                      :: rfactor

    call get_fields_string(reaction_text, " ", max_components, n_found, i0, i1)

    left_side = .true.
    n_in      = 0
    n_out     = 0
    rfactor   = 1.0_dp

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

  !> An inefficient routine to replace *^+- characters in a string
  subroutine to_simple_ascii(text, simple, charge)
    character(len=*), intent(in)    :: text
    character(len=*), intent(inout) :: simple
    integer, intent(out)            :: charge
    integer                         :: n

    charge = 0
    simple = ""

    do n = 1, len_trim(text)
       select case (text(n:n))
       case ('*')
          simple = trim(simple) // "_star"
       case ('+')
          charge = charge + 1
          simple = trim(simple) // "_plus"
       case ('-')
          charge = charge - 1
          simple = trim(simple) // "_min"
       case ('^')
          simple = trim(simple) // "_hat"
       case default
          simple = trim(simple) // text(n:n)
       end select
    end do

    ! Handle some species separately
    if (simple == "e") charge = -1
  end subroutine to_simple_ascii

end module m_chemistry
