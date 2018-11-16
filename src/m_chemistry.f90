#include "../afivo/src/cpp_macros.h"
module m_chemistry
  use m_streamer
  use m_af_all
  use m_lookup_table

  implicit none
  private

  type reaction_t
     integer, allocatable  :: ix_in(:)
     integer, allocatable  :: ix_out(:)
     integer, allocatable  :: multiplicity_out(:)
     integer               :: rate_type
     real(dp), allocatable :: tbl(:, :)
     character(len=50)      :: description
  end type reaction_t

  integer, parameter :: constant_rate        = 1
  integer, parameter :: field_dependent_rate = 2
  integer, parameter :: max_num_species = 100
  integer, parameter :: max_num_reactions = 100

  integer :: n_species = 0
  integer :: n_reactions = 0

  character(len=10)    :: species_list(max_num_species)
  integer              :: species_ix(max_num_species)
  type(reaction_t)     :: reaction_list(max_num_reactions)
  type(lookup_table_t) :: chemtbl

  integer  :: rate_table_size   = 1000
  real(dp) :: rate_min_townsend = 0.0
  real(dp) :: rate_max_townsend = 1e3_dp

  public :: chemistry_init

contains

  subroutine chemistry_init(tree, cfg)
    use m_config
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg
    integer                    :: n
    character(len=ST_slen)     :: reaction_file

    call CFG_add_get(cfg, "chemistry%rate_table_size", rate_table_size, &
         "Size of the lookup table for reaction rates")
    call CFG_add_get(cfg, "chemistry%rate_min_townsend", rate_min_townsend, &
         "Minimal field (in Td) for the rate coeff. lookup table")
    call CFG_add_get(cfg, "chemistry%rate_max_townsend", rate_max_townsend, &
         "Maximal field (in Td) for the rate coeff. lookup table")

    reaction_file = "reactions.txt"
    call CFG_add_get(cfg, "chemistry%reaction_file", reaction_file, &
         "File with a list of reactions")

    call read_reactions(reaction_file)

    chemtbl = LT_create(rate_min_townsend, rate_max_townsend, &
         rate_table_size, n_reactions)

    do n = 1, n_reactions
       select case (reaction_list(n)%rate_type)
       case (constant_rate, field_dependent_rate)
          call LT_set_col(chemtbl, n, &
               reaction_list(n)%tbl(:, 1), reaction_list(n)%tbl(:, 2))
       case default
          error stop "Unknown type of reaction rate"
       end select
    end do

    do n = 1, n_species
       call af_add_cc_variable(tree, trim(species_list(n)), &
            n_copies=2, ix=species_ix(n))
    end do

  end subroutine chemistry_init

  subroutine get_derivatives(box, derivs)
    type(af_t), intent(in) :: box
    real(dp) :: derives(DTIMES(box%n_cell))

  end subroutine get_derivatives

  subroutine read_reactions(filename)
    character(len=*), intent(in) :: filename
    character(len=ST_slen)       :: line, data_value
    character(len=50)            :: reaction, how_to_get
    type(reaction_t)             :: new_reaction
    integer                      :: my_unit
    integer, parameter           :: n_fields = 3
    integer                      :: i0(n_fields), i1(n_fields)
    integer                      :: n_found

    open(newunit=my_unit, file=filename, action="read")

    n_reactions = 0

    do
       read(my_unit, "(A)", end=999) line
       call get_fields_string(line, ",", n_fields, n_found, i0, i1)

       if (n_found /= n_fields) then
          print *, trim(line)
          error stop "Invalid chemistry syntax"
       end if

       reaction = line(i0(1):i1(1))
       how_to_get = line(i0(2):i1(2))
       data_value = line(i0(3):i1(3))

       call get_reaction(trim(reaction), new_reaction)
       new_reaction%description = trim(reaction)
       ! print *, trim(reaction)
       ! print *, species_list(new_reaction%ix_in)
       ! print *, species_list(new_reaction%ix_out)
       ! print *, new_reaction%multiplicity_out

       select case (how_to_get)
       case ("table")
          call read_reaction_table(trim(data_value), new_reaction)
       case ("constant")
          call read_reaction_constant(trim(data_value), new_reaction)
       end select

       n_reactions                = n_reactions + 1
       reaction_list(n_reactions) = new_reaction
    end do
    close(my_unit)

999 continue
  end subroutine read_reactions

  subroutine read_reaction_constant(text, rdata)
    character(len=*), intent(in)    :: text
    type(reaction_t), intent(inout) :: rdata
    real(dp)                        :: tmp

    rdata%rate_type = constant_rate
    read(text, *) tmp

    allocate(rdata%tbl(2, 2))
    rdata%tbl(:, 1) = [rate_min_townsend, rate_max_townsend]
    rdata%tbl(:, 2) = [tmp, tmp]
  end subroutine read_reaction_constant

  subroutine read_reaction_table(filename, rdata)
    character(len=*), intent(in)    :: filename
    type(reaction_t), intent(inout) :: rdata
    integer, parameter              :: max_num_rows = 1000
    real(dp)                        :: tmptbl(2, max_num_rows)
    character(len=ST_slen)          :: line
    integer                         :: n_rows, my_unit

    n_rows = 0

    open(newunit=my_unit, file=filename, action="read")
    do
       read(my_unit, "(A)", end=999) line
       line = adjustl(line)
       if (line(1:1) == '#') cycle
       n_rows = n_rows + 1
       read(line, *) tmptbl(:, n_rows)
    end do
999 close(my_unit)

    rdata%rate_type = field_dependent_rate
    rdata%tbl       = tmptbl(:, 1:n_rows)
  end subroutine read_reaction_table

  subroutine get_reaction(reaction_text, reaction)
    character(len=*), intent(in)  :: reaction_text
    type(reaction_t), intent(out) :: reaction
    integer, parameter            :: max_components = 100
    character(len=10)             :: component
    integer                       :: i, ix, n, n_found, multiplicity
    integer                       :: n_in, n_out
    integer                       :: i0(max_components)
    integer                       :: i1(max_components)
    integer                       :: ix_in(max_components)
    integer                       :: ix_out(max_components)
    integer                       :: multiplicity_out(max_components)
    logical                       :: left_side

    call get_fields_string(reaction_text, " ", max_components, n_found, i0, i1)

    left_side = .true.
    n_in      = 0
    n_out     = 0

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
          print *, multiplicity
          component = component(2:)
       else
          multiplicity = 1
       end if

       ix = species_index(component)
       if (ix == -1) then
          n_species        = n_species + 1
          ix               = n_species
          call to_simple_ascii(trim(component), species_list(ix))
       end if

       if (left_side) then
          do i = 1, multiplicity
             n_in = n_in + 1
             ix_in(n_in) = ix
          end do
       else
          ! Check if species is already present in right-hand side
          do i = 1, n_out
             if (ix == ix_out(i)) then
                multiplicity_out(i) = multiplicity_out(i) + multiplicity
                exit
             end if
          end do

          ! If not already present, add the species
          if (i == n_out+1) then
             n_out = n_out + 1
             ix_out(n_out) = ix
             multiplicity_out(n_out) = multiplicity
          end if
       end if
    end do

    reaction%ix_in = ix_in(1:n_in)
    reaction%ix_out = ix_out(1:n_out)
    reaction%multiplicity_out = multiplicity_out(1:n_out)
  end subroutine get_reaction

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
  subroutine to_simple_ascii(text, simple)
    character(len=*), intent(in)   :: text
    character(len=50), intent(out) :: simple
    integer                        :: n

    simple = ""
    do n = 1, len_trim(text)
       select case (text(n:n))
       case ('*')
          simple = trim(simple) // "_star"
       case ('+')
          simple = trim(simple) // "_plus"
       case ('-')
          simple = trim(simple) // "_min"
       case ('^')
          simple = trim(simple) // "_hat"
       case default
          simple = trim(simple) // text(n:n)
       end select
    end do
  end subroutine to_simple_ascii

end module m_chemistry
