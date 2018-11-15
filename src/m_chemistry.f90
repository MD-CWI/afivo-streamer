#include "../afivo/src/cpp_macros.h"
module m_chemistry
  use m_streamer
  use m_af_all

  implicit none
  private

  type reaction_t
     integer, allocatable :: ix_in(:)
     integer, allocatable :: ix_out(:)
     integer, allocatable :: multiplicity_out(:)
     ! lookup table
  end type reaction_t

  character(len=*), parameter :: undefined_string = "__UNDEFINED"

  integer, parameter :: max_num_species = 100
  integer, parameter :: max_num_reactions = 100

  integer :: n_species = 0
  integer :: n_reactions = 0

  character(len=10) :: species_list(max_num_species)
  type(reaction_t) :: reaction_list(max_num_reactions)


  public :: chemistry_init

contains

  subroutine chemistry_init(tree, cfg)
    use m_config
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg

    character(len=ST_slen) :: reaction_file

    reaction_file = undefined_string
    call CFG_add_get(cfg, "chemistry%reaction_file", reaction_file, &
         "File with a list of reactions")

    if (reaction_file /= undefined_string) then
       call read_reactions(reaction_file)
    else
       n_species = 0
       n_reactions = 0
    end if

    stop

  end subroutine chemistry_init

  subroutine read_reactions(filename)
    character(len=*), intent(in) :: filename
    character(len=ST_slen)       :: line, data_value
    character(len=50)            :: reaction, how_to_get
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

       n_reactions = n_reactions + 1
       reaction = line(i0(1):i1(1))
       how_to_get = line(i0(2):i1(2))
       data_value = line(i0(3):i1(3))

       call add_reaction(trim(reaction))
       print *, trim(reaction)
       print *, trim(how_to_get)
       print *, trim(data_value)
    end do
    close(my_unit)

999 continue
  end subroutine read_reactions

  subroutine add_reaction(reaction)
    character(len=*), intent(in) :: reaction
    integer, parameter           :: max_components = 10
    character(len=10)            :: component
    integer                      :: i, ix, n, n_found, multiplicity
    integer                      :: n_in, n_out
    integer                      :: i0(max_components)
    integer                      :: i1(max_components)
    integer                      :: ix_in(max_components)
    integer                      :: ix_out(max_components)
    integer                      :: multiplicity_out(max_components)
    logical                      :: left_side

    call get_fields_string(reaction, " ", max_components, n_found, i0, i1)

    left_side = .true.
    n_in = 0
    n_out = 0

    do n = 1, n_found
       component = reaction(i0(n):i1(n))

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
          species_list(ix) = trim(component)
       end if

       if (left_side) then
          do i = 1, multiplicity
             n_in = n_in + 1
             ix_in(n_in) = ix
          end do
       else
          n_out = n_out + 1
          ix_out(n_out) = ix
          multiplicity_out(n_out) = multiplicity
       end if
    end do

    print *, reaction
    print *, species_list(ix_in(1:n_in))
    print *, species_list(ix_out(1:n_out)), multiplicity_out(ix_out(1:n_out))
  end subroutine add_reaction

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

end module m_chemistry
