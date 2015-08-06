! Module that allows working with a configuration file.
! Don't use the functions in this module in a print statement!
module m_config

  implicit none
  private

  ! Types
  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: CFG_n_types = 4
  integer, parameter :: CFG_int_type = 1, CFG_real_type = 2, &
       CFG_char_type = 3, CFG_logic_type = 4
  character(len=5), parameter :: CFG_type_names(CFG_n_types) = &
       (/"int  ", "real ", "char ", "logic"/)

  ! String lengths
  integer, parameter :: tiny_len = 20
  integer, parameter :: name_len = 40
  integer, parameter :: line_len = 400

  ! Maximum length of a variable
  integer, parameter :: max_var_len = 20

  type CFG_var_t
     character(len=name_len)              :: p_name
     character(len=line_len)              :: comment
     integer                              :: p_type
     integer                              :: p_size
     logical                              :: dyn_size

     real(dp), allocatable                :: real_data(:)
     integer, allocatable                 :: int_data(:)
     character(len=name_len), allocatable :: char_data(:)
     logical, allocatable                 :: logic_data(:)
  end type CFG_var_t

  type CFG_t
     logical                           :: sorted      = .false.
     integer                           :: n_vars      = 0
     type(CFG_var_t), allocatable      :: vars(:)
  end type CFG_t

  interface CFG_add
     module procedure :: add_real, add_real_array, add_int, add_int_array, &
          add_char, add_char_array, add_logical, add_logical_array
  end interface CFG_add

  interface CFG_get
     module procedure :: get_real, get_int, get_logic, get_string, &
          get_real_array, get_int_array, get_char_array, get_logical_array
  end interface CFG_get

  public :: CFG_add
  public :: CFG_get
  public :: CFG_get_size
  public :: CFG_get_type
  public :: CFG_sort
  public :: CFG_write
  public :: CFG_read_file

  public :: CFG_fget_real
  public :: CFG_fget_int
  public :: CFG_fget_logic
  public :: CFG_fget_string
  public :: CFG_fget_size

  ! Public types
  public :: CFG_t
  public :: CFG_int_type, CFG_real_type, CFG_char_type, CFG_logic_type

contains

  subroutine handle_error(err_string)
    character(len=*), intent(in) :: err_string

    print *, "The following fatal error occured in m_config:"
    print *, trim(err_string)
    ! Gnu extension to get a backtrace
    call abort()
  end subroutine handle_error

  !> Update the variables in the configartion with the values found in 'filename'
  subroutine CFG_read_file(cfg, filename)
    type(CFG_t), intent(inout) :: cfg
    character(len=*), intent(in) :: filename

    integer                      :: ioState, equalSignPos, lineEnd
    integer                      :: n, ix, nL
    integer                      :: startIxs(max_var_len), endIxs(max_var_len), nEntries
    character(len=tiny_len)      :: configFormat
    character(len=name_len)      :: p_name
    character(len=line_len)      :: line, err_string

    open(UNIT = 1, FILE = filename, STATUS = "OLD", ACTION = "READ", ERR = 999, IOSTAT = ioState)
    nL = 0

    ! Set the line format to read, only depends on line_len currently
    write(configFormat, FMT = "(I6)") line_len
    configFormat = "(A" // trim(adjustl(configFormat)) // ")"

    do
       read( UNIT = 1, FMT = configFormat, ERR = 999, end = 666) line; nL = nL + 1
       ! Line format should be "variable = value(s) [#Comment]"

       lineEnd = scan(line, '#') - 1       ! Don't use the part of the line after '#'
       if (lineEnd == -1) lineEnd = line_len
       line = line(1:lineEnd)

       equalSignPos = scan(line, '=')      ! Locate the '=' sign, if there is no such sign then skip the line
       if (equalSignPos == 0) cycle

       p_name = line(1 : equalSignPos - 1)  ! Set variable name
       p_name = adjustl(p_name)              ! Remove leading blanks

       line = line(equalSignPos + 1 : lineEnd)   ! Set line to the values behind the '=' sign
       line = adjustl(line)                      ! Remove leading blanks

       ! Find variable corresponding to name in file
       call get_var_index(cfg, p_name, ix)

       if (ix <= 0) then
          call handle_error("CFG_read_file: variable [" // trim(p_name) // &
               "] from [" //  filename // "] could not be updated")
          return
       end if

       if (cfg%vars(ix)%dyn_size) then
          ! Get the start and end positions of the line content, and the number of entries
          call get_fields_string(line, " ,'"""//char(9), max_var_len, nEntries, startIxs, endIxs)
          if (cfg%vars(ix)%p_size /= nEntries) then
             call resize_storage(cfg, ix, cfg%vars(ix)%p_type, nEntries)
             cfg%vars(ix)%p_size = nEntries
          end if

          do n = 1, nEntries
             select case (cfg%vars(ix)%p_type)
             case (CFG_int_type)
                read(line(startIxs(n):endIxs(n)), *, ERR = 999) cfg%vars(ix)%int_data(n)
             case (CFG_real_type)
                read(line(startIxs(n):endIxs(n)), *, ERR = 999) cfg%vars(ix)%real_data(n)
             case (CFG_char_type)
                cfg%vars(ix)%char_data(n) = trim(line(startIxs(n):endIxs(n)))
             case (CFG_logic_type)
                read(line(startIxs(n):endIxs(n)), *, ERR = 999) cfg%vars(ix)%logic_data(n)
             end select
          end do
       else                  ! Fixed size parameter
          select case (cfg%vars(ix)%p_type)
          case (CFG_int_type)
             read(line, *, ERR = 999, end = 998) cfg%vars(ix)%int_data
          case (CFG_real_type)
             read(line, *, ERR = 999, end = 998) cfg%vars(ix)%real_data
          case (CFG_char_type)
             cfg%vars(ix)%char_data = trim(line)
          case (CFG_logic_type)
             read(line, *, ERR = 999, end = 998) cfg%vars(ix)%logic_data
          end select
       end if
    end do

666 continue ! Routine ends here if the end of "filename" is reached
    close(UNIT = 1, STATUS = "KEEP", ERR = 999, IOSTAT = ioState)
    return

    ! The routine only ends up here through an error
998 continue
    call handle_error("CFG_read_file: Not enough values for variable [" &
         // p_name // "] in " // filename)
    return

999 continue
    write(err_string, *) "ioState = ", ioState, " while reading from ", &
         filename, " at line ", nL
    call handle_error("CFG_read_file:" // err_string)
    return

  end subroutine CFG_read_file

  subroutine resize_storage(cfg, ix, p_type, p_size)
    type(CFG_t), intent(inout) :: cfg
    integer, intent(in) :: ix, p_type, p_size

    select case (p_type)
    case (CFG_int_type)
       deallocate( cfg%vars(ix)%int_data )
       allocate( cfg%vars(ix)%int_data(p_size) )
    case (CFG_logic_type)
       deallocate( cfg%vars(ix)%logic_data )
       allocate( cfg%vars(ix)%logic_data(p_size) )
    case (CFG_real_type)
       deallocate( cfg%vars(ix)%real_data )
       allocate( cfg%vars(ix)%real_data(p_size) )
    case (CFG_char_type)
       deallocate( cfg%vars(ix)%char_data )
       allocate( cfg%vars(ix)%char_data(p_size) )
    end select
  end subroutine resize_storage

  !> Return the index of the variable with name 'p_name', or -1 if not found.
  subroutine get_var_index(cfg, p_name, ix)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: p_name
    integer, intent(inout)       :: ix
    integer                      :: i

    if (cfg%sorted) then
       call bin_search_var(cfg, p_name, ix)
    else
       ix = -1
       do i = 1, cfg%n_vars
          if (cfg%vars(i)%p_name == p_name) then
             ix = i
             exit
          end if
       end do
    end if
  end subroutine get_var_index

  !> Performa a binary search for the variable 'p_name'
  subroutine bin_search_var(cfg, p_name, ix)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: p_name
    integer, intent(inout)       :: ix

    integer                       :: i_min, i_max, i_mid

    i_min = 1
    i_max = cfg%n_vars
    ix    = - 1

    do while (i_min < i_max)
       i_mid = i_min + (i_max - i_min) / 2
       if ( llt(cfg%vars(i_mid)%p_name, p_name) ) then
          i_min = i_mid + 1
       else
          i_max = i_mid
       end if
    end do

    ! If not found, bin_search_var is not set here, and stays -1
    if (i_max == i_min .and. cfg%vars(i_min)%p_name == p_name) then
       ix = i_min
    else
       ix = -1
    end if
  end subroutine bin_search_var

  !> Helper routine to create variables. This is useful because a lot of the same
  !! code is executed for the different types of variables.
  subroutine prepare_store_var(cfg, p_name, p_type, p_size, comment, dyn_size)
    type(CFG_t), intent(inout) :: cfg
    character(len=*), intent(in)  :: p_name, comment
    integer, intent(in)           :: p_type, p_size
    logical, intent(in), optional :: dyn_size
    integer :: n_vars

    call ensure_free_storage(cfg)

    cfg%sorted               = .false.
    n_vars                   = cfg%n_vars + 1
    cfg%n_vars               = n_vars
    cfg%vars(n_vars)%p_name  = p_name
    cfg%vars(n_vars)%comment = comment
    cfg%vars(n_vars)%p_type  = p_type
    cfg%vars(n_vars)%p_size  = p_size

    if (present(dyn_size)) then
       cfg%vars(n_vars)%dyn_size = dyn_size
    else
       cfg%vars(n_vars)%dyn_size = .false.
    end if

    select case (cfg%vars(n_vars)%p_type)
    case (CFG_int_type)
       allocate( cfg%vars(n_vars)%int_data(p_size) )
    case (CFG_real_type)
       allocate( cfg%vars(n_vars)%real_data(p_size) )
    case (CFG_char_type)
       allocate( cfg%vars(n_vars)%char_data(p_size) )
    case (CFG_logic_type)
       allocate( cfg%vars(n_vars)%logic_data(p_size) )
    end select
  end subroutine prepare_store_var

  subroutine prepare_get_var(cfg, p_name, p_type, p_size, ix)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: p_name
    integer, intent(in)          :: p_type, p_size
    integer, intent(out)         :: ix
    character(len=line_len)      :: err_string

    call get_var_index(cfg, p_name, ix)

    if (ix == -1) then
       call handle_error(&
            "CFG_get(...): variable ["//p_name//"] not found")
    else if (cfg%vars(ix)%p_type /= p_type) then
       write(err_string, fmt="(A)") "CFG_get(...): variable [" &
            // p_name // "] has different type (" // &
            trim(CFG_type_names(cfg%vars(ix)%p_type)) // &
            ") than requested (" // trim(CFG_type_names(p_type)) // ")"
       call handle_error(err_string)
    else if (cfg%vars(ix)%p_size /= p_size) then
       write(err_string, fmt="(A,I0,A,I0,A)") "CFG_get(...): variable [" &
            // p_name // "] has different size (", cfg%vars(ix)%p_size, &
            ") than requested (", p_size, ")"
       call handle_error(err_string)
    end if
  end subroutine prepare_get_var

  subroutine add_real(cfg, p_name, real_data, comment)
    type(CFG_t), intent(inout) :: cfg
    character(len=*), intent(in)  :: p_name, comment
    real(dp), intent(in)  :: real_data
    call prepare_store_var(cfg, p_name, CFG_real_type, 1, comment)
    cfg%vars(cfg%n_vars)%real_data(1) = real_data
  end subroutine add_real

  subroutine add_real_array(cfg, p_name, real_data, comment, dyn_size)
    type(CFG_t), intent(inout) :: cfg
    character(len=*), intent(in)  :: p_name, comment
    real(dp), intent(in)  :: real_data(:)
    logical, intent(in), optional :: dyn_size
    call prepare_store_var(cfg, p_name, CFG_real_type, size(real_data), comment, dyn_size)
    cfg%vars(cfg%n_vars)%real_data = real_data
  end subroutine add_real_array

  subroutine add_int(cfg, p_name, int_data, comment)
    type(CFG_t), intent(inout) :: cfg
    character(len=*), intent(in)  :: p_name, comment
    integer, intent(in)  :: int_data
    call prepare_store_var(cfg, p_name, CFG_int_type, 1, comment)
    cfg%vars(cfg%n_vars)%int_data(1) = int_data
  end subroutine add_int

  subroutine add_int_array(cfg, p_name, int_data, comment, dyn_size)
    type(CFG_t), intent(inout) :: cfg
    character(len=*), intent(in)  :: p_name, comment
    integer, intent(in)           :: int_data(:)
    logical, intent(in), optional :: dyn_size
    call prepare_store_var(cfg, p_name, CFG_int_type, size(int_data), comment, dyn_size)
    cfg%vars(cfg%n_vars)%int_data = int_data
  end subroutine add_int_array

  subroutine add_char(cfg, p_name, char_data, comment)
    type(CFG_t), intent(inout) :: cfg
    character(len=*), intent(in)  :: p_name, comment, char_data
    call prepare_store_var(cfg, p_name, CFG_char_type, 1, comment)
    cfg%vars(cfg%n_vars)%char_data(1) = char_data
  end subroutine add_char

  subroutine add_char_array(cfg, p_name, char_data, comment, dyn_size)
    type(CFG_t), intent(inout) :: cfg
    character(len=*), intent(in)  :: p_name, comment, char_data(:)
    logical, intent(in), optional :: dyn_size
    call prepare_store_var(cfg, p_name, CFG_char_type, size(char_data), comment, dyn_size)
    cfg%vars(cfg%n_vars)%char_data = char_data
  end subroutine add_char_array

  subroutine add_logical(cfg, p_name, logic_data, comment)
    type(CFG_t), intent(inout) :: cfg
    character(len=*), intent(in)  :: p_name, comment
    logical, intent(in)           :: logic_data
    call prepare_store_var(cfg, p_name, CFG_logic_type, 1, comment)
    cfg%vars(cfg%n_vars)%logic_data(1) = logic_data
  end subroutine add_logical

  subroutine add_logical_array(cfg, p_name, logic_data, comment, dyn_size)
    type(CFG_t), intent(inout) :: cfg
    character(len=*), intent(in)  :: p_name, comment
    logical, intent(in)           :: logic_data(:)
    logical, intent(in), optional :: dyn_size
    call prepare_store_var(cfg, p_name, CFG_logic_type, size(logic_data), comment, dyn_size)
    cfg%vars(cfg%n_vars)%logic_data = logic_data
  end subroutine add_logical_array

  subroutine get_real_array(cfg, p_name, real_data)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: p_name
    real(dp), intent(inout)      :: real_data(:)
    integer                      :: ix
    call prepare_get_var(cfg, p_name, CFG_real_type, size(real_data), ix)
    real_data = cfg%vars(ix)%real_data
  end subroutine get_real_array

  subroutine get_int_array(cfg, p_name, int_data)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: p_name
    integer, intent(inout)       :: int_data(:)
    integer                      :: ix
    call prepare_get_var(cfg, p_name, CFG_int_type, size(int_data), ix)
    int_data    = cfg%vars(ix)%int_data
  end subroutine get_int_array

  subroutine get_char_array(cfg, p_name, char_data)
    type(CFG_t), intent(in) :: cfg
    character(len=*), intent(in)     :: p_name
    character(len=*), intent(inout)  :: char_data(:)
    integer :: ix
    call prepare_get_var(cfg, p_name, CFG_char_type, size(char_data), ix)
    char_data   = cfg%vars(ix)%char_data
  end subroutine get_char_array

  subroutine get_logical_array(cfg, p_name, logic_data)
    type(CFG_t), intent(in) :: cfg
    character(len=*), intent(in)     :: p_name
    logical, intent(inout)           :: logic_data(:)
    integer :: ix
    call prepare_get_var(cfg, p_name, CFG_logic_type, size(logic_data), ix)
    logic_data   = cfg%vars(ix)%logic_data
  end subroutine get_logical_array

  subroutine get_real(cfg, p_name, res)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: p_name
    real(dp), intent(out)        :: res
    integer                      :: ix
    call prepare_get_var(cfg, p_name, CFG_real_type, 1, ix)
    res = cfg%vars(ix)%real_data(1)
  end subroutine get_real

  subroutine get_int(cfg, p_name, res)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: p_name
    integer, intent(inout)       :: res
    integer                      :: ix
    call prepare_get_var(cfg, p_name, CFG_int_type, 1, ix)
    res = cfg%vars(ix)%int_data(1)
  end subroutine get_int

  subroutine get_logic(cfg, p_name, res)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: p_name
    logical, intent(out)         :: res
    integer                      :: ix
    call prepare_get_var(cfg, p_name, CFG_logic_type, 1, ix)
    res = cfg%vars(ix)%logic_data(1)
  end subroutine get_logic

  subroutine get_string(cfg, p_name, res)
    type(CFG_t), intent(in)       :: cfg
    character(len=*), intent(in)  :: p_name
    character(len=*), intent(out) :: res
    integer                       :: ix
    call prepare_get_var(cfg, p_name, CFG_char_type, 1, ix)
    res = cfg%vars(ix)%char_data(1)
  end subroutine get_string

  subroutine CFG_get_size(cfg, p_name, res)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: p_name
    integer, intent(out)         :: res
    integer                      :: ix

    call get_var_index(cfg, p_name, ix)
    if (ix /= -1) then
       res = cfg%vars(ix)%p_size
    else
       res = -1
       call handle_error("CFG_get_size: variable ["//p_name//"] not found")
    end if
  end subroutine CFG_get_size

  real(dp) function CFG_fget_real(cfg, p_name)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: p_name
    integer                      :: ix
    call prepare_get_var(cfg, p_name, CFG_real_type, 1, ix)
    CFG_fget_real = cfg%vars(ix)%real_data(1)
  end function CFG_fget_real

  integer function CFG_fget_int(cfg, p_name)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: p_name
    integer                      :: ix
    call prepare_get_var(cfg, p_name, CFG_int_type, 1, ix)
    CFG_fget_int = cfg%vars(ix)%int_data(1)
  end function CFG_fget_int

  logical function CFG_fget_logic(cfg, p_name)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: p_name
    integer                      :: ix
    call prepare_get_var(cfg, p_name, CFG_logic_type, 1, ix)
    CFG_fget_logic = cfg%vars(ix)%logic_data(1)
  end function CFG_fget_logic

  character(len=name_len) function CFG_fget_string(cfg, p_name)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: p_name
    integer                      :: ix
    call prepare_get_var(cfg, p_name, CFG_char_type, 1, ix)
    CFG_fget_string = cfg%vars(ix)%char_data(1)
  end function CFG_fget_string

  integer function CFG_fget_size(cfg, p_name)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: p_name
    integer                      :: ix

    call get_var_index(cfg, p_name, ix)
    if (ix /= -1) then
       CFG_fget_size = cfg%vars(ix)%p_size
    else
       CFG_fget_size = -1
       call handle_error("CFG_get_size: variable ["//p_name//"] not found")
    end if
  end function CFG_fget_size

  subroutine CFG_get_type(cfg, p_name, res)
    type(CFG_t), intent(in)      :: cfg
    character(len=*), intent(in) :: p_name
    integer, intent(out)         :: res
    integer                      :: ix

    call get_var_index(cfg, p_name, ix)

    if (ix /= -1) then
       res = cfg%vars(ix)%p_type
    else
       res = -1
       call handle_error("get_type: variable ["//p_name//"] not found")
    end if
  end subroutine CFG_get_type

  subroutine CFG_write(cfg, filename)
    use iso_fortran_env
    type(CFG_t), intent(in) :: cfg
    character(len=*), intent(in) :: filename

    integer                      :: i, j, ioState, myUnit
    character(len=tiny_len)      :: nameFormat
    character(len=line_len)      :: err_string

    write(nameFormat, FMT = "(A,I0,A)") "(A,A", name_len, ",A)"

    if (filename == "stdout") then
       myUnit = output_unit
    else
       myUnit = 333
       open(myUnit, FILE = filename, ACTION = "WRITE", ERR = 999, IOSTAT = ioState)
    end if

    write(myUnit, ERR = 999, FMT = "(A)") " ##############################################"
    write(myUnit, ERR = 999, FMT = "(A)") " ###          Configuration file            ###"
    write(myUnit, ERR = 999, FMT = "(A)") " ##############################################"
    write(myUnit, ERR = 999, FMT = "(A)") ""

    do i = 1, cfg%n_vars
       write(myUnit, ERR = 999, FMT = "(A,A,A)") " # ", trim(cfg%vars(i)%comment), ":"
       write(myUnit, ADVANCE = "NO", ERR = 999, FMT = nameFormat) " ", cfg%vars(i)%p_name, " = "

       select case(cfg%vars(i)%p_type)
       case (CFG_int_type)
          do j = 1, cfg%vars(i)%p_size
             write(myUnit, ADVANCE = "NO", ERR = 999, FMT = "(I10, A) ") cfg%vars(i)%int_data(j), "  "
          end do
       case (CFG_real_type)
          do j = 1, cfg%vars(i)%p_size
             write(myUnit, ADVANCE = "NO", ERR = 999, FMT = "(E10.4, A) ") cfg%vars(i)%real_data(j), "  "
          end do
       case (CFG_char_type)
          do j = 1, cfg%vars(i)%p_size
             write(myUnit, ADVANCE = "NO", ERR = 999, FMT = "(A, A)") &
                  & '"' // trim(cfg%vars(i)%char_data(j)) // '"', "  "
          end do
       case (CFG_logic_type)
          do j = 1, cfg%vars(i)%p_size
             write(myUnit, ADVANCE = "NO", ERR = 999, FMT = "(L10, A) ") cfg%vars(i)%logic_data(j), "  "
          end do
       end select
       write(myUnit, ERR = 999, FMT = "(A)") ""
       write(myUnit, ERR = 999, FMT = "(A)") ""
    end do

    if (myUnit /= output_unit) close(myUnit, ERR = 999, IOSTAT = ioState)
    return

999 continue ! If there was an error, the routine will end here
    write(err_string, *) "CFG_write error: ioState = ", ioState, " while writing to ", filename
    call handle_error(err_string)

  end subroutine CFG_write

  subroutine ensure_free_storage(cfg)
    type(CFG_t), intent(inout) :: cfg
    integer, parameter         :: min_dyn_size = 100
    integer                    :: cur_size, new_size
    type(CFG_var_t), allocatable :: cfg_copy(:)

    if (allocated(cfg%vars)) then
       cur_size = size(cfg%vars)

       if (cur_size - cfg%n_vars < 1) then
          new_size = 2 * cur_size
          allocate(cfg_copy(cur_size))
          cfg_copy = cfg%vars
          deallocate(cfg%vars)
          allocate(cfg%vars(new_size))
          cfg%vars(1:cur_size) = cfg_copy
       end if
    else
       allocate(cfg%vars(min_dyn_size))
    end if

  end subroutine ensure_free_storage

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
    !> On return, startIxs(i) holds the starting point of entry i
    integer, intent(inout)        :: ixs_start(n_max)
    !> On return, endIxs(i) holds the end point of entry i
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

  !> Sort the variables for faster lookup
  subroutine CFG_sort(cfg)
    type(CFG_t), intent(inout) :: cfg
    call qsort_config(cfg%vars(1:cfg%n_vars))
    cfg%sorted = .true.
  end subroutine CFG_sort

  !> Simple implementation of quicksort algorithm to sort the variable list alphabetically.
  recursive subroutine qsort_config(list)
    type(CFG_var_t), intent(inout) :: list(:)
    integer                        :: splitPos

    if (size(list) > 1) then
       call parition_var_list(list, splitPos)
       call qsort_config( list(:splitPos-1) )
       call qsort_config( list(splitPos:) )
    end if
  end subroutine qsort_config

  !> Helper routine for quicksort, to perform partitioning
  subroutine parition_var_list(list, marker)
    type(CFG_var_t), intent(inout) :: list(:)
    integer, intent(out)           :: marker

    integer                        :: left, right, pivot
    type(CFG_var_t)                :: temp
    character(len=name_len)        :: pivotValue

    left       = 0
    right      = size(list) + 1

    ! Take the middle element as pivot
    pivot      = size(list) / 2
    pivotValue = list(pivot)%p_name

    do while (left < right)

       right = right - 1
       do while (lgt(list(right)%p_name, pivotValue))
          right = right - 1
       end do

       left = left + 1
       do while (lgt(pivotValue, list(left)%p_name))
          left = left + 1
       end do

       if (left < right) then
          temp = list(left)
          list(left) = list(right)
          list(right) = temp
       end if
    end do

    if (left == right) then
       marker = left + 1
    else
       marker = left
    end if
  end subroutine parition_var_list

end module m_config
