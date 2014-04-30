! Author: Jannis Teunissen
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

!> Module that allows working with a configuration file.
module m_config

   implicit none
   private

   ! Types
   integer, parameter :: dp = kind(0.0d0)
   integer, parameter :: CFG_int_type = 0, CFG_real_type = 1, CFG_char_type = 2, CFG_logic_type = 3

   ! String lengths
   integer, parameter :: tiny_len = 20
   integer, parameter :: name_len = 40
   integer, parameter :: line_len = 400

   ! Maximum length of a variable
   integer, parameter :: max_var_len = 20

   type CFG_t_var
      character(len=name_len)          :: pName
      character(len=line_len)          :: comment
      integer                          :: pType
      integer                          :: pSize
      logical                          :: varSize

      real(dp), pointer                :: real_data(:)
      integer, pointer                 :: int_data(:)
      character(len=name_len), pointer :: char_data(:)
      logical, pointer                 :: logic_data(:)
   end type CFG_t_var

   logical  :: CFG_sorted = .false.
   integer  :: CFG_num_vars = 0
   integer  :: CFG_max_vars = 0

   type(CFG_t_var), pointer :: CFG_varList(:) => null()

   interface
      subroutine p_err_handler(err_string)
         character(len=*), intent(in) :: err_string
      end subroutine p_err_handler
   end interface

   procedure(p_err_handler), pointer :: CFG_err_handler => null()

   !> Overloaded interface to add config variables
   interface CFG_add
      module procedure CFG_addReal, CFG_addRealArray, CFG_addInt, CFG_addIntArray, &
           CFG_addChar, CFG_addCharArray, CFG_addlogical, CFG_addlogicalArray
   end interface CFG_add

   !> Overloaded interface to get variables
   interface CFG_get
      module procedure CFG_getReal, CFG_getRealArray, CFG_getInt, CFG_getIntArray, &
           CFG_getChar, CFG_getCharArray, CFG_getlogical, CFG_getlogicalArray
   end interface CFG_get

   ! Public types
   public :: CFG_int_type, CFG_real_type, CFG_char_type, CFG_logic_type

   ! Public procedures
   public :: CFG_destruct
   public :: CFG_add, CFG_read_file
   public :: CFG_get_int, CFG_get_real, CFG_get_logic, CFG_get_string
   public :: CFG_get_size, CFG_get_type, CFG_get
   public :: CFG_sort, CFG_write

contains

   subroutine CFG_set_error_handler(err_handler)
      procedure(p_err_handler) :: err_handler
      CFG_err_handler => err_handler
   end subroutine CFG_set_error_handler

   subroutine CFG_handle_error(err_string)
      character(len=*), intent(in) :: err_string

      if (associated(CFG_err_handler)) then
         call CFG_err_handler(trim(err_string))
      else
         print *, "The following fatal error occured in m_config:"
         print *, trim(err_string)
         stop
      end if
   end subroutine CFG_handle_error

   !> Destroy the configuration and deallocate the space used
   subroutine CFG_destruct()
      integer :: i

      do i = 1, CFG_num_vars
         select case (CFG_varList(i)%pType)
         case (CFG_real_type)
            deallocate(CFG_varList(i)%real_data)
         case (CFG_int_type)
            deallocate(CFG_varList(i)%int_data)
         case (CFG_char_type)
            deallocate(CFG_varList(i)%char_data)
         case (CFG_logic_type)
            deallocate(CFG_varList(i)%logic_data)
         end select
      end do

      deallocate(CFG_varList)
   end subroutine CFG_destruct

   !> Update the variables in the configartion with the values found in 'filename'
   subroutine CFG_read_file(filename)
      character(len=*), intent(in) :: filename

      integer                      :: ioState, equalSignPos, lineEnd
      integer                      :: n, ix, nL
      integer                      :: startIxs(max_var_len), endIxs(max_var_len), nEntries
      character(len=tiny_len)      :: configFormat
      character(len=name_len)      :: pName
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

         pName = line(1 : equalSignPos - 1)  ! Set variable name
         pName = adjustl(pName)              ! Remove leading blanks

         line = line(equalSignPos + 1 : lineEnd)   ! Set line to the values behind the '=' sign
         line = adjustl(line)                      ! Remove leading blanks

         ! Find variable corresponding to name in file
         ix = get_var_index(pName)

         if (ix > 0) then

            if (CFG_varList(ix)%varSize) then
               ! Get the start and end positions of the line content, and the number of entries
               call split_line_by_delims(line, " ,'"""//char(9), max_var_len, nEntries, startIxs, endIxs)

               if (CFG_varList(ix)%pSize /= nEntries) then
                  call resize_storage(ix, CFG_varList(ix)%pType, nEntries)
                  CFG_varList(ix)%pSize = nEntries
               end if

               do n = 1, nEntries
                  select case (CFG_varList(ix)%pType)
                  case (CFG_int_type)
                     read(line(startIxs(n):endIxs(n)), *, ERR = 999) CFG_varList(ix)%int_data(n)
                  case (CFG_real_type)
                     read(line(startIxs(n):endIxs(n)), *, ERR = 999) CFG_varList(ix)%real_data(n)
                  case (CFG_char_type)
                     read(line(startIxs(n):endIxs(n)), *, ERR = 999) CFG_varList(ix)%char_data(n)
                  case (CFG_logic_type)
                     read(line(startIxs(n):endIxs(n)), *, ERR = 999) CFG_varList(ix)%logic_data(n)
                  end select
               end do
            else                  ! Fixed size parameter
               select case (CFG_varList(ix)%pType)
               case (CFG_int_type)
                  read(line, *, ERR = 999, end = 998) CFG_varList(ix)%int_data
               case (CFG_real_type)
                  read(line, *, ERR = 999, end = 998) CFG_varList(ix)%real_data
               case (CFG_char_type)
                  read(line, *, ERR = 999, end = 998) CFG_varList(ix)%char_data
               case (CFG_logic_type)
                  read(line, *, ERR = 999, end = 998) CFG_varList(ix)%logic_data
               end select
            end if

         else
            call CFG_handle_error("CFG_read_file: variable [" // trim(pName) // &
                 "] from " //  filename // " could not be updated")
         end if

         ! Go to the next iteration
         cycle
      end do

666   continue ! Routine ends here if the end of "filename" is reached
      close(UNIT = 1, STATUS = "KEEP", ERR = 999, IOSTAT = ioState)
      return

      ! The routine only ends up here through an error
998   continue
call CFG_handle_error("CFG_read_file: Not enough values for variable [" // trim(pName) // "] in " // filename)
      return

999   continue
      write(err_string, *) "ioState = ", ioState, " while reading from ", filename, " at line ", nL
      call CFG_handle_error("CFG_read_file:" // err_string)
      return

   end subroutine CFG_read_file

   subroutine resize_storage(ix, pType, pSize)
      integer, intent(in) :: ix, pType, pSize

      select case (pType)
      case (CFG_int_type)
         deallocate( CFG_varList(ix)%int_data )
         allocate( CFG_varList(ix)%int_data(pSize) )
      case (CFG_logic_type)
         deallocate( CFG_varList(ix)%logic_data )
         allocate( CFG_varList(ix)%logic_data(pSize) )
      case (CFG_real_type)
         deallocate( CFG_varList(ix)%real_data )
         allocate( CFG_varList(ix)%real_data(pSize) )
      case (CFG_char_type)
         deallocate( CFG_varList(ix)%char_data )
         allocate( CFG_varList(ix)%char_data(pSize) )
      end select
   end subroutine resize_storage

   !> Return the index of the variable with name 'pName', or -1 if not found.
   integer function get_var_index(pName)
      character(len=*), intent(in)  :: pName
      integer                       :: i

      if (CFG_sorted) then
         get_var_index = bin_search_var(pName)
      else
         get_var_index = -1
         do i = 1, CFG_num_vars
            if (CFG_varList(i)%pName == pName) then
               get_var_index = i
               exit
            end if
         end do
      end if

   end function get_var_index

   !> Performa a binary search for the variable 'pName'
   integer function bin_search_var(pName)
      character(len=*), intent(in)  :: pName

      integer                       :: iMin, iMax, iMiddle

      iMin           = 1
      iMax           = CFG_num_vars
      bin_search_var = - 1

      do while (iMin <= iMax)
         iMiddle = iMin + (iMax - iMin) / 2

         if ( lgt(CFG_varList(iMiddle)%pName, pName) ) then
            iMax = iMiddle - 1
         else if ( llt(CFG_varList(iMiddle)%pName, pName) ) then
            iMin = iMiddle + 1
         else
            bin_search_var = iMiddle
            exit
         end if
      end do
      ! If not found, bin_search_var is not set here, and stays -1
   end function bin_search_var

   !> Helper routine to create variables. This is useful because a lot of the same
   !! code is executed for the different types of variables.
   subroutine prepareStoreVar(pName, pType, pSize, comment, varSize)
      character(len=*), intent(in)  :: pName, comment
      integer, intent(in)           :: pType, pSize
      logical, intent(in), optional :: varSize

      call ensure_free_storage()

      CFG_sorted                     = .false.
      CFG_num_vars                      = CFG_num_vars + 1
      CFG_varList(CFG_num_vars)%pName   = pName
      CFG_varList(CFG_num_vars)%comment = comment
      CFG_varList(CFG_num_vars)%pType   = pType
      CFG_varList(CFG_num_vars)%pSize   = pSize

      if (present(varSize)) then
         CFG_varList(CFG_num_vars)%varSize = varSize
      else
         CFG_varList(CFG_num_vars)%varSize = .false.
      end if

      select case (CFG_varList(CFG_num_vars)%pType)
      case (CFG_int_type)
         allocate( CFG_varList(CFG_num_vars)%int_data(pSize) )
      case (CFG_real_type)
         allocate( CFG_varList(CFG_num_vars)%real_data(pSize) )
      case (CFG_char_type)
         allocate( CFG_varList(CFG_num_vars)%char_data(pSize) )
      case (CFG_logic_type)
         allocate( CFG_varList(CFG_num_vars)%logic_data(pSize) )
      end select

   end subroutine prepareStoreVar

   integer function prepareGetVar(pName, pType, pSize)
      character(len=*), intent(in)  :: pName
      integer, intent(in)           :: pType, pSize

      integer                       :: ix
      logical                       :: correctType, correctSize

      ix            = get_var_index(pName)
      prepareGetVar = ix
      correctType   = .false.
      correctSize   = .false.

      if (ix /= -1) then
         if (CFG_varList(ix)%pType == pType) correctType = .true.
         if (CFG_varList(ix)%pSize == pSize) correctSize = .true.
      end if

      if (ix == -1) then
         call CFG_handle_error("CFG_get(something): variable ["//trim(pName)//"] not found")
      else if (.not. correctType) then
         call CFG_handle_error("CFG_get(something): variable ["//trim(pName)//"] has different type")
      else if (.not. correctSize) then
         call CFG_handle_error("CFG_get(something): variable["//trim(pName)//"] has different size")
      end if

   end function prepareGetVar

   subroutine CFG_addReal(pName, real_data, comment)
      character(len=*), intent(in)  :: pName, comment
      real(dp), intent(in)  :: real_data
      call prepareStoreVar(pName, CFG_real_type, 1, comment)
      CFG_varList(CFG_num_vars)%real_data(1) = real_data
   end subroutine CFG_addReal

   subroutine CFG_addRealArray(pName, real_data, comment, varSize)
      character(len=*), intent(in)  :: pName, comment
      real(dp), intent(in)  :: real_data(:)
      logical, intent(in), optional :: varSize
      call prepareStoreVar(pName, CFG_real_type, size(real_data), comment, varSize)
      CFG_varList(CFG_num_vars)%real_data = real_data
   end subroutine CFG_addRealArray

   subroutine CFG_addInt(pName, int_data, comment)
      character(len=*), intent(in)  :: pName, comment
      integer, intent(in)  :: int_data
      call prepareStoreVar(pName, CFG_int_type, 1, comment)
      CFG_varList(CFG_num_vars)%int_data(1) = int_data
   end subroutine CFG_addInt

   subroutine CFG_addIntArray(pName, int_data, comment, varSize)
      character(len=*), intent(in)  :: pName, comment
      integer, intent(in)           :: int_data(:)
      logical, intent(in), optional :: varSize
      call prepareStoreVar(pName, CFG_int_type, size(int_data), comment, varSize)
      CFG_varList(CFG_num_vars)%int_data = int_data
   end subroutine CFG_addIntArray

   subroutine CFG_addChar(pName, char_data, comment)
      character(len=*), intent(in)  :: pName, comment, char_data
      call prepareStoreVar(pName, CFG_char_type, 1, comment)
      CFG_varList(CFG_num_vars)%char_data(1) = char_data
   end subroutine CFG_addChar

   subroutine CFG_addCharArray(pName, char_data, comment, varSize)
      character(len=*), intent(in)  :: pName, comment, char_data(:)
      logical, intent(in), optional :: varSize
      call prepareStoreVar(pName, CFG_char_type, size(char_data), comment, varSize)
      CFG_varList(CFG_num_vars)%char_data = char_data
   end subroutine CFG_addCharArray

   subroutine CFG_addlogical(pName, logic_data, comment)
      character(len=*), intent(in)  :: pName, comment
      logical, intent(in)           :: logic_data
      call prepareStoreVar(pName, CFG_logic_type, 1, comment)
      CFG_varList(CFG_num_vars)%logic_data(1) = logic_data
   end subroutine CFG_addlogical

   subroutine CFG_addlogicalArray(pName, logic_data, comment, varSize)
      character(len=*), intent(in)  :: pName, comment
      logical, intent(in)           :: logic_data(:)
      logical, intent(in), optional :: varSize
      call prepareStoreVar(pName, CFG_logic_type, size(logic_data), comment, varSize)
      CFG_varList(CFG_num_vars)%logic_data = logic_data
   end subroutine CFG_addlogicalArray

   subroutine CFG_getReal(pName, real_data)
      character(len=*), intent(in)     :: pName
      real(dp), intent(inout)  :: real_data
      integer :: ix
      ix          = prepareGetVar(pName, CFG_real_type, 1)
      real_data  = CFG_varList(ix)%real_data(1)
   end subroutine CFG_getReal

   subroutine CFG_getRealArray(pName, real_data)
      character(len=*), intent(in)     :: pName
      real(dp), intent(inout)  :: real_data(:)
      integer :: ix
      ix          = prepareGetVar(pName, CFG_real_type, size(real_data))
      real_data  = CFG_varList(ix)%real_data
   end subroutine CFG_getRealArray

   subroutine CFG_getInt(pName, int_data)
      character(len=*), intent(in)     :: pName
      integer, intent(inout)           :: int_data
      integer :: ix
      ix          = prepareGetVar(pName, CFG_int_type, 1)
      int_data    = CFG_varList(ix)%int_data(1)
   end subroutine CFG_getInt

   subroutine CFG_getIntArray(pName, int_data)
      character(len=*), intent(in)     :: pName
      integer, intent(inout)           :: int_data(:)
      integer :: ix
      ix          = prepareGetVar(pName, CFG_int_type, size(int_data))
      int_data    = CFG_varList(ix)%int_data
   end subroutine CFG_getIntArray

   subroutine CFG_getChar(pName, char_data)
      character(len=*), intent(in)     :: pName
      character(len=*), intent(inout)  :: char_data
      integer :: ix
      ix          = prepareGetVar(pName, CFG_char_type, 1)
      char_data   = CFG_varList(ix)%char_data(1)
   end subroutine CFG_getChar

   subroutine CFG_getCharArray(pName, char_data)
      character(len=*), intent(in)     :: pName
      character(len=*), intent(inout)  :: char_data(:)
      integer :: ix
      ix          = prepareGetVar(pName, CFG_char_type, size(char_data))
      char_data   = CFG_varList(ix)%char_data
   end subroutine CFG_getCharArray

   subroutine CFG_getlogical(pName, logic_data)
      character(len=*), intent(in)     :: pName
      logical, intent(inout)           :: logic_data
      integer :: ix
      ix          = prepareGetVar(pName, CFG_logic_type, 1)
      logic_data   = CFG_varList(ix)%logic_data(1)
   end subroutine CFG_getlogical

   subroutine CFG_getlogicalArray(pName, logic_data)
      character(len=*), intent(in)     :: pName
      logical, intent(inout)           :: logic_data(:)
      integer :: ix
      ix          = prepareGetVar(pName, CFG_logic_type, size(logic_data))
      logic_data   = CFG_varList(ix)%logic_data
   end subroutine CFG_getlogicalArray

   subroutine prepare_get_var(pName, pType, argIndex, ix, pIndex)
      character(len=*), intent(in)  :: pName
      integer, intent(in)           :: pType
      integer, intent(inout)        :: ix, pIndex
      integer, intent(in), optional :: argIndex

      if (present(argIndex)) then
         pIndex = argIndex
      else
         pIndex = 1
      end if

      ix = get_var_index(pName)

      if (ix == -1) then
         call CFG_handle_error("CFG_get(something): variable ["//trim(pName)//"] not found")
      else if (CFG_varList(ix)%pType /= pType) then
         call CFG_handle_error("CFG_get(something): variable ["//trim(pName)//"] has different type")
      else if ((.not. present(argIndex) .and. CFG_varList(ix)%pSize > 1) .or. &
           pIndex > CFG_varList(ix)%pSize) then
         call CFG_handle_error("CFG_get(something): variable ["//trim(pName)//"] wrong or no index specified")
      end if

   end subroutine prepare_get_var

   real(dp) function CFG_get_real(pName, argIndex)
      character(len=*), intent(in)  :: pName
      integer, intent(in), optional :: argIndex
      integer :: ix, pIndex
      call prepare_get_var(pName, CFG_real_type, argIndex, ix, pIndex)
      CFG_get_real = CFG_varList(ix)%real_data(pIndex)
   end function CFG_get_real

   integer function CFG_get_int(pName, argIndex)
      character(len=*), intent(in)  :: pName
      integer, intent(in), optional :: argIndex
      integer :: ix, pIndex
      call prepare_get_var(pName, CFG_int_type, argIndex, ix, pIndex)
      CFG_get_int = CFG_varList(ix)%int_data(pIndex)
   end function CFG_get_int

   logical function CFG_get_logic(pName, argIndex)
      character(len=*), intent(in)  :: pName
      integer, intent(in), optional :: argIndex
      integer :: ix, pIndex
      call prepare_get_var(pName, CFG_logic_type, argIndex, ix, pIndex)
      CFG_get_logic = CFG_varList(ix)%logic_data(pIndex)
   end function CFG_get_logic

   character(len=name_len) function CFG_get_string(pName, argIndex)
      character(len=*), intent(in)  :: pName
      integer, intent(in), optional :: argIndex
      integer :: ix, pIndex
      call prepare_get_var(pName, CFG_char_type, argIndex, ix, pIndex)
      CFG_get_string = CFG_varList(ix)%char_data(pIndex)
   end function CFG_get_string

   integer function CFG_get_size(pName)
      character(len=*), intent(in)  :: pName

      integer                       :: ix

      ix = get_var_index(pName)
      if (ix /= -1) then
         CFG_get_size = CFG_varList(ix)%pSize
      else
         CFG_get_size = -1
         call CFG_handle_error("CFG_get_size: variable ["//trim(pName)//"] not found")
      end if
   end function CFG_get_size

   integer function CFG_get_type(pName)
      character(len=*), intent(in)  :: pName

      integer                       :: ix

      ix = get_var_index(pName)
      if (ix /= -1) then
         CFG_get_type = CFG_varList(ix)%pType
      else
         CFG_get_type = -1
         call CFG_handle_error("CFG_get_type: variable ["//trim(pName)//"] not found")
      end if
   end function CFG_get_type

   subroutine CFG_write(filename)
      use iso_fortran_env
      character(len=*), intent(in) :: filename

      integer                      :: i, j, ioState, myUnit
      character(len=tiny_len)      :: nameFormat
      character(len=line_len)      :: err_string

      write(nameFormat, FMT = "(I6)") name_len
      nameFormat = "(A,A" // trim(adjustl(nameFormat)) // ",A)"

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

      do i = 1, CFG_num_vars
         write(myUnit, ERR = 999, FMT = "(A,A,A)") " # ", trim(CFG_varList(i)%comment), ":"
         write(myUnit, ADVANCE = "NO", ERR = 999, FMT = nameFormat) " ", CFG_varList(i)%pName, " = "

         select case(CFG_varList(i)%pType)
         case (CFG_int_type)
            do j = 1, CFG_varList(i)%pSize
               write(myUnit, ADVANCE = "NO", ERR = 999, FMT = "(I10, A) ") CFG_varList(i)%int_data(j), "  "
            end do
         case (CFG_real_type)
            do j = 1, CFG_varList(i)%pSize
               write(myUnit, ADVANCE = "NO", ERR = 999, FMT = "(E10.4, A) ") CFG_varList(i)%real_data(j), "  "
            end do
         case (CFG_char_type)
            do j = 1, CFG_varList(i)%pSize
               write(myUnit, ADVANCE = "NO", ERR = 999, FMT = "(A, A)") &
                    & '"' // trim(CFG_varList(i)%char_data(j)) // '"', "  "
            end do
         case (CFG_logic_type)
            do j = 1, CFG_varList(i)%pSize
               write(myUnit, ADVANCE = "NO", ERR = 999, FMT = "(L10, A) ") CFG_varList(i)%logic_data(j), "  "
            end do
         end select
         write(myUnit, ERR = 999, FMT = "(A)") ""
         write(myUnit, ERR = 999, FMT = "(A)") ""
      end do

      if (myUnit /= output_unit) close(myUnit, STATUS = "KEEP", ERR = 999, IOSTAT = ioState)

      return

999   continue ! If there was an error, the routine will end here
      write(err_string, *) "CFG_write error: ioState = ", ioState, " while writing to ", filename
      call CFG_handle_error(err_string)

   end subroutine CFG_write

   subroutine ensure_free_storage()
      integer, parameter         :: min_var_size = 100
      integer                    :: cur_size, new_size
      type(CFG_t_var), pointer :: cfg_copy(:)

      if (associated(CFG_varList)) then
         cur_size = size(CFG_varList)
         if (cur_size - CFG_num_vars < 1) then
            new_size = 2 * cur_size
            allocate(cfg_copy(new_size))
            cfg_copy(1:CFG_num_vars) = CFG_varList
            deallocate(CFG_varList)
            CFG_varList => cfg_copy
         end if
      else
         allocate(CFG_varList(min_var_size))
      end if

   end subroutine ensure_free_storage

   !> Routine to help read a variable number of entries from a line.
   ! In arguments:
   !  line           The line from which we want to read
   !  delims         A string with delimiters. For example delims = " ,'"""//char(9)
   !                 " ,'"""//char(9) = space, comma, single/real quotation marks, tab
   !  nEntriesMax    Maximum number of entries to read in
   !
   ! Out arguments:
   !  nEntries       Number of entries found
   !  startIxs       On return, startIxs(i) holds the starting point of entry i
   !  endIxs         On return, endIxs(i) holds the end point of entry i
   subroutine split_line_by_delims(line, delims, nEntriesMax, nEntries, startIxs, endIxs)
      character(len=*), intent(in)  :: line, delims
      integer, intent(in)           :: nEntriesMax
      integer, intent(inout)        :: nEntries, startIxs(nEntriesMax), endIxs(nEntriesMax)

      integer                       :: ix, prevIx

      prevIx   = 0
      nEntries = 0

      do while (nEntries < nEntriesMax)

         ! Find the starting point of the next entry (a non-delimiter value)
         ix                   = verify(line(prevIx+1:), delims)
         if (ix == 0) exit                      ! No more entries

         nEntries             = nEntries + 1
         startIxs(nEntries)   = prevIx + ix     ! This is the absolute position in 'line'

         ! Get the end point of the current entry (next delimiter index minus one)
         ix = scan(line(startIxs(nEntries)+1:), delims) - 1

         if (ix == -1) then                     ! If there is no last delimiter,
            endIxs(nEntries)  = len(line)       ! the end of the line is the endpoint
         else
            endIxs(nEntries)  = startIxs(nEntries) + ix
         end if

         prevIx = endIxs(nEntries)              ! We continue to search from here
      end do

   end subroutine split_line_by_delims

   !> Sort the variables for faster lookup
   subroutine CFG_sort()
      call qsort_config(CFG_varList(1:CFG_num_vars))
      CFG_sorted = .true.
   end subroutine CFG_sort

   !> Simple implementation of quicksort algorithm to sort the variable list alphabetically.
   recursive subroutine qsort_config(list)
      type(CFG_t_var), intent(inout) :: list(:)
      integer                          :: splitPos

      if (size(list) > 1) then
         call parition_var_list(list, splitPos)
         call qsort_config( list(:splitPos-1) )
         call qsort_config( list(splitPos:) )
      end if
   end subroutine qsort_config

   !> Helper routine for quicksort, to perform partitioning
   subroutine parition_var_list(list, marker)
      type(CFG_t_var), intent(inout) :: list(:)
      integer, intent(out)             :: marker

      integer                 :: left, right, pivot
      type(CFG_t_var)       :: temp
      character(len=name_len)  :: pivotValue

      left       = 0
      right      = size(list) + 1

      ! Take the middle element as pivot
      pivot      = size(list) / 2
      pivotValue = list(pivot)%pName

      do while (left < right)

         right = right - 1
         do while (lgt(list(right)%pName, pivotValue))
            right = right - 1
         end do

         left = left + 1
         do while (lgt(pivotValue, list(left)%pName))
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
