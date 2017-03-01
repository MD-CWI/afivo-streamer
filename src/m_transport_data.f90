!> Module that provides routines for reading in arbritrary transport data

module m_transport_data

   implicit none
   private

   integer, parameter :: dp = kind(0.0d0)

   ! The maximum number of rows per entry
   integer, parameter :: max_num_rows   = 500

   ! The maximum number columns in a 2D transport data table
   integer, parameter :: max_num_cols   = 50

   integer, parameter :: lineLen        = 200

   public :: TD_get_from_file
   public :: TD_get_from_file_2D

contains

   !> Routine to read in transport data from a file
   !! Searches 'file_name' for transport 'data_name' concerning 'gas_name'
   subroutine TD_get_from_file(file_name, gas_name, data_name, x_data, y_data)
      character(len=*), intent(in)       :: file_name, gas_name, data_name
      real(dp), allocatable, intent(out) :: x_data(:), y_data(:)

      ! Temporary variables
      integer                :: gas_name_len
      integer                :: ioState, nL
      integer                :: n_rows
      integer, parameter     :: my_unit = 333
      character(LEN=40)      :: line_fmt
      character(LEN=lineLen) :: line, prev_line
      real(dp)               :: temp_table(2, max_num_rows)

      nL           = 0 ! Set the number of lines to 0
      gas_name_len = len_trim(gas_name)

      ! Set the line format to read, only depends on lineLen currently
      write(line_fmt, FMT = "(I6)") lineLen
      line_fmt = "(A" // trim(adjustl(line_fmt)) // ")"

      ! Open 'file_name' (with error checking)
      open(my_unit, file = trim(file_name), action = "read", &
           err = 999, iostat = ioState, status="old")

      ! Look for collision processes with the correct gas name in the file,
      ! which should contains entries like below:

      !     Efield[V/m]_vs_energy[eV]     [description of the type of transport data]
      !     AIR                           [the gas (mixture) name]
      !     COMMENT: Compiled by xxxxx    [possibly comments, these are optional]
      !     UPDATED: 2010-06-24 15:04:36
      !     ------------------            [at least 5 dashes]
      !     xxx       xxx                 [transport data in two column format]
      !     ...       ...
      !     xxx       xxx
      !     ------------------

      ! The outer DO loop, running until the end of the file is reached
      do
         ! Search for 'gas_name' in the file
         line = ' '
         do
            prev_line = line
            read(my_unit, FMT = line_fmt, ERR = 999, end = 888) line; nL = nL+1
            line = adjustl(line)
            if ( line(1:gas_name_len) == trim(gas_name) ) exit
         end do

         ! prev_line holds the type of transport data, see the formatting above
         ! cycle if this is not the right entry
         if (trim(data_name) /= trim(prev_line)) cycle

         ! Now we can check whether there is a comment, while scanning lines
         ! until dashes are found, which indicate the start of the transport
         ! data
         do
            read(my_unit, FMT = line_fmt, ERR = 999, end = 777) line; nL = nL+1
            line = adjustl(line)
            if ( line(1:8) == "COMMENT:" ) cycle ! Do nothing for now
            if ( line(1:5) == "-----" ) exit
         end do

         ! Read the transport data into a temporary array
         n_rows = 0
         do
            read(my_unit, FMT = line_fmt, ERR = 999, end = 777) line; nL = nL+1
            line = adjustl(line)
            if ( line(1:5) == "-----" ) then
               exit  ! Dashes mark the end of the data
            else if (trim(line) == "" .or. line(1:1) == "#") then
               cycle ! Ignore whitespace or comments
            else if (n_rows < max_num_rows) then
               n_rows = n_rows + 1
               read(line, FMT = *, ERR = 999, end = 777) temp_table(:, n_rows)
            else
               print *, "CS_read_file error: too many rows in ", &
                        file_name, " at line ", nL
            end if
         end do

         ! Store the data in the actual table
         if (allocated(x_data)) deallocate(x_data)
         if (allocated(y_data)) deallocate(y_data)
         allocate(x_data(n_rows))
         allocate(y_data(n_rows))

         x_data = temp_table(1, 1:n_rows)
         y_data = temp_table(2, 1:n_rows)

         exit                   ! Done
      end do

      close(my_unit)
      return

777   continue ! If the end of the file is reached after finding data
      print *, "TD_get_td_from_file unexpectedly reached end of " // trim(file_name)
      print *, "searching '" // trim(gas_name) // "': '" // trim(data_name) // "'"
      stop

888   continue ! If the end of the file is reached without finding data
      print *, "TD_get_td_from_file: no data in " // trim(file_name)
      print *, "searching '" // trim(gas_name) // "': '" // trim(data_name) // "'"
      stop

999   continue ! If there was an input error, the routine will end here
      print *, "TD_get_td_from_file error at line", nL
      print *, "ioState = ", ioState, " in ", file_name
      print *, "searching '" // trim(gas_name) // "': '" // trim(data_name) // "'"
      stop

    end subroutine TD_get_from_file


    !> Routine to read in two-dimensional transport data from a file
    !! Searches 'file_name' for transport 'data_name' concerning 'gas_name'
    subroutine TD_get_from_file_2D(file_name, gas_name, data_name, &
         x1_data, x2_data, y_data)
      character(len=*), intent(in)       :: file_name, gas_name, data_name
      real(dp), allocatable, intent(out) :: x1_data(:), x2_data(:), y_data(:,:)

      ! Temporary variables
      integer                :: gas_name_len
      integer                :: ioState, nL
      integer                :: n_rows, n_cols
      integer, parameter     :: my_unit = 333
      character(LEN=40)      :: line_fmt
      character(LEN=lineLen) :: line, prev_line
      real(dp)               :: temp_table(max_num_cols, max_num_rows)

      nL           = 0 ! Set the number of lines to 0
      gas_name_len = len_trim(gas_name)

      ! Set the line format to read, only depends on lineLen currently
      write(line_fmt, FMT = "(I6)") lineLen
      line_fmt = "(A" // trim(adjustl(line_fmt)) // ")"

      ! Open 'file_name' (with error checking)
      open(my_unit, file = trim(file_name), action = "read", &
           err = 999, iostat = ioState, status="old")

      ! Look for collision processes with the correct gas name in the file,
      ! which should contains entries like below:

      !     Efield[V/m]_vs_energy[eV]     [description of the type of transport data]
      !     AIR                           [the gas (mixture) name]
      !     COMMENT: Compiled by xxxxx    [possibly comments, these are optional]
      !     UPDATED: 2010-06-24 15:04:36
      !     ------------------            [at least 5 dashes]
      !     xxx       xxx                 [transport data in two column format]
      !     ...       ...
      !     xxx       xxx
      !     ------------------

      ! The outer DO loop, running until the end of the file is reached
      do
         ! Search for 'gas_name' in the file
         line = ' '
         do
            prev_line = line
            read(my_unit, FMT = line_fmt, ERR = 999, end = 888) line; nL = nL+1
            line = adjustl(line)
            if ( line(1:gas_name_len) == trim(gas_name) ) exit
         end do

         ! prev_line holds the type of transport data, see the formatting above
         ! cycle if this is not the right entry
         if (trim(data_name) /= trim(prev_line)) cycle

         ! Now we can check whether there is a comment, while scanning lines
         ! until dashes are found, which indicate the start of the transport
         ! data
         do
            read(my_unit, FMT = line_fmt, ERR = 999, end = 777) line; nL = nL+1
            line = adjustl(line)
            if ( line(1:8) == "COMMENT:" ) then
               cycle ! Do nothing for now
            else if ( line(1:12) == "NUM_COLUMNS:" ) then
               read(line(13:), *) n_cols
               if (allocated(x2_data)) deallocate(x2_data)
               allocate(x2_data(n_cols))
            else if ( line(1:11) == "COL_VALUES:" ) then
               read(line(12:), *) x2_data
            else if ( line(1:5) == "-----" ) then
               exit
            end if
         end do

         ! Read the transport data into a temporary array
         n_rows = 0
         do
            read(my_unit, FMT = line_fmt, ERR = 999, end = 777) line; nL = nL+1
            line = adjustl(line)
            if ( line(1:5) == "-----" ) then
               exit  ! Dashes mark the end of the data
            else if (trim(line) == "" .or. line(1:1) == "#") then
               cycle ! Ignore whitespace or comments
            else if (n_rows < max_num_rows) then
               n_rows = n_rows + 1
               read(line, *, ERR = 999, end = 777) &
                    temp_table(1:n_cols+1, n_rows)
            else
               print *, "CS_read_file error: too many rows in ", &
                    file_name, " at line ", nL
            end if
         end do

         ! Store the data in the actual table
         if (allocated(x1_data)) deallocate(x1_data)
         allocate(x1_data(n_rows))

         if (allocated(y_data)) deallocate(y_data)
         allocate(y_data(n_cols, n_rows))

         x1_data = temp_table(1, 1:n_rows)
         y_data = temp_table(2:n_cols+1, 1:n_rows)

         exit                   ! Done
      end do

      close(my_unit)
      return

777   continue ! If the end of the file is reached after finding data
      print *, "TD_get_td_from_file unexpectedly reached end of " // trim(file_name)
      print *, "searching '" // trim(gas_name) // "': '" // trim(data_name) // "'"
      stop

888   continue ! If the end of the file is reached without finding data
      print *, "TD_get_td_from_file: no data in " // trim(file_name)
      print *, "searching '" // trim(gas_name) // "': '" // trim(data_name) // "'"
      stop

999   continue ! If there was an input error, the routine will end here
      print *, "TD_get_td_from_file error at line", nL
      print *, "ioState = ", ioState, " in ", file_name
      print *, "searching '" // trim(gas_name) // "': '" // trim(data_name) // "'"
      stop

    end subroutine TD_get_from_file_2D

end module m_transport_data
