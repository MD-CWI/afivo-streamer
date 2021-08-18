module m_npy
  use iso_fortran_env

  implicit none
  private

  ! Suffix for temporary .npy files
  character(len=*), parameter :: npy_suffix = '.npy'

  interface save_npy
    module procedure write_int64_vec, write_int64_mtx, &
         write_int32_vec, write_int32_mtx, write_int32_3d, &
         write_int16_vec, write_int16_mtx, &
         write_int8_vec, write_int8_mtx, write_int8_3d, &
         write_dbl_vec, write_dbl_mtx, &
         write_sng_vec, write_sng_mtx, &
         write_cmplx_sgn_vec, write_cmplx_sgn_mtx, &
         write_cmplx_dbl_vec, write_cmplx_dbl_mtx, &
         write_sng_3dT, write_dbl_3dT, &
         write_sng_4dT, write_dbl_4dT, &
         write_dbl_5dT, &
         write_cmplx_dbl_3dT, &
         write_cmplx_dbl_4dT, &
         write_cmplx_dbl_5dT, &
         write_cmplx_dbl_6dT
  end interface save_npy

  public :: save_npy
  public :: remove_file
  public :: add_to_zip

contains

  subroutine run_sys(cmd, stat)
    character(len=*), intent(in) :: cmd
    integer(int32), intent(out)  :: stat

    call execute_command_line(cmd, wait=.true., exitstat=stat)
  end subroutine run_sys

  ! Add npy file to a zip file and remove it
  subroutine add_to_zip(zipfile, filename, keep_file, custom_name)
    character(len=*), intent(in)           :: zipfile     ! Name of zip file
    character(len=*), intent(in)           :: filename    ! Name of file to add
    logical, intent(in)                    :: keep_file   ! Whether to keep 'filename'
    character(len=*), intent(in), optional :: custom_name ! Custom name
    integer(int32)                         :: stat

    ! Be quiet while zipping
    character(len=*), parameter  :: zip_command = "zip -q0"

    call run_sys(zip_command//" "//trim(zipfile)//" "//&
         trim(filename), stat)
    if (stat /= 0) then
       print *, zip_command//" "//trim(zipfile)//" "// trim(filename)
       error stop "add_to_zip: Can't execute zip command"
    endif

    if (present(custom_name)) then
       call run_sys('printf "@ '//trim(filename)//'\n@='//&
            trim(custom_name)//'\n" | zipnote -w '//trim(zipfile), stat)
       if (stat /= 0) then
          error stop "add_to_zip: Failed to rename to custom_name"
       endif
    end if

    if (.not. keep_file) then
       call remove_file(filename)
    end if
  end subroutine add_to_zip

  subroutine remove_file(filename)
    character(len=*), intent(in) :: filename
    integer                      :: p_un, stat

    open(newunit=p_un, iostat=stat, file=filename, status='old')
    if (stat == 0) close(p_un, status='delete')
  end subroutine remove_file

  function dict_str(var_type, var_shape) result(str)
    character(len=*), intent(in)  :: var_type
    integer(int32), intent(in)    :: var_shape(:)
    character(len=:), allocatable :: str
    character(len=1024)           :: buffer
    integer(int32)                :: total_size, my_size

    ! https://numpy.org/devdocs/reference/generated/numpy.lib.format.html

    ! The first 6 bytes are a magic string: exactly \x93NUMPY.

    ! The next 1 byte is an unsigned byte: the major version number of the file
    ! format, e.g. \x01.

    ! The next 1 byte is an unsigned byte: the minor version number of the file
    ! format, e.g. \x00. Note: the version of the file format is not tied to the
    ! version of the numpy package.

    ! The next 2 bytes form a little-endian unsigned short int: the length of
    ! the header data HEADER_LEN.

    ! The next HEADER_LEN bytes form the header data describing the array’s
    ! format. It is an ASCII string which contains a Python literal expression
    ! of a dictionary. It is terminated by a newline (\n) and padded with spaces
    ! (\x20) to make the total of len(magic string) + 2 + len(length) +
    ! HEADER_LEN be evenly divisible by 64 for alignment purposes.

    buffer = "{'descr': '"//var_type// &
         "', 'fortran_order': True, 'shape': ("// &
         shape_str(var_shape)//"), }"

    ! len(magic string) + 2 + len(length) + ending newline =
    ! 6 + 2 + 4 + 1 = 13 bytes
    total_size = len_trim(buffer) + 13

    ! ensure total_size is divisible by 16 bytes
    total_size = ((total_size + 15)/16) * 16

    ! Size of dict_str includes the ending newline (so -12 instead of -13)
    my_size = total_size - 12

    ! End with newline
    buffer(my_size:my_size) = achar(10)
    str = buffer(1:my_size)
  end function dict_str

  function shape_str(var_shape) result(fin_str)
    integer(int32), intent(in)    :: var_shape(:)
    character(len=:), allocatable :: str, small_str, fin_str
    integer(int32)                :: i, length, start, halt

    length = 14*size(var_shape)
    allocate (character(length) :: str)
    allocate (character(14)     :: small_str)
    str = " "

    do i = 1, size(var_shape)
      start = (i - 1)*length + 1
      halt = i*length + 1
      write (small_str, "(I13,A)") var_shape(i), ","
      str = trim(str)//adjustl(small_str)
    enddo

    fin_str = trim(str)
  end function shape_str

  subroutine write_header(p_un, var_type, var_shape)
    integer(int32), intent(in)   :: p_un
    character(len=*), intent(in) :: var_type
    integer(int32), intent(in)   :: var_shape(:)
    integer(int32)               :: header_len

    ! Magic number hex x93 is 147 (unsigned), signed this is -109
    integer(int8), parameter    :: magic_num = int(-109, int8)
    character(len=*), parameter :: magic_str = "NUMPY"
    integer(int8), parameter    :: major = 2_int8 ! major *.npy version
    integer(int8), parameter    :: minor = 0_int8 ! minor *.npy version

    header_len = len(dict_str(var_type, var_shape))
    write (p_un) magic_num, magic_str, major, minor
    write (p_un) header_len
    write (p_un) dict_str(var_type, var_shape)
  end subroutine write_header

  subroutine write_cmplx_sgn_mtx(filename, mtx)
    character(len=*), intent(in) :: filename
    complex(4), intent(in)       :: mtx(:, :)
    character(len=*), parameter  :: var_type = "<c8"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(mtx))
    write (p_un) mtx
    close (unit=p_un)
  end subroutine write_cmplx_sgn_mtx

  subroutine write_cmplx_sgn_vec(filename, vec)
    character(len=*), intent(in) :: filename
    complex(4), intent(in)       :: vec(:)
    character(len=*), parameter  :: var_type = "<c8"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(vec))
    write (p_un) vec
    close (unit=p_un)
  end subroutine write_cmplx_sgn_vec

  subroutine write_cmplx_dbl_6dT(filename, tensor)
    character(len=*), intent(in) :: filename
    complex(8), intent(in)       :: tensor(:, :, :, :, :, :)
    character(len=*), parameter  :: var_type = "<c16"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(tensor))
    write (p_un) tensor
    close (unit=p_un)
  end subroutine write_cmplx_dbl_6dT

  subroutine write_cmplx_dbl_5dT(filename, tensor)
    character(len=*), intent(in) :: filename
    complex(8), intent(in)       :: tensor(:, :, :, :, :)
    character(len=*), parameter  :: var_type = "<c16"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(tensor))
    write (p_un) tensor
    close (unit=p_un)
  end subroutine write_cmplx_dbl_5dT

  subroutine write_cmplx_dbl_4dT(filename, tensor)
    character(len=*), intent(in) :: filename
    complex(8), intent(in)       :: tensor(:, :, :, :)
    character(len=*), parameter  :: var_type = "<c16"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(tensor))
    write (p_un) tensor
    close (unit=p_un)
  end subroutine write_cmplx_dbl_4dT

  subroutine write_cmplx_dbl_3dT(filename, tensor)
    character(len=*), intent(in) :: filename
    complex(8), intent(in)       :: tensor(:, :, :)
    character(len=*), parameter  :: var_type = "<c16"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(tensor))
    write (p_un) tensor
    close (unit=p_un)
  end subroutine write_cmplx_dbl_3dT

  subroutine write_cmplx_dbl_mtx(filename, mtx)
    character(len=*), intent(in) :: filename
    complex(8), intent(in)       :: mtx(:, :)
    character(len=*), parameter  :: var_type = "<c16"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(mtx))
    write (p_un) mtx
    close (unit=p_un)
  end subroutine write_cmplx_dbl_mtx

  subroutine write_cmplx_dbl_vec(filename, vec)
    character(len=*), intent(in) :: filename
    complex(8), intent(in)       :: vec(:)
    character(len=*), parameter  :: var_type = "<c16"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(vec))
    write (p_un) vec
    close (unit=p_un)
  end subroutine write_cmplx_dbl_vec

  subroutine write_sng_3dT(filename, tensor)
    character(len=*), intent(in) :: filename
    real(real32), intent(in)     :: tensor(:, :, :)
    character(len=*), parameter  :: var_type = "<f4"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(tensor))
    write (p_un) tensor
    close (unit=p_un)
  end subroutine write_sng_3dT

  subroutine write_sng_4dT(filename, tensor)
    character(len=*), intent(in) :: filename
    real(real32), intent(in)     :: tensor(:, :, :, :)
    character(len=*), parameter  :: var_type = "<f4"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(tensor))
    write (p_un) tensor
    close (unit=p_un)
  end subroutine write_sng_4dT

  subroutine write_sng_mtx(filename, mtx)
    character(len=*), intent(in) :: filename
    real(real32), intent(in)     :: mtx(:, :)
    character(len=*), parameter  :: var_type = "<f4"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(mtx))
    write (p_un) mtx
    close (unit=p_un)
  end subroutine write_sng_mtx

  subroutine write_sng_vec(filename, vec)
    character(len=*), intent(in) :: filename
    real(real32), intent(in)          :: vec(:)
    character(len=*), parameter  :: var_type = "<f4"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(vec))
    write (p_un) vec
    close (unit=p_un)
  end subroutine write_sng_vec

  subroutine write_dbl_3dT(filename, tensor)
    character(len=*), intent(in) :: filename
    real(real64), intent(in)     :: tensor(:, :, :)
    character(len=*), parameter  :: var_type = "<f8"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(tensor))
    write (p_un) tensor
    close (unit=p_un)
  end subroutine write_dbl_3dT

  subroutine write_dbl_4dT(filename, tensor4)
    character(len=*), intent(in) :: filename
    real(real64), intent(in)     :: tensor4(:, :, :, :)
    character(len=*), parameter  :: var_type = "<f8"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(tensor4))
    write (p_un) tensor4
    close (unit=p_un)
  end subroutine write_dbl_4dT

  subroutine write_dbl_5dT(filename, tensor5)
    character(len=*), intent(in) :: filename
    real(real64), intent(in)     :: tensor5(:, :, :, :, :)
    character(len=*), parameter  :: var_type = "<f8"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(tensor5))
    write (p_un) tensor5
    close (unit=p_un)
  end subroutine write_dbl_5dT

  subroutine write_dbl_mtx(filename, mtx)
    character(len=*), intent(in) :: filename
    real(real64), intent(in)     :: mtx(:, :)
    character(len=*), parameter  :: var_type = "<f8"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(mtx))
    write (p_un) mtx
    close (unit=p_un)
  end subroutine write_dbl_mtx

  subroutine write_dbl_vec(filename, vec)
    character(len=*), intent(in) :: filename
    real(real64), intent(in)     :: vec(:)
    character(len=*), parameter  :: var_type = "<f8"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(vec))
    write (p_un) vec
    close (unit=p_un)
  end subroutine write_dbl_vec

  subroutine write_int64_mtx(filename, mtx)
    character(len=*), intent(in) :: filename
    integer(int64), intent(in)   :: mtx(:, :)
    character(len=*), parameter  :: var_type = "<i8"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(mtx))
    write (p_un) mtx
    close (unit=p_un)
  end subroutine write_int64_mtx

  subroutine write_int64_vec(filename, vec)
    character(len=*), intent(in) :: filename
    integer(int64), intent(in)   :: vec(:)
    character(len=*), parameter  :: var_type = "<i8"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(vec))
    write (p_un) vec
    close (unit=p_un)
  end subroutine write_int64_vec

  subroutine write_int32_mtx(filename, mtx)
    character(len=*), intent(in) :: filename
    integer(int32), intent(in)   :: mtx(:, :)
    character(len=*), parameter  :: var_type = "<i4"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(mtx))
    write (p_un) mtx
    close (unit=p_un)
  end subroutine write_int32_mtx

  subroutine write_int32_3d(filename, mtx)
    character(len=*), intent(in) :: filename
    integer(int32), intent(in)   :: mtx(:,:,:)
    character(len=*), parameter  :: var_type = "<i4"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(mtx))
    write (p_un) mtx
    close (unit=p_un)
  end subroutine write_int32_3d

  subroutine write_int32_vec(filename, vec)
    character(len=*), intent(in) :: filename
    integer(int32), intent(in)   :: vec(:)
    character(len=*), parameter  :: var_type = "<i4"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(vec))
    write (p_un) vec
    close (unit=p_un)
  end subroutine write_int32_vec

  subroutine write_int16_mtx(filename, mtx)
    character(len=*), intent(in) :: filename
    integer(int16), intent(in)   :: mtx(:, :)
    character(len=*), parameter  :: var_type = "<i2"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(mtx))
    write (p_un) mtx
    close (unit=p_un)
  end subroutine write_int16_mtx

  subroutine write_int16_vec(filename, vec)
    character(len=*), intent(in) :: filename
    integer(int16), intent(in)   :: vec(:)
    character(len=*), parameter  :: var_type = "<i2"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(vec))
    write (p_un) vec
    close (unit=p_un)
  end subroutine write_int16_vec

  subroutine write_int8_mtx(filename, mtx)
    character(len=*), intent(in) :: filename
    integer(int8), intent(in)    :: mtx(:, :)
    character(len=*), parameter  :: var_type = "<i1"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(mtx))
    write (p_un) mtx
    close (unit=p_un)
  end subroutine write_int8_mtx

  subroutine write_int8_3d(filename, mtx)
    character(len=*), intent(in) :: filename
    integer(int8), intent(in)    :: mtx(:,:,:)
    character(len=*), parameter  :: var_type = "<i1"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(mtx))
    write (p_un) mtx
    close (unit=p_un)
  end subroutine write_int8_3d

  subroutine write_int8_vec(filename, vec)
    character(len=*), intent(in) :: filename
    integer(int8), intent(in)    :: vec(:)
    character(len=*), parameter  :: var_type = "<i1"
    integer(int32)               :: p_un
    open (newunit=p_un, file=filename, form="unformatted", access="stream")
    call write_header(p_un, var_type, shape(vec))
    write (p_un) vec
    close (unit=p_un)
  end subroutine write_int8_vec

end module m_npy
