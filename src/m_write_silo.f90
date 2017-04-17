!> This module contains wrapper functions to simplify writing Silo files
!> \todo Document the functions in this module
module m_write_silo

  implicit none
  private

  include 'silo_f9x.inc'

  integer, parameter :: dp       = kind(0.0d0)
  integer, parameter :: line_len = 200
  integer, parameter :: DB_TYPE  = DB_PDB

  public :: SILO_mkdir
  public :: SILO_create_file
  public :: SILO_open_file
  public :: SILO_close_file
  public :: SILO_set_time_varying
  public :: SILO_add_grid
  public :: SILO_add_var
  public :: SILO_set_mmesh_grid
  public :: SILO_set_mmesh_var

contains

  subroutine SILO_create_file(filename, dbix)
    character(len=*), intent(in) :: filename
    integer, intent(out)         :: dbix
    integer                      :: ierr
    character(len=line_len)      :: fileinfo

    fileinfo = "A silo file"
    ierr = dbcreate(trim(filename), len_trim(filename), DB_CLOBBER, &
         DB_LOCAL, fileinfo, len_trim(fileinfo), DB_TYPE, dbix)
    if (ierr /= 0) then
       print *, "Error creating file", trim(filename)
       error stop
    end if
  end subroutine SILO_create_file

  subroutine SILO_open_file(filename, dbix)
    character(len=*), intent(in) :: filename
    integer :: dbix, ierr

    ierr = dbopen(trim(filename), len_trim(filename), DB_TYPE, &
         DB_APPEND, dbix)
    if (ierr /= 0) then
       print *, "Error opening file", trim(filename)
       error stop
    end if
  end subroutine SILO_open_file

  subroutine SILO_close_file(dbix)
    integer, intent(in) :: dbix
    integer :: ierr

    ierr = dbclose(dbix)
    if (ierr /= 0) then
       print *, "Error closing file with index", dbix
       error stop
    end if
  end subroutine SILO_close_file

  subroutine SILO_mkdir(dbix, dirname)
    character(len=*), intent(in) :: dirname
    integer, intent(in) :: dbix
    integer :: ierr, iostat

    ierr = dbmkdir(dbix, trim(dirname), len_trim(dirname), iostat)
    if (ierr /= 0) then
       print *, "Error creating directory ", dirname
       error stop
    end if
  end subroutine SILO_mkdir

  !> Write two entries to the Silo file so that Visit treats it as a
  !> time-varying database
  subroutine SILO_set_time_varying(dbix)
    integer, intent(in)         :: dbix
    integer                     :: ierr
    integer, parameter          :: int_bool(1) = 1
    integer, parameter          :: dims(1) = 1
    character(len=*), parameter :: name1       = "/ConnectivityIsTimeVarying"
    character(len=*), parameter :: name2       = "/MetadataIsTimeVarying"

    interface
       integer function dbwrite(dbid, varname, lvarname, var, dims, &
            ndims, datatype)
         use, intrinsic :: iso_c_binding
         integer(c_int) :: dbid, var(*), lvarname, dims(*), ndims, datatype
         character(kind=c_char) :: varname(*)
       end function dbwrite
    end interface

    ierr = DBWrite(dbix, name1, len(name1), int_bool, dims, 1, DB_INT);
    if (ierr /= 0) print *, "Error writing ", trim(name1)
    ierr = DBWrite(dbix, name2, len(name2), int_bool, dims, 1, DB_INT);
    if (ierr /= 0) print *, "Error writing ", trim(name2)
  end subroutine SILO_set_time_varying

  subroutine SILO_add_grid(dbix, gridname, n_dim, N_r, r_min, dr, &
       lo_offset, hi_offset)
    character(len=*), intent(in) :: gridname
    integer, intent(in)          :: dbix, n_dim, N_r(:)
    integer, intent(in)          :: lo_offset(n_dim), hi_offset(n_dim)
    real(dp), intent(in)         :: r_min(:), dr(:)
    real(dp), allocatable        :: x_coords(:), y_coords(:), z_coords(:)
    integer                      :: i, ierr, iostat, dboptix

    interface
       function dbputqm(dbid, name, lname, xname, lxname, yname, &
            lyname, zname, lzname, x, y, z, dims, ndims, &
            datatype, coordtype, optlist_id, status)
         use, intrinsic :: iso_c_binding
         integer(c_int) :: dbid, lname, lxname, lyname, lzname, dims(*), ndims
         integer(c_int) :: datatype, coordtype, optlist_id, status, dbputqm
         real(c_double) :: x(*), y(*), z(*)
         character(kind=c_char) :: name(*), xname(*), yname(*), zname(*)
       end function dbputqm

       integer (c_int) function dbaddiopt(optlist_id, option, ivalue)
         use, intrinsic :: iso_c_binding
         integer(c_int), intent(in) :: optlist_id, option, ivalue(*)
       end function dbaddiopt
    end interface

    if (n_dim < 1 .or. n_dim > 3) then
       print *, "Cannot add grid for which n_dim < 1 or n_dim > 3"
       return
    end if

    allocate(x_coords(N_r(1)))
    do i = 1, N_r(1)
       x_coords(i) = r_min(1) + (i-1) * dr(1)
    end do

    if (n_dim > 1) then
       allocate(y_coords(N_r(2)))
       do i = 1, N_r(2)
          y_coords(i) = r_min(2) + (i-1) * dr(2)
       end do
    else
       allocate(y_coords(0))
    end if

    if (n_dim > 2) then
       allocate(z_coords(N_r(3)))
       do i = 1, N_r(3)
          z_coords(i) = r_min(3) + (i-1) * dr(3)
       end do
    else
       allocate(z_coords(0))
    end if

    ! Make option list
    ierr = dbmkoptlist(20, dboptix)
    if (ierr /= 0) print *, &
         "Error creating options list in SILO_add_grid ", dboptix

    ! Set integer options
    ierr = dbaddiopt(dboptix, DBOPT_NSPACE, [n_dim])
    if (ierr /= 0) print *, &
            "Error dbaddiopt in SILO_add_grid: DBOPT_NSPACE", ierr

    ierr = dbaddiopt(dboptix, DBOPT_LO_OFFSET, lo_offset)
    if (ierr /= 0) print *, &
         "Error dbaddiopt in SILO_add_grid: DBOPT_LO_OFFSET", ierr

    ierr = dbaddiopt(dboptix, DBOPT_HI_OFFSET, hi_offset)
    if (ierr /= 0) print *, &
         "Error dbaddiopt in SILO_add_grid: DBOPT_HI_OFFSET", ierr

    ! Write the grid structure
    ierr = dbputqm(dbix, trim(gridname), len_trim(gridname), &
         'x', 1, 'y', 1, 'z', 1, x_coords, y_coords, z_coords, &
         N_r, n_dim, DB_DOUBLE, DB_COLLINEAR, dboptix, iostat)
    if (ierr /= 0) print *, &
            "Error dbputqm is SILO_add_grid", ierr

    ierr = dbfreeoptlist(dboptix)
    if (ierr /= 0) print *, &
            "Error dbfreeoptlist is SILO_add_grid", ierr
  end subroutine SILO_add_grid

  subroutine SILO_add_var(dbix, dataname, gridname, &
       d_packed, d_shape)
    character(len=*), intent(in) :: gridname, dataname
    real(dp), intent(in)         :: d_packed(:)
    integer, intent(in)          :: dbix, d_shape(:)

    integer                      :: dboptix, ierr, iostat
    real(dp)                     :: dummy(1)

    interface
       function dbputqv1(dbid, name, lname, meshname, lmeshname, &
            var, dims, ndims, mixvar, mixlen, datatype, &
            centering, optlist_id, status)
         use, intrinsic :: iso_c_binding
         integer(c_int) :: dbid, lname, lmeshname, dims(*), ndims, mixlen
         integer(c_int) :: centering, optlist_id, status, datatype, dbputqv1
         real(c_double) :: var(*), mixvar(*)
         character(kind=c_char) :: name(*), meshname(*)
       end function dbputqv1
    end interface

    if (size(d_packed) /= product(d_shape)) then
       error stop "Error: d_packed does not correspond to d_shape"
    end if

    if (size(d_shape) < 1 .or. size(d_shape) > 3) then
       error stop "Error: size(d_shape) < 1 or size(d_shape) > 3"
    end if

    ierr = dbmkoptlist(10, dboptix)
    if (ierr /= 0) then
       error stop "Error creating options list in SILO_add_var"
    end if

    ierr = dbaddiopt(dboptix, DBOPT_HIDE_FROM_GUI, 1)

    ! Write the data to the grid
    ierr = dbputqv1(dbix, trim(dataname), len_trim(dataname), &
         trim(gridname), len_trim(gridname), d_packed, d_shape, &
         size(d_shape), dummy, 0, DB_DOUBLE, DB_ZONECENT, dboptix, iostat)
    if (ierr /= 0) then
       print *, "Error dbputqv1 in SILO_add_var", ierr
       error stop
    end if

    ierr = dbfreeoptlist(dboptix)
  end subroutine SILO_add_var

  subroutine SILO_set_mmesh_grid(dbix, mmname, gridnames, n_cycle, time)
    character(len=*), intent(in)   :: mmname, gridnames(:)
    integer, intent(in)            :: dbix
    integer, intent(in), optional  :: n_cycle
    real(dp), intent(in), optional :: time

    integer                        :: i, ierr,length
    integer                        :: dboptix, iostat, old_str_len
    integer                        :: n_grids, name_len, total_len
    integer, allocatable           :: m_types(:), name_lengths(:)
    character(len=:), allocatable  :: mnames

    interface
       function dbputmmesh(dbid, name, lname, nmesh, meshnames, lmeshnames, &
            meshtypes, optlist_id, status)
         use, intrinsic :: iso_c_binding
         integer(c_int) :: dbputmmesh, lmeshnames(*)
         integer(c_int) :: dbid, lname, nmesh, meshtypes(*), optlist_id, status
         character(kind=c_char) :: name(*), meshnames(*)
       end function dbputmmesh
    end interface

    n_grids = size(gridnames)
    if (n_grids < 1) then
       error stop "SILO_set_mmesh_grid: error too few grids (<1)"
    end if

    name_len  = len(gridnames(1))
    total_len = name_len * n_grids
    allocate(character(total_len) :: mnames)
    allocate(name_lengths(n_grids))
    allocate(m_types(n_grids))

    do i = 1, n_grids
       mnames((i-1)*name_len+1:i*name_len) = trim(gridnames(i)) // char(0)
    end do

    old_str_len  = dbset2dstrlen(name_len)
    m_types      = DB_QUAD_RECT
    name_lengths = name_len

    ierr = dbmkoptlist(10, dboptix)

    if (present(n_cycle)) then
       ierr = dbaddiopt(dboptix, DBOPT_CYCLE, n_cycle)
    end if

    if (present(time)) then
       ierr = dbaddiopt(dboptix, DBOPT_DTIME, time)
    end if

    ierr = dbputmmesh(dbix, trim(mmname), len_trim(mmname), n_grids, &
         mnames(1:total_len), name_lengths, m_types, dboptix, iostat)
    if (ierr /= 0) then
       error stop "Error calling dbputmmesh"
    end if

    ierr = dbfreeoptlist(dboptix)
    length = dbset2dstrlen(old_str_len)
  end subroutine SILO_set_mmesh_grid

  subroutine SILO_set_mmesh_var(dbix, mvname, mmname, datanames, n_cycle, time)
    character(len=*), intent(in)   :: mvname, mmname, datanames(:)
    integer, intent(in)            :: dbix
    integer, intent(in), optional  :: n_cycle
    real(dp), intent(in), optional :: time


    integer                        :: i, ierr, dboptix, iostat,length
    integer                        :: old_str_len, n_grids, name_len, total_len
    integer, allocatable           :: m_types(:), name_lengths(:)
    character(:), allocatable      :: dnames

    interface
       function dbputmvar(dbid, name, lname, nlevels, meshnames, &
            lmnames, meshtypes, optlist_id, status)
         use, intrinsic :: iso_c_binding
         integer(c_int) :: dbputmvar, lmnames(*)
         integer(c_int) :: dbid, lname, nlevels, meshtypes(*)
         integer(c_int) :: optlist_id, status
         character(kind=c_char) :: name(*), meshnames(*)
       end function dbputmvar
    end interface

    n_grids = size(datanames)
    if (n_grids < 1) then
       error stop "SILO_set_mmesh_var: error too few grids (<1)"
    end if

    name_len = len(datanames(1))
    total_len = name_len * n_grids
    allocate(character(total_len) :: dnames)
    allocate(name_lengths(n_grids))
    allocate(m_types(n_grids))

    do i = 1, n_grids
       dnames((i-1)*name_len+1:i*name_len) = trim(datanames(i)) // char(0)
    end do
    old_str_len  = dbset2dstrlen(name_len)
    m_types      = DB_QUADVAR
    name_lengths = name_len

    ierr = dbmkoptlist(10, dboptix)
    if (ierr /= 0) then
       error stop "Error creating options list in SILO_set_mmesh_var"
    end if

    if (present(n_cycle)) then
       ierr = dbaddiopt(dboptix, DBOPT_CYCLE, n_cycle)
    end if

    if (present(time)) then
       ierr = dbaddiopt(dboptix, DBOPT_DTIME, time)
    end if

    ierr = dbaddcopt(dboptix, DBOPT_MMESH_NAME, &
         trim(mmname), len_trim(mmname))
    if (ierr /= 0) print *, &
            "Error dbaddiopt is SILO_set_mmesh_var: DBOPT_MMESH_NAME", ierr

    ierr = dbputmvar(dbix, trim(mvname), len_trim(mvname), n_grids, &
         dnames(1:total_len), name_lengths, m_types, dboptix, iostat)
    if (ierr /= 0) then
       error stop "Error calling dbputmvar"
    end if

    ierr = dbfreeoptlist(dboptix)
    length = dbset2dstrlen(old_str_len)
  end subroutine SILO_set_mmesh_var

end module m_write_silo
