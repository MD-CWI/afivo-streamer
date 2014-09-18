module m_write_silo

   implicit none
   private

   include 'silo_f9x.inc'

   integer, parameter :: dp       = kind(0.0d0)
   integer, parameter :: line_len = 200
   integer, parameter :: DB_TYPE  = DB_PDB

   public :: SILO_create_file
   public :: SILO_add_grid
   public :: SILO_add_var
   public :: SILO_set_mmesh_grid
   public :: SILO_set_mmesh_var

contains

   subroutine SILO_create_file(filename)
     character(len=*), intent(in) :: filename
     integer                      :: ierr, dbix
     character(len=line_len)      :: fileinfo

     fileinfo = "A silo file"
      ierr = dbcreate(trim(filename), len_trim(filename), DB_CLOBBER, &
           DB_LOCAL, fileinfo, len_trim(fileinfo), DB_TYPE, dbix)
      if (ierr /= 0) print *, "Error creating file", trim(filename)
      call close_file(dbix)
      print *, "Created file ", trim(filename)
   end subroutine SILO_create_file

   function open_file(filename) result(dbix)
      character(len=*), intent(in) :: filename
      integer :: dbix, ierr

      ierr = dbopen(trim(filename), len_trim(filename), DB_TYPE, &
           DB_APPEND, dbix)
      if (ierr /= 0) print *, "Error opening file", trim(filename)
   end function open_file

   subroutine close_file(dbix)
      integer, intent(in) :: dbix
      integer :: ierr

      ierr = dbclose(dbix)
      if (ierr /= 0) print *, "Error closing file with index", dbix
   end subroutine close_file

   subroutine set_dir_create(dbix, dirname)
      character(len=*), intent(in) :: dirname
      integer, intent(in) :: dbix
      integer :: ierr, iostat

      ierr = dbsetdir(dbix, trim(dirname), len_trim(dirname))
      if (ierr /= 0) then
         ierr = dbmkdir(dbix, trim(dirname), len_trim(dirname), iostat)
         if (ierr /= 0) print *, "Error creating directory ", dirname
         ierr = dbsetdir(dbix, trim(dirname), len_trim(dirname))
         if (ierr /= 0) print *, "Could not go to directory ", dirname
      end if
   end subroutine set_dir_create

   subroutine SILO_add_grid(filename, gridname, n_dim, N_r, r_min, dr)
      character(len=*), intent(in) :: filename, gridname
      integer, intent(in)          :: n_dim, N_r(:)
      real(dp), intent(in)         :: r_min(:), dr(:)

      real(dp), allocatable        :: x_coords(:), y_coords(:), z_coords(:)
      integer                      :: i, ierr, iostat, dboptix, dbix

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

      ! Set integer options
      ierr = dbaddiopt(dboptix, DBOPT_NSPACE, n_dim)
      ierr = dbaddiopt(dboptix, DBOPT_HIDE_FROM_GUI, 1)
      ! ierr = dbaddiopt(dboptix, DBOPT_MAJORORDER, 1)

      dbix = open_file(filename)

      ! Write the grid structure
      ierr = dbputqm(dbix, trim(gridname), len_trim(gridname), &
           'x', 1, 'y', 1, 'z', 1, x_coords, y_coords, z_coords, &
           N_r, n_dim, DB_DOUBLE, DB_COLLINEAR, dboptix, iostat)

      call close_file(dbix)
      ierr = dbfreeoptlist(dboptix)
   end subroutine SILO_add_grid

   subroutine SILO_add_var(filename, dataname, gridname, &
        d_packed, d_shape, d_unit)
      character(len=*), intent(in) :: filename, gridname, dataname, d_unit
      real(dp), intent(in)         :: d_packed(:)
      integer, intent(in)          :: d_shape(:)

      integer                      :: dboptix, ierr, iostat, dbix
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
         print *, "Error: d_packed does not correspond to d_shape"
         return
      end if

      if (size(d_shape) < 1 .or. size(d_shape) > 3) then
         print *, "Error: size(d_shape) < 1 or size(d_shape) > 3"
         return
      end if

      ierr = dbmkoptlist(10, dboptix)
      ierr = dbaddcopt(dboptix, DBOPT_UNITS, trim(d_unit), len_trim(d_unit))
      ierr = dbaddiopt(dboptix, DBOPT_HIDE_FROM_GUI, 1)
      !       ierr = dbaddiopt(dboptix, DBOPT_MAJORORDER, 1)

      ! Get ix corresponding to filename
      dbix = open_file(filename)

      ! Write the data to the grid
      ierr = dbputqv1(dbix, trim(dataname), len_trim(dataname), &
           trim(gridname), len_trim(gridname), d_packed, d_shape, &
           size(d_shape), dummy, 0, DB_DOUBLE, DB_ZONECENT, dboptix, iostat)

      ! Close the file and remove option list
      call close_file(dbix)
      ierr = dbfreeoptlist(dboptix)
   end subroutine SILO_add_var

   subroutine SILO_set_mmesh_grid(filename, mmname, gridnames, n_cycle, time)
     character(len=*), intent(in)   :: filename, mmname, gridnames(:)
     integer, intent(in), optional  :: n_cycle
     real(dp), intent(in), optional :: time

     integer, parameter             :: long_len = 5000
     integer                        :: i, ierr
     integer                        :: dbix, dboptix, iostat, old_str_len
     integer                        :: n_grids, name_len, total_len
     integer, allocatable           :: m_types(:), name_lengths(:)
     character(:), allocatable      :: mnames

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
         print *, "Error, too few grids (<1)"
         return
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
      m_types      = DB_QUADMESH
      name_lengths = name_len

      dbix = open_file(filename)
      ierr = dbmkoptlist(10, dboptix)
      if (present(n_cycle)) ierr = dbaddiopt(dboptix, DBOPT_CYCLE, n_cycle)
      if (present(time)) ierr = dbaddiopt(dboptix, DBOPT_DTIME, time)

      ierr = dbputmmesh(dbix, trim(mmname), len_trim(mmname), n_grids, &
           mnames(1:total_len), name_lengths, m_types, dboptix, iostat)
      if (ierr /= 0) print *, "Error calling dbputmmesh", ierr

      call close_file(dbix)
      ierr = dbfreeoptlist(dboptix)
      ierr = dbset2dstrlen(old_str_len)
   end subroutine SILO_set_mmesh_grid

   subroutine SILO_set_mmesh_var(filename, mvname, mmname, &
        datanames, n_cycle, time)
     character(len=*), intent(in)   :: filename, mvname, mmname, datanames(:)
     integer, intent(in), optional  :: n_cycle
     real(dp), intent(in), optional :: time

     integer, parameter             :: long_len = 5000
     integer                        :: i, ierr, dbix, dboptix, iostat
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
         print *, "Error, too few grids (<1)"
         return
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

      dbix = open_file(filename)
      ierr = dbmkoptlist(10, dboptix)
      if (present(n_cycle)) ierr = dbaddiopt(dboptix, DBOPT_CYCLE, n_cycle)
      if (present(time)) ierr = dbaddiopt(dboptix, DBOPT_DTIME, time)
      ierr = dbaddcopt(dboptix, DBOPT_MMESH_NAME, &
           trim(mmname), len_trim(mmname))

      ierr = dbputmvar(dbix, trim(mvname), len_trim(mvname), n_grids, &
           dnames(1:total_len), name_lengths, m_types, dboptix, iostat)
      if (ierr /= 0) print *, "Error calling dbputmvar", ierr

      call close_file(dbix)
      ierr = dbfreeoptlist(dboptix)
      ierr = dbset2dstrlen(old_str_len)
   end subroutine SILO_set_mmesh_var

end module m_write_silo
