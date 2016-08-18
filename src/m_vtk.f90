!> This file is a modification of code found in Lib_VTK_IO
!>
!> For the original code, see https://github.com/szaghi/Lib_VTK_IO. I have
!> extracted the parts that I needed and simplified them a bit. The license for
!> this file is thus also GPLv3.
!>
!> Author Lib_VTK_IO: Stefano Zaghi, modifications: Jannis Teunissen
!> \todo Document this module
module m_vtk

  implicit none
  public

  integer, parameter   :: dp            = kind(0.0d0)
  integer, parameter   :: sp            = kind(0.0)
  integer, parameter   :: buf_len       = 200
  ! Ifort does not yet support new_line("c"), so for now...
  character, parameter :: endl          = char(10)
  integer, parameter   :: indent_spaces = 2
  integer, parameter   :: bytes_i4      = 4
  integer, parameter   :: bytes_r4      = 4
  integer, parameter   :: bytes_r8      = 8

  type vtk_t
     character(len=100) :: filename
     logical            :: is_open
     integer            :: funit ! xml file
     integer            :: sunit ! scratch file
     integer            :: indent
     integer            :: offset
  end type vtk_t

contains

  subroutine vtk_ini_xml(vtkf, filename, vtk_type)
    type(vtk_t), intent(out)     :: vtkf
    character(len=*), intent(in) :: filename, vtk_type

    ! Initialize vtk type
    open(newunit=vtkf%funit, file=trim(filename), form='UNFORMATTED', &
         access='STREAM', status='REPLACE')
    open(newunit=vtkf%sunit, form='UNFORMATTED', &
         access='STREAM', status='SCRATCH')

    vtkf%filename = filename
    vtkf%is_open  = .true.
    vtkf%indent   = 0
    vtkf%offset   = 0

    write(vtkf%funit) '<?xml version="1.0"?>' // endl
    write(vtkf%funit) '<VTKFile type="' // trim(vtk_type) // &
         '" version="0.1" byte_order="LittleEndian">' // endl
    vtkf%indent = vtkf%indent + indent_spaces
  end subroutine vtk_ini_xml

  subroutine vtk_unstr_geo_xml(vtkf, coords, n_nodes, n_cells, n_dim, cycl, time)
    type(vtk_t), intent(inout) :: vtkf
    real(dp), intent(in)       :: coords(:), time
    real(sp), allocatable      :: wr_coords(:)
    integer, intent(in)        :: n_nodes, n_cells, n_dim, cycl
    character(len=buf_len)     :: bufstr
    integer                    :: n_bytes, d

    if (n_dim < 1 .or. n_dim > 3) stop "n_dim should be between 1-3"

    ! Always write 3-d point coordinates
    allocate(wr_coords(3 * n_nodes))
    wr_coords = 0

    call vtk_dat_xml(vtkf, "FieldData", .true.)
    write(vtkf%funit) repeat(' ',vtkf%indent) // '<DataArray type="Float64"'// &
         ' Name="TIME" NumberOfTuples="1" format="ascii">' // endl
    vtkf%indent = vtkf%indent + indent_spaces
    write(bufstr, *) repeat(' ',vtkf%indent), time
    write(vtkf%funit) trim(bufstr) // endl
    call vtk_dat_xml(vtkf, "DataArray", .false.)
    write(vtkf%funit) repeat(' ',vtkf%indent) // '<DataArray type="Int32"'// &
         ' Name="CYCLE" NumberOfTuples="1" format="ascii">' // endl
    vtkf%indent = vtkf%indent + indent_spaces
    write(bufstr, *) repeat(' ',vtkf%indent), cycl
    write(vtkf%funit) trim(bufstr) // endl
    call vtk_dat_xml(vtkf, "DataArray", .false.)
    call vtk_dat_xml(vtkf, "FieldData", .false.)

    do d = 1, n_dim
       wr_coords(d::3) = real(coords(d::n_dim), sp)
    end do

    write(bufstr, fmt="(A,I0,A,I0,A)") repeat(' ',vtkf%indent) // &
         '<Piece NumberOfPoints="', n_nodes, '" NumberOfCells="', n_cells, '">'
    write(vtkf%funit) trim(bufstr) // endl
    vtkf%indent = vtkf%indent + indent_spaces

    call vtk_dat_xml(vtkf, "Points", .true.)

    write(bufstr, fmt="(A,I0,A,I0,A)") repeat(' ',vtkf%indent) // &
         '<DataArray type="Float32" NumberOfComponents="3" Name="Points"' // &
         ' format="appended" offset="', vtkf%offset, '"/>'
    write(vtkf%funit) trim(bufstr) // endl

    ! Write raw data
    n_bytes     = 3 * n_nodes * bytes_r4
    vtkf%offset = vtkf%offset + bytes_i4 + n_bytes
    write(vtkf%sunit) n_bytes, wr_coords

    call vtk_dat_xml(vtkf, "Points", .false.)
  end subroutine vtk_unstr_geo_xml

  subroutine vtk_unstr_geo_xml_close(vtkf)
    type(vtk_t), intent(inout) :: vtkf
    call vtk_dat_xml(vtkf, "Piece", .false.)
  end subroutine vtk_unstr_geo_xml_close

  subroutine vtk_unstr_con_xml(vtkf, connects, offsets, cell_types, n_cells)
    type(vtk_t), intent(inout) :: vtkf
    integer, intent(IN)        :: n_cells       !< Number of cells.
    integer, intent(IN)        :: connects(:)   !< Mesh connectivity.
    integer, intent(IN)        :: offsets(:)    !< Cell offset.
    integer, intent(IN)        :: cell_types(:) !< VTK cell type.
    character(len=buf_len)     :: bufstr
    integer                    :: n_bytes

    call vtk_dat_xml(vtkf, "Cells", .true.)

    write(bufstr, fmt="(A,I0,A)") repeat(' ',vtkf%indent) // &
         '<DataArray type="Int32" Name="connectivity" format="appended" offset="', &
         vtkf%offset, '"/>'
    write(vtkf%funit) trim(bufstr) // endl

    n_bytes = offsets(n_cells) * bytes_i4
    vtkf%offset = vtkf%offset + bytes_i4 + n_bytes
    write(vtkf%sunit) n_bytes, connects

    write(bufstr, fmt="(A,I0,A)") repeat(' ',vtkf%indent) // &
         '<DataArray type="Int32" Name="offsets" format="appended" offset="', &
         vtkf%offset, '"/>'
    write(vtkf%funit) trim(bufstr) // endl

    n_bytes = n_cells * bytes_i4
    vtkf%offset = vtkf%offset + bytes_i4 + n_bytes
    write(vtkf%sunit) n_bytes, offsets

    write(bufstr, fmt="(A,I0,A)") repeat(' ',vtkf%indent) // &
         '<DataArray type="Int32" Name="types" format="appended" offset="', &
         vtkf%offset, '"/>'
    write(vtkf%funit) trim(bufstr) // endl

    n_bytes = n_cells * bytes_i4
    vtkf%offset = vtkf%offset + bytes_i4 + n_bytes
    write(vtkf%sunit) n_bytes, cell_types

    call vtk_dat_xml(vtkf, "Cells", .false.)
  end subroutine vtk_unstr_con_xml

  subroutine vtk_dat_xml(vtkf, xml_name, true_is_open)
    type(vtk_t), intent(inout) :: vtkf
    character(*), intent(IN)   :: xml_name
    logical, intent(in)        :: true_is_open

    if (true_is_open) then
       write(vtkf%funit) &
            repeat(' ', vtkf%indent) // '<' // trim(xml_name) // '>' // endl
       vtkf%indent = vtkf%indent + indent_spaces
    else
       vtkf%indent = vtkf%indent - indent_spaces
       write(vtkf%funit) &
            repeat(' ', vtkf%indent) // '</' // trim(xml_name) // '>' // endl
    end if
  end subroutine vtk_dat_xml

  subroutine vtk_var_r8_xml(vtkf, varname, var, n_data)
    type(vtk_t), intent(inout) :: vtkf
    integer, intent(IN)        :: n_data  !< Number of cells or nodes.
    character(*), intent(IN)   :: varname !< Variable name.
    real(dp),    intent(IN)    :: var(:)  !< Variable to be saved [1:n_cells_NN].
    character(len=buf_len)     :: bufstr
    integer                    :: n_bytes

    write(bufstr, fmt="(A,I0,A)") repeat(' ',vtkf%indent) // &
         '<DataArray type="Float64" Name="' // trim(varname) // &
         '" NumberOfComponents="1" format="appended" offset="', &
         vtkf%offset, '"/>'
    write(vtkf%funit) trim(bufstr) // endl
    n_bytes     = n_data * bytes_r8
    vtkf%offset = vtkf%offset + bytes_i4 + n_bytes
    write(vtkf%sunit) n_bytes, var
  end subroutine vtk_var_r8_xml

  subroutine vtk_end_xml(vtkf)
    use, intrinsic :: iso_c_binding
    type(vtk_t), intent(inout) :: vtkf
    integer :: n_bytes
    integer, allocatable :: buffer(:)

    write(vtkf%funit) &
         repeat(' ',vtkf%indent) // '<AppendedData encoding="raw">' // endl
    inquire(vtkf%sunit, pos=n_bytes)
    n_bytes = n_bytes - 1
    rewind(vtkf%sunit)
    allocate(buffer(n_bytes/bytes_i4))

    read(vtkf%sunit) buffer(:)
    write(vtkf%funit) '_', buffer, endl
    write(vtkf%funit) &
         repeat(' ',vtkf%indent) // '</AppendedData>' // endl
    call vtk_dat_xml(vtkf, "VTKFile", .false.)
    close(vtkf%sunit)
    close(vtkf%funit)
    vtkf%is_open = .false.
  end subroutine vtk_end_xml

end module m_vtk
