program test_unstr_vtk
  implicit none

  integer, parameter :: dp = kind(0.0d0)

  call test_unst()

  contains

    !> Subroutine for testing UnstructuredGrid functions.
    subroutine test_unst()
      use m_vtk

      integer, parameter      :: nn = 100
      integer, parameter      :: n_dim       = 2
      integer, parameter :: n_nodes = nn**2
      integer, parameter :: n_cells  = (nn-1)**2
      real(dp)               :: coords(0:n_dim*n_nodes-1)
      integer            :: cell_types(0:n_cells-1)
      integer            :: offsets(0:n_cells-1)
      integer            :: connects(0:n_cells * 4 - 1)
      real(dp)               :: v(0:n_nodes-1)

      integer :: i, j, ic, ix, row, col
      real :: dr(2) = [0.1, 0.1]

      type(vtk_t) :: vtkf

      do j = 0, nn-1
         do i = 0, nn-1
            ix = j * nn + i
            coords(2 * ix) = 1.2345 + i * dr(1)
            coords(2 * ix + 1) = 2.3456 + j * dr(2)
         end do
      end do

      do ic = 0, n_cells-1
         offsets(ic) = (ic+1) * 4
         row = ic/(nn-1)
         col = modulo(ic, nn-1)
         ix = row * nn + col
         connects(4*ic:4*ic+3) = (/ix, ix+1, ix+nn, ix+nn+1/)
      end do

      cell_types = 8
      v = 0

      call vtk_ini_xml(vtkf, 'test.vtu', 'UnstructuredGrid')
      call vtk_dat_xml(vtkf, "UnstructuredGrid", .true.)
      call vtk_geo_xml(vtkf, coords, n_nodes, n_cells, 2, 0, 0.0_dp)
      call vtk_con_xml(vtkf, connects, offsets, cell_types, n_cells)
      call vtk_dat_xml(vtkf, "PointData", .true.)
      call vtk_var_r8_xml(vtkf, 'scalars', v, n_nodes)
      call vtk_dat_xml(vtkf, "PointData", .false.)
      call vtk_geo_xml_close(vtkf)
      call vtk_dat_xml(vtkf, "UnstructuredGrid", .false.)
      call vtk_end_xml(vtkf)
    endsubroutine test_unst

end program test_unstr_vtk
