! This module contains routines for writing output files with Afivo. The Silo
! format should probably be used for larger files, especially in 3D.
!
! Author: Jannis Teunissen
! License: GPLv3

module m_a$D_output

  use m_a$D_types

  implicit none
  private

  public :: a$D_write_vtk
  public :: a$D_write_silo
  public :: a$D_write_plane

contains

  !> Write data in a plane (2D) to a VTK ASCII file. In 3D, r_min and r_max
  !> should have one identical coordinate (i.e., they differ in two
  !> coordinates).
  subroutine a$D_write_plane(tree, filename, ivs, r_min, r_max, n_pixels, dir)
    use m_a$D_interp, only: a$D_interp1, a$D_interp2
    type(a$D_t), intent(in)       :: tree        !< Tree to write out
    character(len=*), intent(in) :: filename    !< Filename for the vtk file
    integer, intent(in)          :: ivs(:)      !< Variables to write
    real(dp), intent(in)         :: r_min($D)   !< Minimum coordinate of plane
    real(dp), intent(in)         :: r_max($D)   !< Maximum coordinate of plane
    integer, intent(in)          :: n_pixels(2) !< Number of pixels in the plane
    character(len=*), optional, intent(in) :: dir !< Directory to place files in

    integer, parameter    :: my_unit = 100
    character(len=100)    :: fmt_string
    character(len=400)    :: fname
    integer               :: i, j, n_vars, dim_unused, n_points(3)
    real(dp)              :: r($D), dvec(3), origin(3)
    real(dp)              :: v1($D), v2($D)
    real(dp), allocatable :: pixel_data(:, :, :)

    n_vars      = size(ivs)
#if $D == 2
    origin(1:2) = r_min
    origin(3)   = 0
    dvec(1:2)   = r_max(1:2) - r_min(1:2)
    dvec(3)     = 0
    dim_unused  = 3
    n_points    = [n_pixels(1), n_pixels(2), 1]
    v1          = [dvec(1), 0.0_dp] / (n_pixels(1) - 1)
    v2          = [0.0_dp, dvec(2)] / (n_pixels(2) - 1)
#elif $D == 3
    origin      = r_min
    dvec        = r_max - r_min
    dim_unused  = minloc(abs(dvec), 1)

    select case (dim_unused)
    case (1)
       v1 = [0.0_dp, dvec(2), 0.0_dp] / (n_pixels(1) - 1)
       v2 = [0.0_dp, 0.0_dp, dvec(3)] / (n_pixels(2) - 1)
       n_points = [1, n_pixels(1), n_pixels(2)]
    case (2)
       v1 = [dvec(1), 0.0_dp, 0.0_dp] / (n_pixels(1) - 1)
       v2 = [0.0_dp, 0.0_dp, dvec(3)] / (n_pixels(2) - 1)
       n_points = [n_pixels(1), 1, n_pixels(2)]
    case (3)
       v1 = [dvec(1), 0.0_dp, 0.0_dp] / (n_pixels(1) - 1)
       v2 = [0.0_dp, dvec(2), 0.0_dp] / (n_pixels(2) - 1)
       n_points = [n_pixels(1), n_pixels(2), 1]
    end select
#endif

    allocate(pixel_data(n_vars, n_pixels(1), n_pixels(2)))

    do j = 1, n_pixels(2)
       do i = 1, n_pixels(1)
          r = r_min + (i-1) * v1 + (j-1) * v2
          pixel_data(:, i, j) = a$D_interp2(tree, r, ivs, n_vars)
       end do
    end do

    ! Construct format string. Write one row at a time
    write(fmt_string, '(A,I0,A)') '(', n_pixels(1), 'E16.8)'

    ! Construct file name
    fname = trim(filename) // ".vtk"
    if (present(dir)) then
       i = len_trim(dir)
       if (i > 0) then
          if (dir(i:i) == "/") then ! Dir has trailing slash
             fname = trim(dir) // trim(fname)
          else
             fname = trim(dir) // "/" // trim(fname)
          end if
       end if
    end if

    open(my_unit, file=trim(fname), action="write")
    write(my_unit, '(A)') "# vtk DataFile Version 2.0"
    write(my_unit, '(A)') trim(filename)
    write(my_unit, '(A)') "ASCII"
    write(my_unit, '(A)') "DATASET STRUCTURED_POINTS"
    write(my_unit, '(A,3I10)') "DIMENSIONS ", n_points
    write(my_unit, '(A,3E16.8)') "ORIGIN ", r_min
    write(my_unit, '(A,3E16.8)') "SPACING ", v1 + v2
    write(my_unit, '(A,2I0)') "POINT_DATA ", product(n_points)
    do i = 1, n_vars
       write(my_unit, '(A)') "SCALARS " // &
            trim(tree%cc_names(ivs(i))) // " double 1"
       write(my_unit, '(A)') "LOOKUP_TABLE default"
       write(my_unit, trim(fmt_string)) pixel_data(i, :, :)
    end do
    close(my_unit)
  end subroutine a$D_write_plane

  !> Write the cell centered data of a tree to a vtk unstructured file. Only the
  !> leaves of the tree are used
  subroutine a$D_write_vtk(tree, filename, n_cycle, time, ixs_cc, ixs_fc, dir)
    use m_vtk

    type(a$D_t), intent(in)       :: tree        !< Tree to write out
    character(len=*), intent(in)  :: filename    !< Filename for the vtk file
    integer, intent(in), optional :: n_cycle     !< Cycle-number for vtk file (counter)
    real(dp), intent(in), optional :: time        !< Time for output file
    integer, intent(in), optional :: ixs_cc(:)   !< Oncly include these cell variables
    character(len=*), optional, intent(in) :: dir !< Directory to place files in
    integer, intent(in), optional :: ixs_fc(:)   !< Oncly include these face variables

    integer                       :: lvl, bc, bn, n, n_cells, n_nodes
    integer                       :: ig, i, j, id, n_ix, c_ix, n_grids
    integer                       :: cell_ix, node_ix, n_cycle_val
    integer                       :: n_cc, n_fc
    integer, parameter            :: n_ch = a$D_num_children
    integer                       :: nodes_per_box, cells_per_box
    real(dp)                      :: time_val
    real(dp), allocatable         :: coords(:), cc_vars(:,:)
    integer, allocatable          :: offsets(:), connects(:)
    integer, allocatable          :: cell_types(:), icc_val(:), ifc_val(:)
    type(vtk_t)                   :: vtkf
    character(len=400)            :: fname
    character(len=100)            :: tmp_name
    character(len=100), allocatable :: var_names(:)
#if $D == 3
    integer                       :: k, bn2
#endif

    if (.not. tree%ready) stop "Tree not ready"
    time_val = 0.0_dp; if (present(time)) time_val = time
    n_cycle_val = 0; if (present(n_cycle)) n_cycle_val = n_cycle

    if (present(ixs_cc)) then
       if (maxval(ixs_cc) > tree%n_var_cell .or. &
            minval(ixs_cc) < 1) stop "a$D_write_silo: wrong indices given (ixs_cc)"
       icc_val = ixs_cc
    else
       icc_val = [(i, i = 1, tree%n_var_cell)]
    end if

    if (present(ixs_fc)) then
       if (maxval(ixs_fc) > tree%n_var_cell .or. &
            minval(ixs_fc) < 1) stop "a$D_write_silo: wrong indices given (ixs_fc)"
       ifc_val = ixs_fc
    else
       allocate(ifc_val(0))
    end if

    n_cc = size(icc_val)
    n_fc = size(ifc_val)

    allocate(var_names(n_cc + n_fc * $D))
    var_names(1:n_cc) = tree%cc_names(icc_val)

    do i = 1, n_fc
       n = n_cc + (i-1)*$D
       tmp_name = tree%fc_names(ifc_val(i))
#if $D == 2
       if (tree%coord_t == af_cyl) then
          var_names(n+1) = trim(tmp_name) // "_r"
          var_names(n+2) = trim(tmp_name) // "_z"
       else
          var_names(n+1) = trim(tmp_name) // "_x"
          var_names(n+2) = trim(tmp_name) // "_y"
       end if
#elif $D == 3
       var_names(n+1) = trim(tmp_name) // "_x"
       var_names(n+2) = trim(tmp_name) // "_y"
       var_names(n+3) = trim(tmp_name) // "_z"
#endif
    end do

    bc            = tree%n_cell     ! number of Box Cells
    bn            = tree%n_cell + 1 ! number of Box Nodes
    nodes_per_box = bn**$D
    cells_per_box = bc**$D

    n_grids = 0
    do lvl = 1, tree%highest_lvl
       n_grids = n_grids + size(tree%lvls(lvl)%leaves)
    end do
    n_nodes = nodes_per_box * n_grids
    n_cells = cells_per_box * n_grids

    allocate(coords($D * n_nodes))
    allocate(cc_vars(n_cells, n_cc + n_fc * $D))
    allocate(offsets(cells_per_box * n_grids))
    allocate(cell_types(cells_per_box * n_grids))
    allocate(connects(n_ch * cells_per_box * n_grids))

#if $D == 2
    cell_types = 8  ! VTK pixel type
#elif $D       == 3
    bn2        = bn**2
    cell_types = 11 ! VTK voxel type
#endif

    ig = 0
    do lvl = 1, tree%highest_lvl
       do n = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(n)

          ig = ig + 1
          cell_ix = (ig-1) * cells_per_box
          node_ix = (ig-1) * nodes_per_box

#if $D == 2
          do j = 1, bn
             do i = 1, bn
                n_ix = 2 * (node_ix + (j-1) * bn + i)
                coords(n_ix-1:n_ix) = tree%boxes(id)%r_min + &
                     [i-1,j-1] * tree%boxes(id)%dr
             end do
          end do

          do j = 1, bc
             do i = 1, bc
                ! In vtk, indexing starts at 0, so subtract 1
                n_ix                      = node_ix + (j-1) * bn + i - 1
                c_ix                      = cell_ix + (j-1) * bc + i
                cc_vars(c_ix, 1:n_cc)     = tree%boxes(id)%cc(i, j, icc_val)
                cc_vars(c_ix, n_cc+1::$D) = &
                     0.5_dp * sum(tree%boxes(id)%fx(i:i+1, j, ifc_val))
                cc_vars(c_ix, n_cc+2::$D) = &
                     0.5_dp * sum(tree%boxes(id)%fy(i, j:j+1, ifc_val))
                offsets(c_ix)             = a$D_num_children * c_ix
                connects(n_ch*(c_ix-1)+1:n_ch*c_ix) = [n_ix, n_ix+1, n_ix+bn, n_ix+bn+1]
             end do
          end do
#elif $D == 3
          do k = 1, bn
             do j = 1, bn
                do i = 1, bn
                   n_ix = 3 * (node_ix + (k-1) * bn2 + (j-1) * bn + i)
                   coords(n_ix-2:n_ix) = tree%boxes(id)%r_min + &
                        [i-1,j-1,k-1] * tree%boxes(id)%dr
                end do
             end do
          end do

          do k = 1, bc
             do j = 1, bc
                do i = 1, bc
                   ! In vtk, indexing starts at 0, so subtract 1
                   n_ix                      = node_ix + (k-1) * bn2 + &
                        (j-1) * bn + i - 1
                   c_ix                      = cell_ix + (k-1) * bc**2 + &
                        (j-1) * bc + i
                   cc_vars(c_ix, 1:n_cc)     = tree%boxes(id)%cc(i, j, k, icc_val)
                   cc_vars(c_ix, n_cc+1::$D) = &
                        0.5_dp * sum(tree%boxes(id)%fx(i:i+1, j, k, ifc_val))
                   cc_vars(c_ix, n_cc+2::$D) = &
                        0.5_dp * sum(tree%boxes(id)%fy(i, j:j+1, k, ifc_val))
                   cc_vars(c_ix, n_cc+3::$D) = &
                        0.5_dp * sum(tree%boxes(id)%fz(i, j, k:k+1, ifc_val))
                   offsets(c_ix)             = 8 * c_ix
                   connects(n_ch*(c_ix-1)+1:n_ch*c_ix) = &
                        [n_ix, n_ix+1, n_ix+bn, n_ix+bn+1, &
                        n_ix+bn2, n_ix+bn2+1, n_ix+bn2+bn, n_ix+bn2+bn+1]
                end do
             end do
          end do
#endif
       end do
    end do

    fname = trim(filename) // ".vtu"

    if (present(dir)) then
       i = len_trim(dir)
       if (i > 0) then
          if (dir(i:i) == "/") then ! Dir has trailing slash
             fname = trim(dir) // trim(fname)
          else
             fname = trim(dir) // "/" // trim(fname)
          end if
       end if
    end if

    call vtk_ini_xml(vtkf, trim(fname), 'UnstructuredGrid')
    call vtk_dat_xml(vtkf, "UnstructuredGrid", .true.)
    call vtk_unstr_geo_xml(vtkf, coords, n_nodes, n_cells, $D, n_cycle_val, time_val)
    call vtk_unstr_con_xml(vtkf, connects, offsets, cell_types, n_cells)
    call vtk_dat_xml(vtkf, "CellData", .true.)

    do n = 1, n_cc + n_fc * $D
       call vtk_var_r8_xml(vtkf, trim(var_names(n)), cc_vars(:, n), n_cells)
    end do

    call vtk_dat_xml(vtkf, "CellData", .false.)
    call vtk_unstr_geo_xml_close(vtkf)
    call vtk_dat_xml(vtkf, "UnstructuredGrid", .false.)
    call vtk_end_xml(vtkf)
    print *, "a$D_write_vtk: written " // trim(fname)
  end subroutine a$D_write_vtk

  !> Write the cell centered data of a tree to a Silo file. Only the
  !> leaves of the tree are used
  subroutine a$D_write_silo(tree, filename, n_cycle, time, ixs_cc, ixs_fc, dir)
    use m_write_silo

    type(a$D_t), intent(in)       :: tree        !< Tree to write out
    character(len=*)              :: filename    !< Filename for the vtk file
    integer, intent(in), optional :: n_cycle     !< Cycle-number for vtk file (counter)
    real(dp), intent(in), optional :: time        !< Time for output file
    integer, intent(in), optional :: ixs_cc(:)      !< Oncly include these cell variables
    character(len=*), optional, intent(in) :: dir !< Directory to place files in
    integer, intent(in), optional :: ixs_fc(:)      !< Oncly include these face variables

    character(len=*), parameter     :: grid_name = "gg"
    character(len=*), parameter     :: amr_name  = "mesh", meshdir = "data"
    character(len=100)              :: tmp_name
    character(len=100), allocatable :: grid_list(:), var_list(:, :), var_names(:)
    character(len=400)              :: fname
    integer                         :: lvl, i, id, i_grid, iv, nc, n_grids_max
    integer                         :: n, n_vars, i0, j0, dbix, n_cc, n_fc
    integer                         :: nx, ny, nx_prev, ny_prev, ix, iy
    integer                         :: n_cycle_val
    integer, allocatable            :: ids(:), nb_ids(:), icc_val(:), ifc_val(:)
    logical, allocatable            :: box_done(:)
    real(dp)                        :: dr($D), r_min($D), time_val
#if $D == 2
    integer, allocatable            :: box_list(:,:), new_box_list(:, :)
    real(dp), allocatable           :: var_data(:,:,:)
#elif $D == 3
    integer, allocatable            :: box_list(:,:,:), new_box_list(:,:,:)
    real(dp), allocatable           :: var_data(:,:,:,:)
    integer                         :: k0, nz, nz_prev, iz
#endif

    if (.not. tree%ready) stop "Tree not ready"
    time_val = 0.0_dp; if (present(time)) time_val = time
    n_cycle_val = 0; if (present(n_cycle)) n_cycle_val = n_cycle

    if (present(ixs_cc)) then
       if (maxval(ixs_cc) > tree%n_var_cell .or. &
            minval(ixs_cc) < 1) stop "a$D_write_silo: wrong indices given (ixs_cc)"
       icc_val = ixs_cc
    else
       icc_val = [(i, i = 1, tree%n_var_cell)]
    end if

    if (present(ixs_fc)) then
       if (maxval(ixs_fc) > tree%n_var_cell .or. &
            minval(ixs_fc) < 1) stop "a$D_write_silo: wrong indices given (ixs_fc)"
       ifc_val = ixs_fc
    else
       allocate(ifc_val(0))
    end if

    n_cc = size(icc_val)
    n_fc = size(ifc_val)

    allocate(var_names(n_cc + n_fc * $D))
    var_names(1:n_cc) = tree%cc_names(icc_val)

    do i = 1, n_fc
       n = n_cc + (i-1)*$D
       tmp_name = tree%fc_names(ifc_val(i))
#if $D == 2
       if (tree%coord_t == af_cyl) then
          var_names(n+1) = trim(tmp_name) // "_r"
          var_names(n+2) = trim(tmp_name) // "_z"
       else
          var_names(n+1) = trim(tmp_name) // "_x"
          var_names(n+2) = trim(tmp_name) // "_y"
       end if
#elif $D == 3
       var_names(n+1) = trim(tmp_name) // "_x"
       var_names(n+2) = trim(tmp_name) // "_y"
       var_names(n+3) = trim(tmp_name) // "_z"
#endif
    end do

    nc = tree%n_cell
    n_vars = n_cc + n_fc * $D
    n_grids_max = 0
    do lvl = 1, tree%highest_lvl
       n_grids_max = n_grids_max + size(tree%lvls(lvl)%leaves)
    end do

    allocate(grid_list(n_grids_max))
    allocate(var_list(n_vars, n_grids_max))
    allocate(box_done(tree%highest_id))
    box_done = .false.

    fname = trim(filename) // ".silo"

    if (present(dir)) then
       i = len_trim(dir)
       if (i > 0) then
          if (dir(i:i) == "/") then ! Dir has trailing slash
             fname = trim(dir) // trim(fname)
          else
             fname = trim(dir) // "/" // trim(fname)
          end if
       end if
    end if

    call SILO_create_file(trim(fname), dbix)
    call SILO_mkdir(dbix, meshdir)
    i_grid = 0

    do lvl = 1, tree%highest_lvl
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          if (box_done(id)) cycle

          i_grid = i_grid + 1

          ! Find largest rectangular box including id and other leaves that
          ! haven't been written yet
#if $D == 2
          allocate(box_list(1,1))
          box_list(1,1) = id
          box_done(id) = .true.
          nx = 1
          ny = 1

          do
             nx_prev = nx
             ny_prev = ny

             ! Check whether we can extend to the -x direction
             ids = box_list(1, :)
             nb_ids = tree%boxes(ids)%neighbors(a2_neighb_lowx)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(1) < tree%boxes(ids(1))%ix(1) .and. &
                     .not. any(a2_has_children(tree%boxes(nb_ids)))) then
                   nx = nx + 1
                   allocate(new_box_list(nx, ny))
                   new_box_list(1, :) = nb_ids
                   new_box_list(2:, :) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the +x direction
             ids = box_list(nx, :)
             nb_ids = tree%boxes(ids)%neighbors(a2_neighb_highx)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(1) > tree%boxes(ids(1))%ix(1) .and. &
                     .not. any(a2_has_children(tree%boxes(nb_ids)))) then
                   nx = nx + 1
                   allocate(new_box_list(nx, ny))
                   new_box_list(nx, :) = nb_ids
                   new_box_list(1:nx-1, :) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the -y direction
             ids = box_list(:, 1)
             nb_ids = tree%boxes(ids)%neighbors(a2_neighb_lowy)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(2) < tree%boxes(ids(1))%ix(2) .and. &
                     .not. any(a2_has_children(tree%boxes(nb_ids)))) then
                   ny = ny + 1
                   allocate(new_box_list(nx, ny))
                   new_box_list(:, 1) = nb_ids
                   new_box_list(:, 2:) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the +y direction
             ids = box_list(:, ny)
             nb_ids = tree%boxes(ids)%neighbors(a2_neighb_highy)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(2) > tree%boxes(ids(1))%ix(2) .and. &
                     .not. any(a2_has_children(tree%boxes(nb_ids)))) then
                   ny = ny + 1
                   allocate(new_box_list(nx, ny))
                   new_box_list(:, ny) = nb_ids
                   new_box_list(:, 1:ny-1) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             if (nx == nx_prev .and. ny == ny_prev) exit
          end do

          allocate(var_data(nx * nc, ny * nc, n_vars))
          do ix = 1, nx
             do iy = 1, ny
                id = box_list(ix, iy)
                i0 = 1 + (ix-1) * nc
                j0 = 1 + (iy-1) * nc
                var_data(i0:i0+nc-1, j0:j0+nc-1, 1:n_cc) = &
                     tree%boxes(id)%cc(1:nc, 1:nc, icc_val)
                var_data(i0:i0+nc-1, j0:j0+nc-1, n_cc+1::$D) = &
                     0.5_dp * (tree%boxes(id)%fx(1:nc, 1:nc, ifc_val) + &
                     tree%boxes(id)%fx(2:nc+1, 1:nc, ifc_val))
                var_data(i0:i0+nc-1, j0:j0+nc-1, n_cc+2::$D) = &
                     0.5_dp * (tree%boxes(id)%fy(1:nc, 1:nc, ifc_val) + &
                     tree%boxes(id)%fy(1:nc, 2:nc+1, ifc_val))
             end do
          end do

          id = box_list(1, 1)
          dr = tree%boxes(id)%dr
          r_min = tree%boxes(id)%r_min

          write(grid_list(i_grid), "(A,I0)") meshdir // '/' // grid_name, i_grid
          call SILO_add_grid(dbix, grid_list(i_grid), 2, &
               [nx*nc + 1, ny*nc + 1], r_min, dr)
          do iv = 1, n_vars
             write(var_list(iv, i_grid), "(A,I0)") meshdir // '/' // &
                  trim(var_names(iv)) // "_", i_grid
             call SILO_add_var(dbix, var_list(iv, i_grid), grid_list(i_grid), &
                  pack(var_data(:, :, iv), .true.), [nx*nc, ny*nc])
          end do

          deallocate(var_data)
          deallocate(box_list)
#elif $D == 3
          allocate(box_list(1,1,1))
          box_list(1,1,1) = id
          nx = 1
          ny = 1
          nz = 1

          do
             nx_prev = nx
             ny_prev = ny
             nz_prev = nz

             ! Check whether we can extend to the -x direction
             ids = pack(box_list(1, :, :), .true.)
             nb_ids = tree%boxes(ids)%neighbors(a3_neighb_lowx)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(1) < tree%boxes(ids(1))%ix(1) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   nx = nx + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(1, :, :) = reshape(nb_ids, [ny, nz])
                   new_box_list(2:, :, :) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the +x direction
             ids = pack(box_list(nx, :, :), .true.)
             nb_ids = tree%boxes(ids)%neighbors(a3_neighb_highx)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(1) > tree%boxes(ids(1))%ix(1) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   nx = nx + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(nx, :, :) = reshape(nb_ids, [ny, nz])
                   new_box_list(1:nx-1, :, :) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the -y direction
             ids = pack(box_list(:, 1, :), .true.)
             nb_ids = tree%boxes(ids)%neighbors(a3_neighb_lowy)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(2) < tree%boxes(ids(1))%ix(2) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   ny = ny + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(:, 1, :) = reshape(nb_ids, [nx, nz])
                   new_box_list(:, 2:, :) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the +y direction
             ids = pack(box_list(:, ny, :), .true.)
             nb_ids = tree%boxes(ids)%neighbors(a3_neighb_highy)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(2) > tree%boxes(ids(1))%ix(2) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   ny = ny + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(:, ny, :) = reshape(nb_ids, [nx, nz])
                   new_box_list(:, 1:ny-1, :) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the -z direction
             ids = pack(box_list(:, :, 1), .true.)
             nb_ids = tree%boxes(ids)%neighbors(a3_neighb_lowz)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(3) < tree%boxes(ids(1))%ix(3) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   nz = nz + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(:, :, 1) = reshape(nb_ids, [nx, ny])
                   new_box_list(:, :, 2:) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the +z direction
             ids = pack(box_list(:, :, nz), .true.)
             nb_ids = tree%boxes(ids)%neighbors(a3_neighb_highz)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(3) > tree%boxes(ids(1))%ix(3) .and. &
                     .not. any(a3_has_children(tree%boxes(nb_ids)))) then
                   nz = nz + 1
                   allocate(new_box_list(nx, ny, nz))
                   new_box_list(:, :, nz) = reshape(nb_ids, [nx, ny])
                   new_box_list(:, :, 1:nz-1) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             if (nx == nx_prev .and. ny == ny_prev .and. nz == nz_prev) exit
          end do

          allocate(var_data(nx * nc, ny * nc, nz * nc, n_vars))
          do iz = 1, nz
             do ix = 1, nx
                do iy = 1, ny
                   id = box_list(ix, iy, iz)
                   i0 = 1 + (ix-1) * nc
                   j0 = 1 + (iy-1) * nc
                   k0 = 1 + (iz-1) * nc
                   var_data(i0:i0+nc-1, j0:j0+nc-1, k0:k0+nc-1, 1:n_cc) = &
                        tree%boxes(id)%cc(1:nc, 1:nc, 1:nc, icc_val)
                   var_data(i0:i0+nc-1, j0:j0+nc-1, k0:k0+nc-1, n_cc+1::$D) = &
                        0.5_dp * (tree%boxes(id)%fx(1:nc, 1:nc, 1:nc, ifc_val) + &
                        tree%boxes(id)%fx(2:nc+1, 1:nc, 1:nc, ifc_val))
                   var_data(i0:i0+nc-1, j0:j0+nc-1, k0:k0+nc-1, n_cc+2::$D) = &
                        0.5_dp * (tree%boxes(id)%fy(1:nc, 1:nc, 1:nc, ifc_val) + &
                        tree%boxes(id)%fy(1:nc, 2:nc+1, 1:nc, ifc_val))
                   var_data(i0:i0+nc-1, j0:j0+nc-1, k0:k0+nc-1, n_cc+3::$D) = &
                        0.5_dp * (tree%boxes(id)%fz(1:nc, 1:nc, 1:nc, ifc_val) + &
                        tree%boxes(id)%fz(1:nc, 1:nc, 2:nc+1, ifc_val))
                end do
             end do
          end do

          id = box_list(1, 1, 1)
          dr = tree%boxes(id)%dr
          r_min = tree%boxes(id)%r_min

          write(grid_list(i_grid), "(A,I0)") meshdir // '/' // grid_name, i_grid
          call SILO_add_grid(dbix, grid_list(i_grid), 3, &
               [nx*nc + 1, ny*nc + 1, nz*nc + 1], r_min, dr)
          do iv = 1, n_vars
             write(var_list(iv, i_grid), "(A,I0)") meshdir // '/' // &
                  trim(var_names(iv)) // "_", i_grid
             call SILO_add_var(dbix, var_list(iv, i_grid), grid_list(i_grid), &
                  pack(var_data(:, :, :, iv), .true.), [nx*nc, ny*nc, nz*nc])
          end do

          deallocate(var_data)
          deallocate(box_list)
#endif

       end do
    end do

    call SILO_set_mmesh_grid(dbix, amr_name, grid_list(1:i_grid), &
         n_cycle_val, time_val)
    do iv = 1, n_vars
       call SILO_set_mmesh_var(dbix, trim(var_names(iv)), amr_name, &
            var_list(iv, 1:i_grid), n_cycle_val, time_val)
    end do
    call SILO_close_file(dbix)
    print *, "a$D_write_silo: written " // trim(fname)
  end subroutine a$D_write_silo

end module m_a$D_output
