!> This module contains routines for writing output files with Afivo. The Silo
!> format should probably be used for larger files, especially in 3D.
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
    use m_a$D_interp, only: a$D_interp1
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
    integer               :: i, j, n_cc, dim_unused, n_points(3)
    real(dp)              :: r($D), dvec(3)
    real(dp)              :: v1($D), v2($D)
    real(dp), allocatable :: pixel_data(:, :, :)

    n_cc      = size(ivs)
#if $D == 2
    dvec(1:2)   = r_max(1:2) - r_min(1:2)
    dvec(3)     = 0
    dim_unused  = 3
    n_points    = [n_pixels(1), n_pixels(2), 1]
    v1          = [dvec(1), 0.0_dp] / (n_pixels(1) - 1)
    v2          = [0.0_dp, dvec(2)] / (n_pixels(2) - 1)
#elif $D == 3
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

    allocate(pixel_data(n_cc, n_pixels(1), n_pixels(2)))

    do j = 1, n_pixels(2)
       do i = 1, n_pixels(1)
          r = r_min + (i-1) * v1 + (j-1) * v2
          pixel_data(:, i, j) = a$D_interp1(tree, r, ivs, n_cc)
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
#if $D == 2
    write(my_unit, '(A,3E16.8)') "ORIGIN ", [r_min(1), r_min(2), 0.0_dp]
    write(my_unit, '(A,3E16.8)') "SPACING ", &
         [v1(1) + v2(1), v1(2) + v2(2), 0.0_dp]
#elif $D == 3
    write(my_unit, '(A,3E16.8)') "ORIGIN ", r_min
    write(my_unit, '(A,3E16.8)') "SPACING ", v1 + v2
#endif
    write(my_unit, '(A,2I0)') "POINT_DATA ", product(n_points)
    do i = 1, n_cc
       write(my_unit, '(A)') "SCALARS " // &
            trim(tree%cc_names(ivs(i))) // " double 1"
       write(my_unit, '(A)') "LOOKUP_TABLE default"
       write(my_unit, trim(fmt_string)) pixel_data(i, :, :)
    end do
    close(my_unit)
  end subroutine a$D_write_plane

  !> Write the cell centered data of a tree to a vtk unstructured file. Only the
  !> leaves of the tree are used
  subroutine a$D_write_vtk(tree, filename, n_cycle, time, ixs_cc, dir)
    use m_vtk

    type(a$D_t), intent(in)       :: tree        !< Tree to write out
    character(len=*), intent(in)  :: filename    !< Filename for the vtk file
    integer, intent(in), optional :: n_cycle     !< Cycle-number for vtk file (counter)
    real(dp), intent(in), optional :: time        !< Time for output file
    integer, intent(in), optional :: ixs_cc(:)   !< Oncly include these cell variables
    character(len=*), optional, intent(in) :: dir !< Directory to place files in

    integer                       :: lvl, bc, bn, n, n_cells, n_nodes
    integer                       :: ig, i, j, id, n_ix, c_ix, n_grids
    integer                       :: cell_ix, node_ix, n_cycle_val
    integer                       :: n_cc
    integer, parameter            :: n_ch = a$D_num_children
    integer                       :: nodes_per_box, cells_per_box
    real(dp)                      :: time_val
    real(dp), allocatable         :: coords(:), cc_vars(:,:)
    integer, allocatable          :: offsets(:), connects(:)
    integer, allocatable          :: cell_types(:), icc_val(:)
    type(vtk_t)                   :: vtkf
    character(len=400)            :: fname
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

    n_cc = size(icc_val)

    allocate(var_names(n_cc))
    var_names(1:n_cc) = tree%cc_names(icc_val)

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
    allocate(cc_vars(n_cells, n_cc))
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

    do n = 1, n_cc
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
  subroutine a$D_write_silo(tree, filename, n_cycle, time, ixs_cc, dir)
    use m_write_silo

    type(a$D_t), intent(in)       :: tree        !< Tree to write out
    character(len=*)              :: filename    !< Filename for the vtk file
    integer, intent(in), optional :: n_cycle     !< Cycle-number for vtk file (counter)
    real(dp), intent(in), optional :: time        !< Time for output file
    integer, intent(in), optional :: ixs_cc(:)      !< Oncly include these cell variables
    character(len=*), optional, intent(in) :: dir !< Directory to place files in

    character(len=*), parameter     :: grid_name = "gg"
    character(len=*), parameter     :: amr_name  = "mesh", meshdir = "data"
    character(len=100), allocatable :: grid_list(:), var_list(:, :), var_names(:)
    character(len=400)              :: fname
    integer                         :: lvl, i, id, i_grid, iv, nc, n_grids_max
    integer                         :: n_cc, dbix
    integer                         :: nx, ny, nx_prev, ny_prev, ix, iy
    integer                         :: n_cycle_val
    integer                         :: lo($D), hi($D), vlo($D), vhi($D)
    integer                         :: blo($D), bhi($D)
    logical                         :: lo_bnd($D), hi_bnd($D)
    integer, allocatable            :: ids(:), nb_ids(:), icc_val(:)
    logical, allocatable            :: box_done(:)
    real(dp)                        :: dr($D), r_min($D), time_val
#if $D == 2
    integer, allocatable            :: box_list(:,:), new_box_list(:, :)
    real(dp), allocatable           :: var_data(:,:,:)
#elif $D == 3
    integer, allocatable            :: box_list(:,:,:), new_box_list(:,:,:)
    real(dp), allocatable           :: var_data(:,:,:,:)
    integer                         :: nz, nz_prev, iz
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

    n_cc = size(icc_val)

    allocate(var_names(n_cc))
    var_names(1:n_cc) = tree%cc_names(icc_val)

    nc = tree%n_cell
    n_grids_max = 0
    do lvl = 1, tree%highest_lvl
       n_grids_max = n_grids_max + size(tree%lvls(lvl)%leaves)
    end do

    allocate(grid_list(n_grids_max))
    allocate(var_list(n_cc, n_grids_max))
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
    call SILO_set_time_varying(dbix)
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

          ! Check for (periodic) boundaries (this could give problems for
          ! complex geometries, e.g. a triangle block)
          id     = box_list(1, 1)
          lo_bnd = a$D_phys_boundary(tree%boxes, id, a$D_low_neighbs)

          id = box_list(nx, ny)
          hi_bnd = a$D_phys_boundary(tree%boxes, id, a$D_high_neighbs)

          lo(:) = 1
          where (.not. lo_bnd) lo = lo - 1

          hi = [nx, ny] * nc
          where (.not. hi_bnd) hi = hi + 1

          ! Include ghost cells around internal boundaries
          allocate(var_data(lo(1):hi(1), lo(2):hi(2), n_cc))

          do ix = 1, nx
             do iy = 1, ny
                id = box_list(ix, iy)

                ! Include ghost cells on internal block boundaries
                blo = 1
                where ([ix, iy] == 1 .and. .not. lo_bnd) blo = 0

                bhi = nc
                where ([ix, iy] == [nx, ny] .and. .not. hi_bnd) bhi = nc+1

                vlo = blo + ([ix, iy]-1) * nc
                vhi = bhi + ([ix, iy]-1) * nc

                var_data(vlo(1):vhi(1), vlo(2):vhi(2), 1:n_cc) = &
                     tree%boxes(id)%cc(blo(1):bhi(1), blo(2):bhi(2), icc_val)
             end do
          end do

          id = box_list(1, 1)
          dr = tree%boxes(id)%dr
          r_min = tree%boxes(id)%r_min - (1 - lo) * dr

          write(grid_list(i_grid), "(A,I0)") meshdir // '/' // grid_name, i_grid
          call SILO_add_grid(dbix, grid_list(i_grid), 2, &
               hi - lo + 2, r_min, dr, 1-lo, hi - [nx, ny] * nc)

          do iv = 1, n_cc
             write(var_list(iv, i_grid), "(A,I0)") meshdir // '/' // &
                  trim(var_names(iv)) // "_", i_grid
             call SILO_add_var(dbix, var_list(iv, i_grid), grid_list(i_grid), &
                  pack(var_data(:, :, iv), .true.), hi-lo+1)
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

          ! Check for (periodic) boundaries (this could give problems for
          ! complex geometries, e.g. a triangle block)
          id     = box_list(1, 1, 1)
          lo_bnd = a$D_phys_boundary(tree%boxes, id, a$D_low_neighbs)

          id = box_list(nx, ny, nz)
          hi_bnd = a$D_phys_boundary(tree%boxes, id, a$D_high_neighbs)

          lo(:) = 1
          where (.not. lo_bnd) lo = lo - 1

          hi = [nx, ny, nz] * nc
          where (.not. hi_bnd) hi = hi + 1

          ! Include ghost cells around internal boundaries
          allocate(var_data(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), n_cc))

          do iz = 1, nz
             do ix = 1, nx
                do iy = 1, ny
                   id     = box_list(ix, iy, iz)

                   ! Include ghost cells on internal block boundaries
                   blo = 1
                   where ([ix, iy, iz] == 1 .and. .not. lo_bnd) blo = 0

                   bhi = nc
                   where ([ix, iy, iz] == [nx, ny, nz] &
                        .and. .not. hi_bnd) bhi = nc+1

                   vlo = blo + ([ix, iy, iz]-1) * nc
                   vhi = bhi + ([ix, iy, iz]-1) * nc

                   var_data(vlo(1):vhi(1), vlo(2):vhi(2), vlo(3):vhi(3), 1:n_cc) = &
                        tree%boxes(id)%cc(blo(1):bhi(1), blo(2):bhi(2), &
                        blo(3):bhi(3), icc_val)
                end do
             end do
          end do

          id = box_list(1, 1, 1)
          dr = tree%boxes(id)%dr
          r_min = tree%boxes(id)%r_min - (1 - lo) * dr

          write(grid_list(i_grid), "(A,I0)") meshdir // '/' // grid_name, i_grid
          call SILO_add_grid(dbix, grid_list(i_grid), 3, &
               hi - lo + 2, r_min, dr, 1-lo, hi-[nx, ny, nz]*nc)

          do iv = 1, n_cc
             write(var_list(iv, i_grid), "(A,I0)") meshdir // '/' // &
                  trim(var_names(iv)) // "_", i_grid
             call SILO_add_var(dbix, var_list(iv, i_grid), grid_list(i_grid), &
                  pack(var_data(:, :, :, iv), .true.), hi-lo+1)
          end do

          deallocate(var_data)
          deallocate(box_list)
#endif

       end do
    end do

    call SILO_set_mmesh_grid(dbix, amr_name, grid_list(1:i_grid), &
         n_cycle_val, time_val)
    do iv = 1, n_cc
       call SILO_set_mmesh_var(dbix, trim(var_names(iv)), amr_name, &
            var_list(iv, 1:i_grid), n_cycle_val, time_val)
    end do
    call SILO_close_file(dbix)
    print *, "a$D_write_silo: written " // trim(fname)
  end subroutine a$D_write_silo

end module m_a$D_output
