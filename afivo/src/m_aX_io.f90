! This module contains routines for writing output files with Afivo
!
! Author: Jannis Teunissen
! License: GPLv3

module m_a$D_io

  use m_a$D_t

  implicit none
  private

  public :: a$D_write_vtk
  public :: a$D_write_silo

contains

  !> Write the cell centered data of a tree to a vtk unstructured file. Only the
  !> leaves of the tree are used
  subroutine a$D_write_vtk(tree, filename, cc_names, n_cycle, time, ixs_cc, &
       fc_names, ixs_fc, dir)
    use m_vtk
    use m_a$D_core, only: a$D_has_children

    type(a$D_t), intent(in)       :: tree        !< Tree to write out
    character(len=*), intent(in)  :: filename    !< Filename for the vtk file
    character(len=*), intent(in)  :: cc_names(:) !< Names of the cell-centered variables
    integer, intent(in)           :: n_cycle     !< Cycle-number for vtk file (counter)
    real(dp), intent(in)          :: time        !< Time for output file
    integer, intent(in), optional :: ixs_cc(:)   !< Oncly include these cell variables
    character(len=*), optional, intent(in) :: dir !< Directory to place files in
    !> If present, output fluxes with these names
    character(len=*), optional, intent(in) :: fc_names(:)
    integer, intent(in), optional :: ixs_fc(:)   !< Oncly include these face variables

    integer                       :: lvl, bc, bn, n, n_cells, n_nodes
    integer                       :: ig, i, j, id, n_ix, c_ix, n_grids
    integer                       :: cell_ix, node_ix
    integer                       :: n_cc, n_fc
    integer, parameter            :: n_ch = a$D_num_children
    integer                       :: nodes_per_box, cells_per_box
    real(dp), allocatable         :: coords(:), cc_vars(:,:)
    integer, allocatable          :: offsets(:), connects(:)
    integer, allocatable          :: cell_types(:), icc_used(:), ifc_used(:)
    type(vtk_t)                   :: vtkf
    character(len=400)            :: fname
    character(len=100), allocatable :: var_names(:)
#if $D == 3
    integer                       :: k, bn2
#endif

    if (.not. tree%ready) stop "Tree not ready"
    if (present(ixs_cc)) then
       if (maxval(ixs_cc) > tree%n_var_cell .or. &
            minval(ixs_cc) < 1) stop "a$D_write_vtk: wrong indices given (ixs_cc)"
       if (size(ixs_cc) /= size(cc_names)) &
            stop "a$D_write_vtk: size(cc_names) /= size(ixs_cc)"
       icc_used = ixs_cc
    else
       if (size(cc_names) /= tree%n_var_cell) &
            stop "a$D_write_vtk: size(cc_names) /= n_var_cell"
       icc_used = [(i, i = 1, tree%n_var_cell)]
    end if

    if (present(fc_names)) then
       if (.not. present(ixs_fc)) then
          stop "a$D_write_vtk: ixs_fc not present (but fc_names is)"
       else
          if (size(ixs_fc) * $D /= size(fc_names)) then
             stop "a$D_write_vtk: size(fc_names) /= size(ixs_fc) * $D"
          end if
       end if
       ifc_used = ixs_fc
    else
       allocate(ifc_used(0))
    end if

    n_cc = size(icc_used)
    n_fc = size(ifc_used)

    allocate(var_names(n_cc + n_fc * $D))
    var_names(1:n_cc) = cc_names
    if (present(fc_names)) var_names(n_cc+1:) = fc_names

    bc            = tree%n_cell     ! number of Box Cells
    bn            = tree%n_cell + 1 ! number of Box Nodes
    nodes_per_box = bn**$D
    cells_per_box = bc**$D

    n_grids = 0
    do lvl = 1, tree%max_lvl
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
    do lvl = 1, tree%max_lvl
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
                cc_vars(c_ix, 1:n_cc)     = tree%boxes(id)%cc(i, j, icc_used)
                cc_vars(c_ix, n_cc+1::$D) = &
                     0.5_dp * sum(tree%boxes(id)%fx(i:i+1, j, ifc_used))
                cc_vars(c_ix, n_cc+2::$D) = &
                     0.5_dp * sum(tree%boxes(id)%fy(i, j:j+1, ifc_used))
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
                   cc_vars(c_ix, 1:n_cc)     = tree%boxes(id)%cc(i, j, k, icc_used)
                   cc_vars(c_ix, n_cc+1::$D) = &
                        0.5_dp * sum(tree%boxes(id)%fx(i:i+1, j, k, ifc_used))
                   cc_vars(c_ix, n_cc+2::$D) = &
                        0.5_dp * sum(tree%boxes(id)%fy(i, j:j+1, k, ifc_used))
                   cc_vars(c_ix, n_cc+3::$D) = &
                        0.5_dp * sum(tree%boxes(id)%fz(i, j, k:k+1, ifc_used))
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
    call vtk_geo_xml(vtkf, coords, n_nodes, n_cells, $D, n_cycle, time)
    call vtk_con_xml(vtkf, connects, offsets, cell_types, n_cells)
    call vtk_dat_xml(vtkf, "CellData", .true.)

    do n = 1, n_cc + n_fc * $D
       call vtk_var_r8_xml(vtkf, trim(var_names(n)), cc_vars(:, n), n_cells)
    end do

    call vtk_dat_xml(vtkf, "CellData", .false.)
    call vtk_geo_xml_close(vtkf)
    call vtk_dat_xml(vtkf, "UnstructuredGrid", .false.)
    call vtk_end_xml(vtkf)
    print *, "Written ", trim(fname), ", n_grids", n_grids
  end subroutine a$D_write_vtk

  !> Write the cell centered data of a tree to a Silo file. Only the
  !> leaves of the tree are used
  subroutine a$D_write_silo(tree, filename, cc_names, n_cycle, time, ixs_cc, &
       fc_names, ixs_fc, dir)
    use m_write_silo
    use m_a$D_core, only: a$D_has_children

    type(a$D_t), intent(in)       :: tree        !< Tree to write out
    character(len=*)              :: filename    !< Filename for the vtk file
    character(len=*)              :: cc_names(:) !< Names of the cell-centered variables
    integer, intent(in)           :: n_cycle     !< Cycle-number for vtk file (counter)
    real(dp), intent(in)          :: time        !< Time for output file
    integer, intent(in), optional :: ixs_cc(:)      !< Oncly include these cell variables
    character(len=*), optional, intent(in) :: dir !< Directory to place files in
    !> If present, output fluxes with these names
    character(len=*), optional, intent(in) :: fc_names(:)
    integer, intent(in), optional :: ixs_fc(:)      !< Oncly include these face variables

    character(len=*), parameter     :: grid_name = "gg"
    character(len=*), parameter     :: amr_name  = "mesh", meshdir = "data"
    character(len=100), allocatable :: grid_list(:), var_list(:, :), var_names(:)
    character(len=400)              :: fname
    integer                         :: lvl, i, id, i_grid, iv, nc, n_grids_max
    integer                         :: n_vars, i0, j0, dbix, n_cc, n_fc
    integer                         :: nx, ny, nx_prev, ny_prev, ix, iy
    integer, allocatable            :: ids(:), nb_ids(:), icc_used(:), ifc_used(:)
    logical, allocatable            :: box_done(:)
    real(dp)                        :: dr($D), r_min($D)
#if $D == 2
    integer, allocatable            :: box_list(:,:), new_box_list(:, :)
    real(dp), allocatable           :: var_data(:,:,:)
#elif $D == 3
    integer, allocatable            :: box_list(:,:,:), new_box_list(:,:,:)
    real(dp), allocatable           :: var_data(:,:,:,:)
    integer                         :: k0, nz, nz_prev, iz
#endif

    if (present(ixs_cc)) then
       if (maxval(ixs_cc) > tree%n_var_cell .or. &
            minval(ixs_cc) < 1) stop "a$D_write_silo: wrong indices given (ixs_cc)"
       if (size(ixs_cc) /= size(cc_names)) &
            stop "a$D_write_silo: size(cc_names) /= size(ixs_cc)"
       icc_used = ixs_cc
    else
       if (size(cc_names) /= tree%n_var_cell) &
            stop "a$D_write_silo: size(cc_names) /= n_var_cell"
       icc_used = [(i, i = 1, tree%n_var_cell)]
    end if

    if (present(fc_names)) then
       if (.not. present(ixs_fc)) then
          stop "a$D_write_vtk: ixs_fc not present (but fc_names is)"
       else
          if (size(ixs_fc) * $D /= size(fc_names)) then
             stop "a$D_write_vtk: size(fc_names) /= size(ixs_fc) * $D"
          end if
       end if
       ifc_used = ixs_fc
    else
       allocate(ifc_used(0))
    end if

    n_cc = size(icc_used)
    n_fc = size(ifc_used)

    allocate(var_names(n_cc + n_fc * $D))
    var_names(1:n_cc) = cc_names
    if (present(fc_names)) var_names(n_cc+1:) = fc_names

    nc = tree%n_cell
    n_vars = n_cc + n_fc * $D
    n_grids_max = 0
    do lvl = 1, tree%max_lvl
       n_grids_max = n_grids_max + size(tree%lvls(lvl)%leaves)
    end do

    allocate(grid_list(n_grids_max))
    allocate(var_list(n_vars, n_grids_max))
    allocate(box_done(tree%max_id))
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

    do lvl = 1, tree%max_lvl
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
             nb_ids = tree%boxes(ids)%neighbors(a2_nb_lx)
             if (all(nb_ids > a5_no_box)) then
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
             nb_ids = tree%boxes(ids)%neighbors(a2_nb_hx)
             if (all(nb_ids > a5_no_box)) then
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
             nb_ids = tree%boxes(ids)%neighbors(a2_nb_ly)
             if (all(nb_ids > a5_no_box)) then
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
             nb_ids = tree%boxes(ids)%neighbors(a2_nb_hy)
             if (all(nb_ids > a5_no_box)) then
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
                     tree%boxes(id)%cc(1:nc, 1:nc, icc_used)
                var_data(i0:i0+nc-1, j0:j0+nc-1, n_cc+1::$D) = &
                     0.5_dp * (tree%boxes(id)%fx(1:nc, 1:nc, ifc_used) + &
                     tree%boxes(id)%fx(2:nc+1, 1:nc, ifc_used))
                var_data(i0:i0+nc-1, j0:j0+nc-1, n_cc+2::$D) = &
                     0.5_dp * (tree%boxes(id)%fy(1:nc, 1:nc, ifc_used) + &
                     tree%boxes(id)%fy(1:nc, 2:nc+1, ifc_used))
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
             nb_ids = tree%boxes(ids)%neighbors(a3_nb_lx)
             if (all(nb_ids > a5_no_box)) then
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
             nb_ids = tree%boxes(ids)%neighbors(a3_nb_hx)
             if (all(nb_ids > a5_no_box)) then
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
             nb_ids = tree%boxes(ids)%neighbors(a3_nb_ly)
             if (all(nb_ids > a5_no_box)) then
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
             nb_ids = tree%boxes(ids)%neighbors(a3_nb_hy)
             if (all(nb_ids > a5_no_box)) then
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
             nb_ids = tree%boxes(ids)%neighbors(a3_nb_lz)
             if (all(nb_ids > a5_no_box)) then
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
             nb_ids = tree%boxes(ids)%neighbors(a3_nb_hz)
             if (all(nb_ids > a5_no_box)) then
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
                        tree%boxes(id)%cc(1:nc, 1:nc, 1:nc, icc_used)
                   var_data(i0:i0+nc-1, j0:j0+nc-1, k0:k0+nc-1, n_cc+1::$D) = &
                        0.5_dp * (tree%boxes(id)%fx(1:nc, 1:nc, 1:nc, ifc_used) + &
                        tree%boxes(id)%fx(2:nc+1, 1:nc, 1:nc, ifc_used))
                   var_data(i0:i0+nc-1, j0:j0+nc-1, k0:k0+nc-1, n_cc+2::$D) = &
                        0.5_dp * (tree%boxes(id)%fy(1:nc, 1:nc, 1:nc, ifc_used) + &
                        tree%boxes(id)%fy(1:nc, 2:nc+1, 1:nc, ifc_used))
                   var_data(i0:i0+nc-1, j0:j0+nc-1, k0:k0+nc-1, n_cc+3::$D) = &
                        0.5_dp * (tree%boxes(id)%fz(1:nc, 1:nc, 1:nc, ifc_used) + &
                        tree%boxes(id)%fz(1:nc, 1:nc, 2:nc+1, ifc_used))
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

    call SILO_set_mmesh_grid(dbix, amr_name, grid_list(1:i_grid), n_cycle, time)
    do iv = 1, n_vars
       call SILO_set_mmesh_var(dbix, trim(var_names(iv)), amr_name, &
            var_list(iv, 1:i_grid), n_cycle, time)
    end do
    call SILO_close_file(dbix)

    print *, "Written ", trim(fname), ", n_grids", i_grid
  end subroutine a$D_write_silo

end module m_a$D_io
