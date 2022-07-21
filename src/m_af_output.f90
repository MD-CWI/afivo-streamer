#include "cpp_macros.h"
!> This module contains routines for writing output files with Afivo. The Silo
!> format should probably be used for larger files, especially in 3D.
module m_af_output

  use m_af_types

  implicit none
  private

  integer, parameter :: af_dat_file_version = 3

  abstract interface
     subroutine subr_add_vars(box, new_vars, n_var)
       import
       type(box_t), intent(in) :: box
       integer, intent(in)       :: n_var
       real(dp)                  :: &
            new_vars(DTIMES(0:box%n_cell+1), n_var)
     end subroutine subr_add_vars

     ! Subroutine that reads/writes unformatted data
     subroutine subr_other_data(my_unit)
       integer, intent(in) :: my_unit
     end subroutine subr_other_data
  end interface

  public :: af_write_tree
  public :: af_read_tree
  public :: af_tree_copy_variable
  public :: af_write_vtk
  public :: af_write_numpy
#if NDIM > 1
  public :: af_write_plane
#endif
  public :: af_write_silo
  public :: af_write_line

contains

  !> Write full tree in binary format
  subroutine af_write_tree(tree, filename, write_other_data)
    type(af_t), intent(in)                 :: tree     !< Tree to write out
    character(len=*), intent(in)           :: filename !< Filename for the output
    procedure(subr_other_data), optional   :: write_other_data
    integer                                :: my_unit, lvl, id, n
    character(len=400)                     :: fname

    fname = trim(filename) // ".dat"

    open(newunit=my_unit, file=trim(fname), form='unformatted', &
         access='stream', status='replace')

    write(my_unit) af_dat_file_version
    write(my_unit) tree%ready
    write(my_unit) tree%box_limit
    write(my_unit) tree%highest_lvl
    write(my_unit) tree%highest_id
    write(my_unit) tree%n_cell
    write(my_unit) tree%n_var_cell
    write(my_unit) tree%n_var_face
    write(my_unit) tree%coord_t
    write(my_unit) tree%coarse_grid_size
    write(my_unit) tree%periodic
    write(my_unit) tree%r_base
    write(my_unit) tree%dr_base

    write(my_unit) tree%cc_names
    write(my_unit) tree%fc_names
    write(my_unit) tree%cc_num_copies
    write(my_unit) tree%cc_write_output
    write(my_unit) tree%cc_write_binary
    write(my_unit) tree%fc_write_binary

    write(my_unit) tree%n_removed_ids
    write(my_unit) tree%removed_ids(1:tree%n_removed_ids)

    ! Skip methods (these have to be set again)

    do lvl = 1, tree%highest_lvl
       write(my_unit) size(tree%lvls(lvl)%ids)
       write(my_unit) tree%lvls(lvl)%ids

       write(my_unit) size(tree%lvls(lvl)%leaves)
       write(my_unit) tree%lvls(lvl)%leaves

       write(my_unit) size(tree%lvls(lvl)%parents)
       write(my_unit) tree%lvls(lvl)%parents
    end do

    do id = 1, tree%highest_id
       associate (box => tree%boxes(id))
         write(my_unit) box%in_use  !< is the box in use?
         if (.not. box%in_use) cycle

         write(my_unit) box%n_cell  !< number of cells per dimension
         write(my_unit) box%n_bc    !< Number of boundary conditions
         write(my_unit) box%n_stencils !< Number of stencils
         write(my_unit) box%lvl     !< level of the box
         write(my_unit) box%tag     !< for the user
         write(my_unit) box%ix      !< index in the domain
         write(my_unit) box%parent  !< index of parent in box list
         write(my_unit) box%children
         write(my_unit) box%neighbors
         write(my_unit) box%neighbor_mat
         write(my_unit) box%dr      !< width/height of a cell
         write(my_unit) box%r_min   !< min coords. of box
         write(my_unit) box%coord_t !< Coordinate type (e.g. Cartesian)

         do n = 1, tree%n_var_cell
            if (tree%cc_write_binary(n)) then
               write(my_unit) box%cc(DTIMES(:), n)
            end if
         end do

         do n = 1, tree%n_var_face
            if (tree%fc_write_binary(n)) then
               write(my_unit) box%fc(DTIMES(:), :, n)
            end if
         end do

         ! Write boundary conditions
         if (box%n_bc > 0) then
            write(my_unit) box%bc_index_to_nb
            write(my_unit) box%nb_to_bc_index
            write(my_unit) box%bc_type
            write(my_unit) box%bc_val
            write(my_unit) box%bc_coords
         end if

         ! Write stencils
         if (box%n_stencils > 0) then
            do n = 1, box%n_stencils
               associate (stncl => box%stencils(n))
                 write(my_unit) stncl%key
                 write(my_unit) stncl%shape
                 write(my_unit) stncl%stype
                 write(my_unit) stncl%cylindrical_gradient

                 if (allocated(stncl%c)) then
                    write(my_unit) size(stncl%c)
                    write(my_unit) stncl%c
                 else
                    write(my_unit) 0
                 end if

                 if (allocated(stncl%v)) then
                    write(my_unit) size(stncl%v, 1)
                    write(my_unit) stncl%v
                 else
                    write(my_unit) 0
                 end if

                 if (allocated(stncl%f)) then
                    write(my_unit) 1
                    write(my_unit) stncl%f
                 else
                    write(my_unit) 0
                 end if

                 if (allocated(stncl%bc_correction)) then
                    write(my_unit) 1
                    write(my_unit) stncl%bc_correction
                 else
                    write(my_unit) 0
                 end if

                 if (allocated(stncl%sparse_ix)) then
                    write(my_unit) size(stncl%sparse_ix, 2)
                    write(my_unit) stncl%sparse_ix
                 else
                    write(my_unit) 0
                 end if

                 if (allocated(stncl%sparse_v)) then
                    write(my_unit) size(stncl%sparse_v, 1)
                    write(my_unit) stncl%sparse_v
                 else
                    write(my_unit) 0
                 end if
               end associate
            end do
         end if
       end associate
    end do

    write(my_unit) present(write_other_data)

    if (present(write_other_data)) then
       call write_other_data(my_unit)
    end if

    close(my_unit)
    print *, "af_write_tree: written " // trim(fname)
  end subroutine af_write_tree

  !> Read full tree in binary format
  subroutine af_read_tree(tree, filename, read_other_data)
    use m_af_core, only: af_init_box
    type(af_t), intent(inout)    :: tree     !< Tree to read in
    character(len=*), intent(in) :: filename !< Filename for the input
    procedure(subr_other_data), optional :: read_other_data
    integer                      :: my_unit, lvl, n, id, version, k, m, nc
    logical                      :: other_data_present

    open(newunit=my_unit, file=trim(filename), form='unformatted', &
         access='stream', status='old', action='read')

    read(my_unit) version
    if (version /= af_dat_file_version) then
       print *, "File version read:", version
       print *, "File version required:", af_dat_file_version
       error stop "af_read_tree: incompatible file versions"
    end if

    read(my_unit) tree%ready
    read(my_unit) tree%box_limit
    read(my_unit) tree%highest_lvl
    read(my_unit) tree%highest_id
    read(my_unit) tree%n_cell
    read(my_unit) tree%n_var_cell
    read(my_unit) tree%n_var_face
    read(my_unit) tree%coord_t
    read(my_unit) tree%coarse_grid_size
    read(my_unit) tree%periodic
    read(my_unit) tree%r_base
    read(my_unit) tree%dr_base

    nc = tree%n_cell

    ! Skip methods (these have to be set again)
    if (.not. allocated(tree%cc_auto_vars)) &
         allocate(tree%cc_auto_vars(0))
    if (.not. allocated(tree%cc_func_vars)) &
         allocate(tree%cc_func_vars(0))

    read(my_unit) tree%cc_names
    read(my_unit) tree%fc_names
    read(my_unit) tree%cc_num_copies
    read(my_unit) tree%cc_write_output
    read(my_unit) tree%cc_write_binary
    read(my_unit) tree%fc_write_binary

    allocate(tree%removed_ids(tree%box_limit))
    read(my_unit) tree%n_removed_ids
    read(my_unit) tree%removed_ids(1:tree%n_removed_ids)

    do lvl = 1, tree%highest_lvl
       read(my_unit) n
       allocate(tree%lvls(lvl)%ids(n))
       read(my_unit) tree%lvls(lvl)%ids

       read(my_unit) n
       allocate(tree%lvls(lvl)%leaves(n))
       read(my_unit) tree%lvls(lvl)%leaves

       read(my_unit) n
       allocate(tree%lvls(lvl)%parents(n))
       read(my_unit) tree%lvls(lvl)%parents
    end do

    do lvl = tree%highest_lvl+1, af_max_lvl
       allocate(tree%lvls(lvl)%ids(0))
       allocate(tree%lvls(lvl)%leaves(0))
       allocate(tree%lvls(lvl)%parents(0))
    end do

    allocate(tree%boxes(tree%box_limit))

    do id = 1, tree%highest_id
       associate (box => tree%boxes(id))
         read(my_unit) box%in_use  !< is the box in use?
         if (.not. box%in_use) cycle

         read(my_unit) box%n_cell  !< number of cells per dimension
         read(my_unit) box%n_bc    !< Number of boundary conditions
         read(my_unit) box%n_stencils !< Number of stencils
         read(my_unit) box%lvl     !< level of the box
         read(my_unit) box%tag     !< for the user
         read(my_unit) box%ix      !< index in the domain
         read(my_unit) box%parent  !< index of parent in box list
         read(my_unit) box%children
         read(my_unit) box%neighbors
         read(my_unit) box%neighbor_mat
         read(my_unit) box%dr      !< width/height of a cell
         read(my_unit) box%r_min   !< min coords. of box
         read(my_unit) box%coord_t !< Coordinate type (e.g. Cartesian)

         call af_init_box(tree, id)

         do n = 1, tree%n_var_cell
            if (tree%cc_write_binary(n)) then
               read(my_unit) box%cc(DTIMES(:), n)
            end if
         end do

         do n = 1, tree%n_var_face
            if (tree%fc_write_binary(n)) then
               read(my_unit) box%fc(DTIMES(:), :, n)
            end if
         end do

         ! Read boundary conditions
         if (box%n_bc > 0) then
            read(my_unit) box%bc_index_to_nb
            read(my_unit) box%nb_to_bc_index
            read(my_unit) box%bc_type
            read(my_unit) box%bc_val
            read(my_unit) box%bc_coords
         end if

         ! Read stencils
         if (box%n_stencils > 0) then
            allocate(box%stencils(box%n_stencils))

            do n = 1, box%n_stencils
               associate (stncl => box%stencils(n))
                 read(my_unit) stncl%key
                 read(my_unit) stncl%shape
                 read(my_unit) stncl%stype
                 read(my_unit) stncl%cylindrical_gradient

                 read(my_unit) k
                 if (k > 0) then
                    allocate(stncl%c(k))
                    read(my_unit) stncl%c
                 end if

                 read(my_unit) k
                 if (k > 0) then
                    allocate(stncl%v(k, DTIMES(nc)))
                    read(my_unit) stncl%v
                 end if

                 read(my_unit) k
                 if (k > 0) then
                    allocate(stncl%f(DTIMES(nc)))
                    read(my_unit) stncl%f
                 end if

                 read(my_unit) k
                 if (k > 0) then
                    allocate(stncl%bc_correction(DTIMES(nc)))
                    read(my_unit) stncl%bc_correction
                 end if

                 read(my_unit) k
                 if (k > 0) then
                    allocate(stncl%sparse_ix(NDIM, k))
                    read(my_unit) stncl%sparse_ix
                 end if

                 read(my_unit) m
                 if (k > 0 .and. m > 0) then
                    allocate(stncl%sparse_v(m, k))
                    read(my_unit) stncl%sparse_v
                 end if
               end associate
            end do
         end if
       end associate
    end do

    read(my_unit) other_data_present

    if (present(read_other_data)) then
       if (other_data_present) then
          call read_other_data(my_unit)
       else
          error stop "af_read_tree: other data is not present"
       end if
    end if

    close(my_unit)
  end subroutine af_read_tree

  subroutine af_tree_copy_variable(tree_from, ivs_from, tree_to, ivs_to)
    use m_af_interp
    type(af_t), intent(in)    :: tree_from   !< Copy from this grid
    integer, intent(in)       :: ivs_from(:) !< From these variable
    type(af_t), intent(inout) :: tree_to     !< Copy to this grid
    integer, intent(in)       :: ivs_to(:)   !< To these variable
    integer                   :: lvl, id, n, nc, i, j, k
    real(dp)                  :: rr(NDIM)
    logical                   :: success

    !$omp parallel private(lvl, id, n, nc, i, j, k, rr)
    do lvl = 1, tree_to%highest_lvl
       !$omp do
       do n = 1, size(tree_to%lvls(lvl)%leaves)
          id = tree_to%lvls(lvl)%leaves(n)
          nc = tree_to%boxes(id)%n_cell

          do KJI_DO(1, nc)
             rr = af_r_cc(tree_to%boxes(id), [IJK])
             tree_to%boxes(id)%cc(IJK, ivs_to) = &
                  af_interp1(tree_from, rr, [ivs_from], success)
             if (.not. success) error stop "af_tree_copy_variable: interpolation error"
          end do; CLOSE_DO
       end do
       !$omp end do
    end do
    !$omp end parallel

  end subroutine af_tree_copy_variable

  !> Write line data in a text file
  subroutine af_write_line(tree, filename, ivs, r_min, r_max, n_points)
    use m_af_interp, only: af_interp1
    type(af_t), intent(in)       :: tree        !< Tree to write out
    character(len=*), intent(in) :: filename    !< Filename for the vtk file
    integer, intent(in)          :: ivs(:)      !< Variables to write
    real(dp), intent(in)         :: r_min(NDIM)   !< Minimum coordinate of line
    real(dp), intent(in)         :: r_max(NDIM)   !< Maximum coordinate of line
    integer, intent(in)          :: n_points !< Number of points along line

    integer, parameter    :: my_unit = 100
    character(len=400)    :: fname
    integer               :: i, n_cc
    real(dp)              :: r(NDIM), dr_vec(NDIM)
    real(dp), allocatable :: line_data(:, :)
    logical               :: success

    n_cc = size(ivs)
    dr_vec = (r_max - r_min) / max(1, n_points-1)

    allocate(line_data(n_cc+NDIM, n_points))

    !$omp parallel do private(r)
    do i = 1, n_points
       r = r_min + (i-1) * dr_vec
       line_data(1:NDIM, i) = r
       line_data(NDIM+1:NDIM+n_cc, i) = af_interp1(tree, r, ivs, success)
       if (.not. success) error stop "af_write_line: interpolation error"
    end do
    !$omp end parallel do

    fname = trim(filename) // ".txt"

    ! Write header
    open(my_unit, file=trim(fname), action="write")
#if NDIM == 1
    write(my_unit, '(A)', advance="no") "# x"
#elif NDIM == 2
    write(my_unit, '(A)', advance="no") "# x y"
#elif NDIM == 3
    write(my_unit, '(A)', advance="no") "# x y z"
#endif
    do i = 1, n_cc
       write(my_unit, '(A)', advance="no") " "//trim(tree%cc_names(ivs(i)))
    end do
    write(my_unit, '(A)') ""

    ! Write data
    do i = 1, n_points
       write(my_unit, *) line_data(:, i)
    end do

    close(my_unit)
  end subroutine af_write_line

#if NDIM > 1
  !> Write data in a plane (2D) to a VTK ASCII file. In 3D, r_min and r_max
  !> should have one identical coordinate (i.e., they differ in two
  !> coordinates).
  subroutine af_write_plane(tree, filename, ivs, r_min, r_max, n_pixels)
    use m_af_interp, only: af_interp1
    type(af_t), intent(in)       :: tree        !< Tree to write out
    character(len=*), intent(in) :: filename    !< Filename for the vtk file
    integer, intent(in)          :: ivs(:)      !< Variables to write
    real(dp), intent(in)         :: r_min(NDIM)   !< Minimum coordinate of plane
    real(dp), intent(in)         :: r_max(NDIM)   !< Maximum coordinate of plane
    integer, intent(in)          :: n_pixels(2) !< Number of pixels in the plane

    integer, parameter    :: my_unit = 100
    character(len=100)    :: fmt_string
    character(len=400)    :: fname
    integer               :: i, j, n_cc, dim_unused, n_points(3)
    real(dp)              :: r(NDIM), dvec(3)
    real(dp)              :: v1(NDIM), v2(NDIM)
    real(dp), allocatable :: pixel_data(:, :, :)
    logical               :: success

    n_cc      = size(ivs)
#if NDIM == 2
    dvec(1:2)   = r_max(1:2) - r_min(1:2)
    dvec(3)     = 0
    dim_unused  = 3
    n_points    = [n_pixels(1), n_pixels(2), 1]
    v1          = [dvec(1), 0.0_dp] / (n_pixels(1) - 1)
    v2          = [0.0_dp, dvec(2)] / (n_pixels(2) - 1)
#elif NDIM == 3
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

    !$omp parallel do private(i, r)
    do j = 1, n_pixels(2)
       do i = 1, n_pixels(1)
          r = r_min + (i-1) * v1 + (j-1) * v2
          pixel_data(:, i, j) = af_interp1(tree, r, ivs, success)
          if (.not. success) error stop "af_write_plane: interpolation error"
       end do
    end do
    !$omp end parallel do

    ! Construct format string. Write one row at a time
    write(fmt_string, '(A,I0,A)') '(', n_pixels(1), 'E20.8)'

    ! Construct file name
    fname = trim(filename) // ".vtk"

    open(my_unit, file=trim(fname), action="write")
    write(my_unit, '(A)') "# vtk DataFile Version 2.0"
    write(my_unit, '(A)') trim(filename)
    write(my_unit, '(A)') "ASCII"
    write(my_unit, '(A)') "DATASET STRUCTURED_POINTS"
    write(my_unit, '(A,3I10)') "DIMENSIONS ", n_points
#if NDIM == 2
    write(my_unit, '(A,3E20.8)') "ORIGIN ", [r_min(1), r_min(2), 0.0_dp]
    write(my_unit, '(A,3E20.8)') "SPACING ", &
         [v1(1) + v2(1), v1(2) + v2(2), 0.0_dp]
#elif NDIM == 3
    write(my_unit, '(A,3E20.8)') "ORIGIN ", r_min
    write(my_unit, '(A,3E20.8)') "SPACING ", v1 + v2
#endif
    write(my_unit, '(A,2I0)') "POINT_DATA ", product(n_points)
    do i = 1, n_cc
       write(my_unit, '(A)') "SCALARS " // &
            trim(tree%cc_names(ivs(i))) // " double 1"
       write(my_unit, '(A)') "LOOKUP_TABLE default"
       write(my_unit, trim(fmt_string)) pixel_data(i, :, :)
    end do
    close(my_unit)
  end subroutine af_write_plane
#endif

  !> Write the cell centered data of a tree to a vtk unstructured file. Only the
  !> leaves of the tree are used
  subroutine af_write_vtk(tree, filename, n_cycle, time, ixs_cc, &
       add_vars, add_names)
    use m_vtk

    type(af_t), intent(in)       :: tree        !< Tree to write out
    character(len=*), intent(in)  :: filename    !< Filename for the vtk file
    integer, intent(in), optional :: n_cycle     !< Cycle-number for vtk file (counter)
    real(dp), intent(in), optional :: time        !< Time for output file
    integer, intent(in), optional :: ixs_cc(:)   !< Only include these cell variables
    procedure(subr_add_vars), optional :: add_vars !< Optional routine to add extra variables
    character(len=*), intent(in), optional :: add_names(:) !< Names of extra variables

    integer                       :: lvl, bc, bn, n, n_cells, n_nodes
    integer                       :: ig, IJK, id, n_ix, c_ix, n_grids
    integer                       :: cell_ix, node_ix, n_cycle_val
    integer                       :: n_cc, n_add
    integer, parameter            :: n_ch = af_num_children
    integer                       :: nodes_per_box, cells_per_box
    real(dp)                      :: time_val
    real(dp), allocatable         :: coords(:), cc_vars(:,:)
    integer, allocatable          :: offsets(:), connects(:)
    integer, allocatable          :: cell_types(:), icc_val(:)
    type(vtk_t)                   :: vtkf
    character(len=400)            :: fname
    character(len=100), allocatable :: var_names(:)
    real(dp), allocatable         :: cc(DTIMES(:), :)
#if NDIM == 3
    integer                       :: bn2
#endif

    if (.not. tree%ready) error stop "Tree not ready"
    time_val = 0.0_dp; if (present(time)) time_val = time
    n_cycle_val = 0; if (present(n_cycle)) n_cycle_val = n_cycle
    n_add = 0; if (present(add_names)) n_add = size(add_names)

    if (present(add_names) .neqv. present(add_vars)) &
         stop "af_write_vtk: both arguments (add_names, add_vars) needed"

    if (present(ixs_cc)) then
       if (maxval(ixs_cc) > tree%n_var_cell .or. &
            minval(ixs_cc) < 1) stop "af_write_vtk: wrong indices given (ixs_cc)"
       allocate(icc_val(size(ixs_cc)))
       icc_val = ixs_cc
    else
       call get_output_vars(tree, icc_val)
    end if

    n_cc               = size(icc_val)

    allocate(var_names(n_cc+n_add))
    var_names(1:n_cc) = tree%cc_names(icc_val)

    if (present(add_names)) then
       var_names(n_cc+1:n_cc+n_add) = add_names(:)
    end if

    bc            = tree%n_cell     ! number of Box Cells
    bn            = tree%n_cell + 1 ! number of Box Nodes
    nodes_per_box = bn**NDIM
    cells_per_box = bc**NDIM

    allocate(cc(DTIMES(0:bc+1), n_cc + n_add))

    n_grids = 0
    do lvl = 1, tree%highest_lvl
       n_grids = n_grids + size(tree%lvls(lvl)%leaves)
    end do
    n_nodes = nodes_per_box * n_grids
    n_cells = cells_per_box * n_grids

    allocate(coords(NDIM * n_nodes))
    allocate(cc_vars(n_cells, n_cc+n_add))
    allocate(offsets(cells_per_box * n_grids))
    allocate(cell_types(cells_per_box * n_grids))
    allocate(connects(n_ch * cells_per_box * n_grids))

#if NDIM == 1
    cell_types = 3  ! VTK line type
#elif NDIM == 2
    cell_types = 8  ! VTK pixel type
#elif NDIM       == 3
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

#if NDIM == 1
          cc(:, 1:n_cc) = tree%boxes(id)%cc(:, icc_val)

          if (present(add_vars)) then
             call add_vars(tree%boxes(id), &
                  cc(:, n_cc+1:n_cc+n_add), n_add)
          end if

          do i = 1, bn
             n_ix = node_ix + i
             coords(n_ix) = tree%boxes(id)%r_min(1) + &
                  (i-1) * tree%boxes(id)%dr(1)
          end do

          do i = 1, bc
             ! In vtk, indexing starts at 0, so subtract 1
             n_ix                      = node_ix + i - 1
             c_ix                      = cell_ix + i
             cc_vars(c_ix, :)          = cc(i, :)
             offsets(c_ix)             = af_num_children * c_ix
             connects(n_ch*(c_ix-1)+1:n_ch*c_ix) = [n_ix, n_ix+1]
          end do
#elif NDIM == 2
          cc(:, :, 1:n_cc) = tree%boxes(id)%cc(:, :, icc_val)

          if (present(add_vars)) then
             call add_vars(tree%boxes(id), &
                  cc(:, :, n_cc+1:n_cc+n_add), n_add)
          end if

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
                cc_vars(c_ix, :)          = cc(i, j, :)
                offsets(c_ix)             = af_num_children * c_ix
                connects(n_ch*(c_ix-1)+1:n_ch*c_ix) = [n_ix, n_ix+1, n_ix+bn, n_ix+bn+1]
             end do
          end do
#elif NDIM == 3
          cc(:, :, :, 1:n_cc) = tree%boxes(id)%cc(:, :, :, icc_val)

          if (present(add_vars)) then
             call add_vars(tree%boxes(id), &
                  cc(:, :, :, n_cc+1:n_cc+n_add), n_add)
          end if

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
                   cc_vars(c_ix, :)          = cc(i, j, k, :)
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

    call vtk_ini_xml(vtkf, trim(fname), 'UnstructuredGrid')
    call vtk_dat_xml(vtkf, "UnstructuredGrid", .true.)
    call vtk_unstr_geo_xml(vtkf, coords, n_nodes, n_cells, NDIM, n_cycle_val, time_val)
    call vtk_unstr_con_xml(vtkf, connects, offsets, cell_types, n_cells)
    call vtk_dat_xml(vtkf, "CellData", .true.)

    do n = 1, n_cc + n_add
       call vtk_var_r8_xml(vtkf, trim(var_names(n)), cc_vars(:, n), n_cells)
    end do

    call vtk_dat_xml(vtkf, "CellData", .false.)
    call vtk_unstr_geo_xml_close(vtkf)
    call vtk_dat_xml(vtkf, "UnstructuredGrid", .false.)
    call vtk_end_xml(vtkf)
    print *, "af_write_vtk: written " // trim(fname)
  end subroutine af_write_vtk

  !> Write uniform data interpolated from a region to a .npy or .npz numpy file.
  !> The format is determined based on the extension of filename
  subroutine af_write_numpy(tree, filename, r_min, r_max, n_points, &
       n_cycle, time, ixs_cc)
    use m_npy
    use m_af_interp, only: af_interp1

    type(af_t), intent(inout)      :: tree           !< Tree to save
    !> Filename, possible extensions: .npz, .npy
    character(len=*), intent(in)   :: filename
    integer, intent(in), optional  :: n_points(NDIM) !< Number of points to use
    real(dp), intent(in), optional :: r_min(NDIM)    !< Minimum coordinates
    real(dp), intent(in), optional :: r_max(NDIM)    !< Maximum coordinates
    integer, intent(in), optional  :: n_cycle        !< Cycle-number (counter)
    real(dp), intent(in), optional :: time           !< Time
    integer, intent(in), optional  :: ixs_cc(:)      !< Only include these cell variables

    integer                             :: n, IJK, id, lvl, id_guess, nx(NDIM)
    integer                             :: n_cycle_val, n_cc, nc
    logical                             :: success
    real(dp)                            :: time_val, dr(NDIM), r(NDIM)
    real(dp)                            :: rmin(NDIM), rmax(NDIM)
    integer, allocatable                :: icc_val(:)
    character(len=af_nlen), allocatable :: var_names(:)
    character(len=400)                  :: fname
    real(dp), allocatable               :: cc(DTIMES(:), :)

    if (.not. tree%ready) error stop "Tree not ready"
    nc = tree%n_cell
    time_val = 0.0_dp; if (present(time)) time_val = time
    n_cycle_val = 0; if (present(n_cycle)) n_cycle_val = n_cycle
    rmin = tree%r_base
    if (present(r_min)) rmin = r_min
    rmax = tree%r_base + tree%dr_base * tree%coarse_grid_size
    if (present(r_max)) rmax = r_max

    if (present(ixs_cc)) then
       if (maxval(ixs_cc) > tree%n_var_cell .or. &
            minval(ixs_cc) < 1) stop "af_write_vtk: wrong indices given (ixs_cc)"
       allocate(icc_val(size(ixs_cc)))
       icc_val = ixs_cc
    else
       call get_output_vars(tree, icc_val)
    end if
    n_cc = size(icc_val)

    if (present(n_points)) then
       nx = n_points
    else
       ! Determine highest fully refined grid level
       outer: do lvl = 1, tree%highest_lvl-1
          do i = 1, size(tree%lvls(lvl)%leaves)
             id = tree%lvls(lvl)%leaves(i)

             ! If there is a leave in the region, the highest is found
             if (all(tree%boxes(id)%r_min <= rmax .and. &
                  tree%boxes(id)%r_min + nc * &
                  tree%boxes(id)%dr >= rmin)) exit outer
          end do
       end do outer

       ! Use same grid spacing as the grid
       nx = nint((rmax - rmin)/af_lvl_dr(tree, lvl))
    end if

    dr = (rmax - rmin)/nx

    allocate(var_names(n_cc))
    var_names(1:n_cc) = tree%cc_names(icc_val)
    allocate(cc(DINDEX(nx), n_cc))

    id_guess = -1
    !$omp parallel do private(IJK, r) firstprivate(id_guess)
#if NDIM == 1
    do i = 1, nx(1)
       r = rmin + ([IJK] - 0.5_dp) * dr
       cc(IJK, :) = af_interp1(tree, r, icc_val, success, id_guess)
       if (.not. success) error stop "af_write_plane: interpolation error"
    end do
#elif NDIM == 2
    do j = 1, nx(2)
       do i = 1, nx(1)
          r = rmin + ([IJK] - 0.5) * dr
          cc(IJK, :) = af_interp1(tree, r, icc_val, success, id_guess)
          if (.not. success) error stop "af_write_plane: interpolation error"
       end do
    end do
#elif NDIM == 3
    do k = 1, nx(3)
       do j = 1, nx(2)
          do i = 1, nx(1)
             r = rmin + ([IJK] - 0.5) * dr
             cc(IJK, :) = af_interp1(tree, r, icc_val, success, id_guess)
             if (.not. success) error stop "af_write_plane: interpolation error"
          end do
       end do
    end do
#endif
    !$omp end parallel do

    n = len_trim(filename)

    ! Check last four characters of filename
    select case (filename(n-3:n))
    case ('.npy')
       call save_npy(filename, cc)
    case ('.npz')
       call remove_file(filename) ! Clear npz file

       ! Add variables separately
       do n = 1, n_cc
          fname = filename(:n-4)//trim(var_names(n))//'.npy'
          call save_npy(fname, cc(DTIMES(:), n))
          call add_to_zip(filename, fname, .false., var_names(n))
       end do

       fname = filename(:n-4)//'nx.npy'
       call save_npy(fname, nx)
       call add_to_zip(filename, fname, .false., 'nx')

       fname = filename(:n-4)//'r_min.npy'
       call save_npy(fname, rmin)
       call add_to_zip(filename, fname, .false., 'r_min')

       fname = filename(:n-4)//'r_max.npy'
       call save_npy(fname, rmax)
       call add_to_zip(filename, fname, .false., 'r_max')

       fname = filename(:n-4)//'dr.npy'
       call save_npy(fname, dr)
       call add_to_zip(filename, fname, .false., 'dr')

       fname = filename(:n-4)//'coord_t.npy'
       call save_npy(fname, [tree%coord_t])
       call add_to_zip(filename, fname, .false., 'coord_t')

       fname = filename(:n-4)//'time_cycle.npy'
       call save_npy(fname, [time_val, real(n_cycle_val, dp)])
       call add_to_zip(filename, fname, .false., 'time_cycle')
    case default
       error stop "Unknown file extension"
    end select

    print *, "af_write_numpy: written " // trim(filename)
  end subroutine af_write_numpy

#if NDIM == 1
  subroutine af_write_silo(tree, filename, n_cycle, time, ixs_cc, &
       add_vars, add_names)
    use m_write_silo
    use m_mrgrnk

    type(af_t), intent(in)                 :: tree         !< Tree to write out
    character(len=*)                       :: filename     !< Filename for the vtk file
    integer, intent(in), optional          :: n_cycle      !< Cycle-number for vtk file (counter)
    real(dp), intent(in), optional         :: time         !< Time for output file
    integer, intent(in), optional          :: ixs_cc(:)    !< Only include these cell variables
    procedure(subr_add_vars), optional     :: add_vars     !< Optional routine to add extra variables
    character(len=*), intent(in), optional :: add_names(:) !< Names of extra variables
    integer                                :: n, i, id, ix, lvl, n_boxes
    integer                                :: dbix, n_cc, n_leaves, nc
    integer, allocatable                   :: icc_val(:), id_leaves(:)
    integer, allocatable                   :: id_sorted_ix(:)
    real(dp), allocatable                  :: xmin_leaves(:)
    real(dp), allocatable                  :: xdata(:), ydata(:, :)
    character(len=af_nlen)                 :: yname
    character(len=400)                     :: fname

    if (present(ixs_cc)) then
       if (maxval(ixs_cc) > tree%n_var_cell .or. &
            minval(ixs_cc) < 1) stop "af_write_silo: wrong indices given (ixs_cc)"
       allocate(icc_val(size(ixs_cc)))
       icc_val = ixs_cc
    else
       call get_output_vars(tree, icc_val)
    end if

    if (present(add_vars) .or. present(add_names)) &
         error stop "add_vars / add_names not implemented in 1D"

    n_cc     = size(icc_val)
    n_leaves = af_num_leaves_used(tree)
    nc       = tree%n_cell

    allocate(ydata(n_leaves * nc, n_cc))
    allocate(xdata(n_leaves * nc))
    allocate(id_leaves(n_leaves))
    allocate(xmin_leaves(n_leaves))
    allocate(id_sorted_ix(n_leaves))

    ! Store ids of all leaves
    n = 0
    do lvl = 1, tree%highest_lvl
       n_boxes = size(tree%lvls(lvl)%leaves)
       id_leaves(n+1:n+n_boxes) = tree%lvls(lvl)%leaves
       n = n + n_boxes
    end do

    ! Get minimum coordinate of each leave
    xmin_leaves = tree%boxes(id_leaves)%r_min(1)

    ! Sort by xmin
    call mrgrnk(xmin_leaves, id_sorted_ix)

    do n = 1, n_leaves
       ix = (n-1) * nc
       id = id_leaves(id_sorted_ix(n))
       ydata(ix+1:ix+nc, :) = tree%boxes(id)%cc(1:nc, icc_val)
       do i = 1, nc
          xdata(ix+i:ix+i) = af_r_cc(tree%boxes(id), [i])
       end do
    end do

    fname = trim(filename) // ".silo"
    call SILO_create_file(trim(fname), dbix)
    call SILO_set_time_varying(dbix)

    do n = 1, n_cc
       yname = tree%cc_names(icc_val(n))
       call SILO_add_curve(dbix, yname, xdata, ydata(:, n))
    end do

    call SILO_close_file(dbix)
    print *, "af_write_silo: written " // trim(fname)
  end subroutine af_write_silo
#else
  !> Write the cell centered data of a tree to a Silo file. Only the
  !> leaves of the tree are used
  !>
  !> Note: a 1D version is present below, and seems to work, but its output
  !> cannot be visualized by Visit. That's why we use curve output for 1D cases.
  subroutine af_write_silo(tree, filename, n_cycle, time, ixs_cc, &
       add_vars, add_names, add_curve_names, add_curve_dat)
    use m_write_silo

    type(af_t), intent(in)       :: tree        !< Tree to write out
    character(len=*)              :: filename    !< Filename for the vtk file
    integer, intent(in), optional :: n_cycle     !< Cycle-number for vtk file (counter)
    real(dp), intent(in), optional :: time        !< Time for output file
    integer, intent(in), optional :: ixs_cc(:)      !< Only include these cell variables
    procedure(subr_add_vars), optional :: add_vars !< Optional routine to add extra variables
    character(len=*), intent(in), optional :: add_names(:), add_curve_names(:) !< Names of extra variables or curves
    real(dp), intent(in), optional  :: add_curve_dat(:, :, :) !< Data for additional curves (#curves, x/y-dimensions, #points)

    character(len=*), parameter     :: grid_name = "gg", block_prefix = "blk_"
    character(len=*), parameter     :: amr_name  = "mesh", meshdir = "data"
    character(len=100), allocatable :: grid_list(:), grid_list_block(:)
    character(len=100), allocatable :: var_list(:, :), var_names(:)
    character(len=400)              :: fname
    integer                         :: lvl, i, id, i_grid, iv, nc, n_grids_max
    integer                         :: n_cc, n_add, dbix
    integer                         :: nx, nx_prev, ix
    integer                         :: n_cycle_val, ncurve, ncurve_add
    integer                         :: lo(NDIM), hi(NDIM), vlo(NDIM), vhi(NDIM)
    integer                         :: blo(NDIM), bhi(NDIM)
    logical                         :: lo_bnd(NDIM), hi_bnd(NDIM)
    integer, allocatable            :: ids(:), nb_ids(:), icc_val(:)
    logical, allocatable            :: box_done(:)
    real(dp)                        :: dr(NDIM), r_min(NDIM), time_val
    integer, allocatable            :: box_list(DTIMES(:)), new_box_list(DTIMES(:))
    real(dp), allocatable           :: var_data(DTIMES(:),:), cc(DTIMES(:), :)
#if NDIM > 1
    integer                         :: ny, ny_prev, iy
#endif
#if NDIM > 2
    integer                         :: nz, nz_prev, iz
#endif

    if (.not. tree%ready) stop "Tree not ready"
    time_val = 0.0_dp; if (present(time)) time_val = time
    n_cycle_val = 0; if (present(n_cycle)) n_cycle_val = n_cycle
    n_add = 0; if (present(add_names)) n_add = size(add_names)

    if (present(add_names) .neqv. present(add_vars)) &
         stop "af_write_silo: both arguments (add_names, add_vars) needed"

    if (present(add_curve_names) .neqv. present(add_curve_dat)) &
         stop "af_write_silo: both arguments (add_curve_names, add_curve_dat) needed"

    if (present(add_curve_names)) then
       if (size(add_curve_names, 1) .ne. size(add_curve_dat, 1)) &
            stop "af_write_silo: number of curve names and data do not agree"
    end if

    if (present(ixs_cc)) then
       if (maxval(ixs_cc) > tree%n_var_cell .or. &
            minval(ixs_cc) < 1) stop "af_write_silo: wrong indices given (ixs_cc)"
       allocate(icc_val(size(ixs_cc)))
       icc_val = ixs_cc
    else
       call get_output_vars(tree, icc_val)
    end if

    n_cc = size(icc_val)

    allocate(var_names(n_cc+n_add))
    var_names(1:n_cc) = tree%cc_names(icc_val)

    if (present(add_names)) then
       var_names(n_cc+1:n_cc+n_add) = add_names(:)
    end if

    nc = tree%n_cell
    n_grids_max = 0
    do lvl = 1, tree%highest_lvl
       n_grids_max = n_grids_max + size(tree%lvls(lvl)%leaves)
    end do

    allocate(grid_list(n_grids_max))
    allocate(grid_list_block(n_grids_max))
    allocate(var_list(n_cc+n_add, n_grids_max))
    allocate(box_done(tree%highest_id))
    box_done = .false.

    allocate(cc(DTIMES(0:nc+1), n_cc + n_add))

    fname = trim(filename) // ".silo"
    call SILO_create_file(trim(fname), dbix)
    call SILO_set_time_varying(dbix)
    call SILO_mkdir(dbix, meshdir)

    ! Adding additional curve objects
    if (present(add_curve_names)) then
      ncurve_add = size(add_curve_dat, 1)
      do ncurve = 1, ncurve_add
        call SILO_add_curve(dbix, add_curve_names(ncurve), &
          add_curve_dat(ncurve, 1, :), add_curve_dat(ncurve, 2, :))
      end do
    end if

    i_grid = 0

    do lvl = 1, tree%highest_lvl
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          if (box_done(id)) cycle

          i_grid = i_grid + 1

          ! Find largest rectangular box including id and other leaves that
          ! haven't been written yet
#if NDIM == 1
          allocate(box_list(1))
          box_list(1) = id
          box_done(id) = .true.
          nx = 1

          do
             nx_prev = nx

             ! Check whether we can extend to the -x direction
             ids = box_list(1:1)
             nb_ids = tree%boxes(ids)%neighbors(af_neighb_lowx)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(1) < tree%boxes(ids(1))%ix(1) .and. &
                     .not. any(af_has_children(tree%boxes(nb_ids)))) then
                   nx = nx + 1
                   allocate(new_box_list(nx))
                   new_box_list(1:1) = nb_ids
                   new_box_list(2:) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             ! Check whether we can extend to the +x direction
             ids = box_list(nx:nx)
             nb_ids = tree%boxes(ids)%neighbors(af_neighb_highx)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(1) > tree%boxes(ids(1))%ix(1) .and. &
                     .not. any(af_has_children(tree%boxes(nb_ids)))) then
                   nx = nx + 1
                   allocate(new_box_list(nx))
                   new_box_list(nx:nx) = nb_ids
                   new_box_list(1:nx-1) = box_list
                   box_list = new_box_list
                   box_done(nb_ids) = .true.
                   deallocate(new_box_list)
                end if
             end if

             if (nx == nx_prev) exit
          end do

          ! Check for (periodic) boundaries (this could give problems for
          ! complex geometries, e.g. a triangle block)
          id     = box_list(1)
          lo_bnd = af_is_phys_boundary(tree%boxes, id, af_low_neighbs)

          id = box_list(nx)
          hi_bnd = af_is_phys_boundary(tree%boxes, id, af_high_neighbs)

          lo(:) = 1
          where (.not. lo_bnd) lo = lo - 1

          hi = [nx] * nc
          where (.not. hi_bnd) hi = hi + 1

          ! Include ghost cells around internal boundaries
          allocate(var_data(lo(1):hi(1), n_cc+n_add))

          do ix = 1, nx
             id = box_list(ix)

             cc(:, 1:n_cc) = tree%boxes(id)%cc(:, icc_val)
             if (present(add_vars)) then
                call add_vars(tree%boxes(id), &
                     cc(:, n_cc+1:n_cc+n_add), n_add)
             end if

             ! Include ghost cells on internal block boundaries
             blo = 1
             where ([ix] == 1 .and. .not. lo_bnd) blo = 0

             bhi = nc
             where ([ix] == [nx] .and. .not. hi_bnd) bhi = nc+1

             vlo = blo + ([ix]-1) * nc
             vhi = bhi + ([ix]-1) * nc

             var_data(vlo(1):vhi(1), :) = cc(blo(1):bhi(1), :)
          end do

          id = box_list(1)
          dr = tree%boxes(id)%dr
          r_min = tree%boxes(id)%r_min - (1 - lo) * dr

          write(grid_list(i_grid), "(A,I0)") meshdir // '/' // grid_name, i_grid
          call SILO_add_grid(dbix, grid_list(i_grid), 1, &
               hi - lo + 2, r_min, dr, 1-lo, hi - [nx] * nc)
          write(grid_list_block(i_grid), "(A,I0)") meshdir // '/' // block_prefix &
               // grid_name, i_grid
          call SILO_add_grid(dbix, grid_list_block(i_grid), 1, [nx+1], &
               tree%boxes(id)%r_min, nc*dr, [0], [0])

          do iv = 1, n_cc+n_add
             write(var_list(iv, i_grid), "(A,I0)") meshdir // '/' // &
                  trim(var_names(iv)) // "_", i_grid
             call SILO_add_var(dbix, var_list(iv, i_grid), grid_list(i_grid), &
                  pack(var_data(:, iv), .true.), hi-lo+1)
          end do

          deallocate(var_data)
          deallocate(box_list)
#elif NDIM == 2
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
             nb_ids = tree%boxes(ids)%neighbors(af_neighb_lowx)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(1) < tree%boxes(ids(1))%ix(1) .and. &
                     .not. any(af_has_children(tree%boxes(nb_ids)))) then
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
             nb_ids = tree%boxes(ids)%neighbors(af_neighb_highx)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(1) > tree%boxes(ids(1))%ix(1) .and. &
                     .not. any(af_has_children(tree%boxes(nb_ids)))) then
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
             nb_ids = tree%boxes(ids)%neighbors(af_neighb_lowy)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(2) < tree%boxes(ids(1))%ix(2) .and. &
                     .not. any(af_has_children(tree%boxes(nb_ids)))) then
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
             nb_ids = tree%boxes(ids)%neighbors(af_neighb_highy)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(2) > tree%boxes(ids(1))%ix(2) .and. &
                     .not. any(af_has_children(tree%boxes(nb_ids)))) then
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
          lo_bnd = af_is_phys_boundary(tree%boxes, id, af_low_neighbs)

          id = box_list(nx, ny)
          hi_bnd = af_is_phys_boundary(tree%boxes, id, af_high_neighbs)

          lo(:) = 1
          where (.not. lo_bnd) lo = lo - 1

          hi = [nx, ny] * nc
          where (.not. hi_bnd) hi = hi + 1

          ! Include ghost cells around internal boundaries
          allocate(var_data(lo(1):hi(1), lo(2):hi(2), n_cc+n_add))

          do ix = 1, nx
             do iy = 1, ny
                id = box_list(ix, iy)

                cc(:, :, 1:n_cc) = tree%boxes(id)%cc(:, :, icc_val)
                if (present(add_vars)) then
                   call add_vars(tree%boxes(id), &
                        cc(:, :, n_cc+1:n_cc+n_add), n_add)
                end if

                ! Include ghost cells on internal block boundaries
                blo = 1
                where ([ix, iy] == 1 .and. .not. lo_bnd) blo = 0

                bhi = nc
                where ([ix, iy] == [nx, ny] .and. .not. hi_bnd) bhi = nc+1

                vlo = blo + ([ix, iy]-1) * nc
                vhi = bhi + ([ix, iy]-1) * nc

                var_data(vlo(1):vhi(1), vlo(2):vhi(2), :) = &
                     cc(blo(1):bhi(1), blo(2):bhi(2), :)
             end do
          end do

          id = box_list(1, 1)
          dr = tree%boxes(id)%dr
          r_min = tree%boxes(id)%r_min - (1 - lo) * dr

          write(grid_list(i_grid), "(A,I0)") meshdir // '/' // grid_name, i_grid
          call SILO_add_grid(dbix, grid_list(i_grid), 2, &
               hi - lo + 2, r_min, dr, 1-lo, hi - [nx, ny] * nc, n_cycle_val)
          write(grid_list_block(i_grid), "(A,I0)") meshdir // '/' // block_prefix &
               // grid_name, i_grid
          call SILO_add_grid(dbix, grid_list_block(i_grid), 2, [nx+1, ny+1], &
               tree%boxes(id)%r_min, nc*dr, [0, 0], [0, 0], n_cycle_val)

          do iv = 1, n_cc+n_add
             write(var_list(iv, i_grid), "(A,I0)") meshdir // '/' // &
                  trim(var_names(iv)) // "_", i_grid
             call SILO_add_var(dbix, var_list(iv, i_grid), grid_list(i_grid), &
                  pack(var_data(:, :, iv), .true.), hi-lo+1, n_cycle_val)
          end do

          deallocate(var_data)
          deallocate(box_list)
#elif NDIM == 3
          allocate(box_list(1,1,1))
          box_list(1,1,1) = id
          box_done(id) = .true.
          nx = 1
          ny = 1
          nz = 1

          do
             nx_prev = nx
             ny_prev = ny
             nz_prev = nz

             ! Check whether we can extend to the -x direction
             ids = pack(box_list(1, :, :), .true.)
             nb_ids = tree%boxes(ids)%neighbors(af_neighb_lowx)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(1) < tree%boxes(ids(1))%ix(1) .and. &
                     .not. any(af_has_children(tree%boxes(nb_ids)))) then
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
             nb_ids = tree%boxes(ids)%neighbors(af_neighb_highx)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(1) > tree%boxes(ids(1))%ix(1) .and. &
                     .not. any(af_has_children(tree%boxes(nb_ids)))) then
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
             nb_ids = tree%boxes(ids)%neighbors(af_neighb_lowy)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(2) < tree%boxes(ids(1))%ix(2) .and. &
                     .not. any(af_has_children(tree%boxes(nb_ids)))) then
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
             nb_ids = tree%boxes(ids)%neighbors(af_neighb_highy)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(2) > tree%boxes(ids(1))%ix(2) .and. &
                     .not. any(af_has_children(tree%boxes(nb_ids)))) then
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
             nb_ids = tree%boxes(ids)%neighbors(af_neighb_lowz)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(3) < tree%boxes(ids(1))%ix(3) .and. &
                     .not. any(af_has_children(tree%boxes(nb_ids)))) then
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
             nb_ids = tree%boxes(ids)%neighbors(af_neighb_highz)
             if (all(nb_ids > af_no_box)) then
                if (.not. any(box_done(nb_ids)) .and. &
                     tree%boxes(nb_ids(1))%ix(3) > tree%boxes(ids(1))%ix(3) .and. &
                     .not. any(af_has_children(tree%boxes(nb_ids)))) then
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
          lo_bnd = af_is_phys_boundary(tree%boxes, id, af_low_neighbs)

          id = box_list(nx, ny, nz)
          hi_bnd = af_is_phys_boundary(tree%boxes, id, af_high_neighbs)

          lo(:) = 1
          where (.not. lo_bnd) lo = lo - 1

          hi = [nx, ny, nz] * nc
          where (.not. hi_bnd) hi = hi + 1

          ! Include ghost cells around internal boundaries
          allocate(var_data(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3), n_cc+n_add))

          do iz = 1, nz
             do ix = 1, nx
                do iy = 1, ny
                   id     = box_list(ix, iy, iz)
                   cc(:, :, :, 1:n_cc) = tree%boxes(id)%cc(:, :, :, icc_val)
                   if (present(add_vars)) then
                      call add_vars(tree%boxes(id), &
                           cc(:, :, :, n_cc+1:n_cc+n_add), n_add)
                   end if

                   ! Include ghost cells on internal block boundaries
                   blo = 1
                   where ([ix, iy, iz] == 1 .and. .not. lo_bnd) blo = 0

                   bhi = nc
                   where ([ix, iy, iz] == [nx, ny, nz] &
                        .and. .not. hi_bnd) bhi = nc+1

                   vlo = blo + ([ix, iy, iz]-1) * nc
                   vhi = bhi + ([ix, iy, iz]-1) * nc

                   var_data(vlo(1):vhi(1), vlo(2):vhi(2), vlo(3):vhi(3), :) = &
                        cc(blo(1):bhi(1), blo(2):bhi(2), blo(3):bhi(3), :)
                end do
             end do
          end do

          id = box_list(1, 1, 1)
          dr = tree%boxes(id)%dr
          r_min = tree%boxes(id)%r_min - (1 - lo) * dr

          write(grid_list(i_grid), "(A,I0)") meshdir // '/' // grid_name, i_grid
          call SILO_add_grid(dbix, grid_list(i_grid), 3, &
               hi - lo + 2, r_min, dr, 1-lo, hi-[nx, ny, nz]*nc, n_cycle_val)
          write(grid_list_block(i_grid), "(A,I0)") meshdir // '/' // block_prefix // &
               grid_name, i_grid
          call SILO_add_grid(dbix, grid_list_block(i_grid), 3, [nx+1, ny+1, nz+1], &
               tree%boxes(id)%r_min, nc*dr, [0, 0, 0], [0, 0, 0], n_cycle_val)

          do iv = 1, n_cc+n_add
             write(var_list(iv, i_grid), "(A,I0)") meshdir // '/' // &
                  trim(var_names(iv)) // "_", i_grid
             call SILO_add_var(dbix, var_list(iv, i_grid), grid_list(i_grid), &
                  pack(var_data(:, :, :, iv), .true.), hi-lo+1, n_cycle_val)
          end do

          deallocate(var_data)
          deallocate(box_list)
#endif

       end do
    end do

    call SILO_set_mmesh_grid(dbix, amr_name, grid_list(1:i_grid), &
         n_cycle_val, time_val)
    do iv = 1, n_cc+n_add
       call SILO_set_mmesh_var(dbix, trim(var_names(iv)), amr_name, &
            var_list(iv, 1:i_grid), n_cycle_val, time_val)
    end do

    call SILO_set_mmesh_grid(dbix, block_prefix // amr_name, &
         grid_list_block(1:i_grid), n_cycle_val, time_val)
    call SILO_close_file(dbix)
    print *, "af_write_silo: written " // trim(fname)
  end subroutine af_write_silo
#endif

  subroutine get_output_vars(tree, ix_out)
    type(af_t), intent(in)              :: tree
    integer, allocatable, intent(inout) :: ix_out(:)
    integer                             :: n, i

    n = count(tree%cc_write_output(1:tree%n_var_cell))
    allocate(ix_out(n))

    n = 0
    do i = 1, tree%n_var_cell
       if (tree%cc_write_output(i)) then
          n = n + 1
          ix_out(n) = i
       end if
    end do
  end subroutine get_output_vars

end module m_af_output
