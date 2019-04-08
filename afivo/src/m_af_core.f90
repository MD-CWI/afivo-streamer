#include "cpp_macros.h"
!> This module contains the core routines of Afivo, namely those that deal with
!> initializing and changing the quadtree/octree mesh.
module m_af_core
  use m_af_types

  implicit none
  private

  public :: af_add_cc_variable
  public :: af_add_fc_variable
  public :: af_find_cc_variable
  public :: af_find_fc_variable
  public :: af_init
  public :: af_set_cc_methods
  public :: af_init_box
  public :: af_destroy
  public :: af_adjust_refinement
  public :: af_consistent_fluxes

contains

  !> Add cell-centered variable
  subroutine af_add_cc_variable(tree, name, write_out, n_copies, &
       max_lvl, ix, write_binary)
    !> Tree to add variable to
    type(af_t), intent(inout)      :: tree
    !> Name of the variable
    character(len=*), intent(in)   :: name
    !> Include variable in output
    logical, intent(in), optional  :: write_out
    !> Include variable in binary output (for restarting)
    logical, intent(in), optional  :: write_binary
    !> How many copies of variable to store (default: 1)
    integer, intent(in), optional  :: n_copies
    !> Store variable up to this refinement level (default: af_max_lvl)
    integer, intent(in), optional  :: max_lvl
    !> On output: index of variable
    integer, intent(out), optional :: ix

    integer :: n
    integer :: ncpy, maxlvl
    logical :: writeout, writebin

    ncpy = 1; if (present(n_copies)) ncpy = n_copies
    writeout = .true.; if (present(write_out)) writeout = write_out
    writebin = .true.; if (present(write_binary)) writebin = write_binary
    maxlvl = 100; if (present(max_lvl)) maxlvl = max_lvl

    if (ncpy < 1) error stop "af_add_cc_variable: n_copies < 1"
    if (maxlvl < 1) error stop "af_add_cc_variable: max_lvl < 1"

    do n = 1, ncpy
       tree%n_var_cell = tree%n_var_cell + 1
       if (n == 1) then
          if (present(ix)) ix = tree%n_var_cell
          tree%cc_names(tree%n_var_cell)        = name
          tree%cc_write_output(tree%n_var_cell) = writeout
          tree%cc_write_binary(tree%n_var_cell) = writebin
          tree%cc_max_level(tree%n_var_cell)    = maxlvl
          tree%cc_num_copies(tree%n_var_cell)   = ncpy
       else
          write(tree%cc_names(tree%n_var_cell), "(A,I0)") &
               trim(name) // '_', n
          tree%cc_write_output(tree%n_var_cell) = .false.
          tree%cc_write_binary(tree%n_var_cell) = .false.
          tree%cc_max_level(tree%n_var_cell)    = maxlvl
          tree%cc_num_copies(tree%n_var_cell)   = 1
       end if
    end do

  end subroutine af_add_cc_variable

  !> Add face-centered variable
  subroutine af_add_fc_variable(tree, name, ix, write_binary)
    !> Tree to add variable to
    type(af_t), intent(inout)      :: tree
    !> Name of the variable
    character(len=*), intent(in)   :: name
    !> On output: index of variable
    integer, intent(out), optional :: ix
    !> Include variable in binary output
    logical, intent(in), optional  :: write_binary
    logical                        :: writebin

    writebin = .true.; if (present(write_binary)) writebin = write_binary

    tree%n_var_face                       = tree%n_var_face + 1
    tree%fc_names(tree%n_var_face)        = name
    tree%fc_write_binary(tree%n_var_face) = writebin
    if (present(ix)) ix = tree%n_var_face
  end subroutine af_add_fc_variable

  !> Find index of cell-centered variable
  integer function af_find_cc_variable(tree, name)
    type(af_t), intent(in)       :: tree
    character(len=*), intent(in) :: name
    integer                      :: n

    do n = 1, tree%n_var_cell
       if (tree%cc_names(n) == name) exit
    end do

    if (n == tree%n_var_cell+1) then
       print *, "variable name: ", trim(name)
       error stop "af_find_cc_variable: variable not found"
    end if

    af_find_cc_variable = n
  end function af_find_cc_variable

  !> Find index of face-centered variable
  integer function af_find_fc_variable(tree, name)
    type(af_t), intent(in)       :: tree
    character(len=*), intent(in) :: name
    integer                      :: n

    do n = 1, tree%n_var_face
       if (tree%fc_names(n) == name) exit
    end do

    if (n == tree%n_var_face+1) then
       print *, "variable name: ", trim(name)
       error stop "af_find_fc_variable: variable not found"
    end if

    af_find_fc_variable = n
  end function af_find_fc_variable

  !> Initialize a NDIM-d octree/quadtree grid
  subroutine af_init(tree, n_cell, r_max, grid_size, periodic, r_min, coord, mem_limit_gb)
    type(af_t), intent(inout)      :: tree            !< The tree to initialize
    integer, intent(in)            :: n_cell          !< Boxes have n_cell^dim cells
    real(dp), intent(in)           :: r_max(NDIM)     !< Maximal coordinates of the domain
    integer, intent(in)            :: grid_size(NDIM) !< Size of the coarse grid
    logical, intent(in), optional  :: periodic(NDIM)  !< True for periodic dimensions
    real(dp), intent(in), optional :: r_min(NDIM)     !< Lowest coordinate, default is (0., 0., 0.)
    integer, intent(in), optional  :: coord           !< Select coordinate type
    real(dp), intent(in), optional :: mem_limit_gb    !< Memory limit in GByte

    real(dp)                       :: r_min_a(NDIM), gb_limit
    integer                        :: lvl, coord_a, box_bytes

    ! Set default arguments if not present
    r_min_a = 0.0_dp;  if (present(r_min)) r_min_a = r_min
    coord_a = af_xyz;  if (present(coord)) coord_a = coord
    gb_limit = 4;      if (present(mem_limit_gb)) gb_limit = mem_limit_gb

    if (tree%ready)       stop "af_init: tree was already initialized"
    if (n_cell < 2)       stop "af_init: n_cell should be >= 2"
    if (btest(n_cell, 0)) stop "af_init: n_cell should be even"
    if (gb_limit <= 0)    stop "af_init: mem_limit_gb should be > 0"
#if NDIM == 3
    if (coord_a == af_cyl) stop "af_init: cannot have 3d cyl coords"
#endif
    if (tree%n_var_cell <= 0) stop "af_init: no cell-centered variables present"

    do lvl = af_min_lvl, af_max_lvl
       allocate(tree%lvls(lvl)%ids(0))
       allocate(tree%lvls(lvl)%leaves(0))
       allocate(tree%lvls(lvl)%parents(0))
    end do

    tree%n_cell      = n_cell
    tree%r_base      = r_min_a
    tree%dr_base     = (r_max - r_min_a) / grid_size
    tree%highest_id  = 0
    tree%highest_lvl = 0
    tree%coord_t     = coord_a

    ! Calculate size of a box
    box_bytes = af_box_bytes(n_cell, tree%n_var_cell, tree%n_var_face)
    tree%box_limit = nint(gb_limit * 2.0_dp**30 / box_bytes)

    ! Allocate the full list of boxes
    allocate(tree%boxes(tree%box_limit))

    ! This list can probably be a bit smaller
    tree%n_removed_ids = 0
    allocate(tree%removed_ids(tree%box_limit))

    ! Initialize list of cell-centered variables with methods
    if (.not. allocated(tree%cc_method_vars)) &
         allocate(tree%cc_method_vars(0))

    call af_set_coarse_grid(tree, grid_size, periodic)

  end subroutine af_init

  !> Create the coarse grid
  subroutine af_set_coarse_grid(tree, coarse_grid_size, periodic_dims)
    !> Tree for which we set the base
    type(af_t), intent(inout)     :: tree
    !> Size of coarse grid (in cells)
    integer, intent(in)           :: coarse_grid_size(NDIM)
    !> Whether dimensions are periodic (default: false)
    logical, intent(in), optional :: periodic_dims(NDIM)
    logical                       :: periodic(NDIM)
    integer                       :: nx(NDIM), ix(NDIM), IJK, id, n_boxes, nb
    integer, allocatable          :: id_array(DTIMES(:))

    if (tree%highest_id > 0) &
         error stop "af_set_coarse_grid: this tree already has boxes"
    if (.not. allocated(tree%boxes)) &
         error stop "af_set_coarse_grid: tree not initialized"
    if (any(coarse_grid_size < tree%n_cell)) &
         error stop "af_set_coarse_grid: coarse_grid_size < tree%n_cell"
    if (any(modulo(coarse_grid_size, tree%n_cell) /= 0)) &
         error stop "af_set_coarse_grid: coarse_grid_size not divisible by tree%n_cell"

    periodic(:) = .false.; if (present(periodic_dims)) periodic = periodic_dims

    tree%coarse_grid_size(1:NDIM) = coarse_grid_size
    tree%periodic(1:NDIM)         = periodic
    tree%ready                    = .true.

    nx      = coarse_grid_size / tree%n_cell
    n_boxes = product(nx)

    ! For easy lookup of box neighbors
    call create_index_array(nx, periodic, id_array)

    ! Check if we have enough space
    if (n_boxes > size(tree%boxes(:))) &
         error stop "Not enough memory available for coarse grid"

    ! Create level 1
    deallocate(tree%lvls(1)%ids)
    allocate(tree%lvls(1)%ids(n_boxes))

    ! The ids are simply 1, 2, 3, ..., N
    call get_free_ids(tree, tree%lvls(1)%ids)
    tree%lvls(1)%leaves = tree%lvls(1)%ids

    ! Loop over the boxes and set their neighbors
#if NDIM == 2
    do j = 1, nx(2)
       do i = 1, nx(1)
          id                     = id_array(IJK)
          tree%boxes(id)%lvl     = 1
          tree%boxes(id)%ix      = [IJK]
          tree%boxes(id)%dr      = tree%dr_base
          tree%boxes(id)%r_min   = tree%r_base + &
               (tree%boxes(id)%ix - 1) * tree%dr_base * tree%n_cell
          tree%boxes(id)%n_cell  = tree%n_cell
          tree%boxes(id)%coord_t = tree%coord_t

          tree%boxes(id)%parent      = af_no_box
          tree%boxes(id)%children(:) = af_no_box

          ! Connectivity
          do nb = 1, af_num_neighbors
             ix = [IJK] + af_neighb_dix(:, nb)
             tree%boxes(id)%neighbors(nb) = id_array(ix(1), ix(2))
          end do
          tree%boxes(id)%neighbor_mat = id_array(i-1:i+1, j-1:j+1)

          call af_init_box(tree%boxes(id), tree%boxes(id)%n_cell, &
               tree%n_var_cell, tree%n_var_face)
       end do
    end do
#elif NDIM == 3
    do k = 1, nx(3)
       do j = 1, nx(2)
          do i = 1, nx(1)
             id                     = id_array(IJK)
             tree%boxes(id)%lvl     = 1
             tree%boxes(id)%ix      = [IJK]
             tree%boxes(id)%dr      = tree%dr_base
             tree%boxes(id)%r_min   = tree%r_base + &
                  (tree%boxes(id)%ix - 1) * tree%dr_base * tree%n_cell
             tree%boxes(id)%n_cell  = tree%n_cell
             tree%boxes(id)%coord_t = tree%coord_t

             tree%boxes(id)%parent      = af_no_box
             tree%boxes(id)%children(:) = af_no_box

             ! Connectivity
             do nb = 1, af_num_neighbors
                ix = [IJK] + af_neighb_dix(:, nb)
                tree%boxes(id)%neighbors(nb) = id_array(ix(1), ix(2), ix(3))
             end do
             tree%boxes(id)%neighbor_mat = id_array(i-1:i+1, j-1:j+1, k-1:k+1)

             call af_init_box(tree%boxes(id), tree%boxes(id)%n_cell, &
                  tree%n_var_cell, tree%n_var_face)
          end do
       end do
    end do
#endif

    tree%highest_lvl = 1

  end subroutine af_set_coarse_grid

  !> Set the methods for a cell-centered variable
  subroutine af_set_cc_methods(tree, iv, bc, rb, prolong, restrict)
    use m_af_ghostcell, only: af_gc_interp
    use m_af_prolong, only: af_prolong_linear
    use m_af_restrict, only: af_restrict_box
    type(af_t), intent(inout)             :: tree     !< Tree to operate on
    integer, intent(in)                   :: iv       !< Index of variable
    procedure(af_subr_bc)                 :: bc       !< Boundary condition method
    procedure(af_subr_rb), optional       :: rb       !< Refinement boundary method
    procedure(af_subr_prolong), optional  :: prolong  !< Prolongation method
    procedure(af_subr_restrict), optional :: restrict !< Restriction method
    integer                               :: i, n

    tree%cc_methods(iv)%bc => bc

    if (present(rb)) then
       tree%cc_methods(iv)%rb => rb
    else
       tree%cc_methods(iv)%rb => af_gc_interp
    end if

    if (present(prolong)) then
       tree%cc_methods(iv)%prolong => prolong
    else
       tree%cc_methods(iv)%prolong => af_prolong_linear
    end if

    if (present(restrict)) then
       tree%cc_methods(iv)%restrict => restrict
    else
       tree%cc_methods(iv)%restrict => af_restrict_box
    end if

    tree%has_cc_method(iv) = .true.

    if (.not. allocated(tree%cc_method_vars)) &
         allocate(tree%cc_method_vars(0))

    n = size(tree%cc_method_vars)
    tree%cc_method_vars = [(tree%cc_method_vars(i), i=1,n), iv]

  end subroutine af_set_cc_methods

  !> "Destroy" the data in a tree. Since we don't use pointers, you can also
  !> just let a tree get out of scope
  subroutine af_destroy(tree)
    type(af_t), intent(inout) :: tree
    integer                   :: lvl

    if (.not. tree%ready) stop "af_destroy: Tree not fully initialized"
    deallocate(tree%boxes)
    deallocate(tree%removed_ids)
    deallocate(tree%cc_method_vars)

    do lvl = af_min_lvl, af_max_lvl
       deallocate(tree%lvls(lvl)%ids)
       deallocate(tree%lvls(lvl)%leaves)
       deallocate(tree%lvls(lvl)%parents)
    end do

    tree%highest_lvl = 0
    tree%highest_id  = 0
    tree%n_var_cell  = 0
    tree%n_var_face  = 0
    tree%n_removed_ids = 0
    tree%ready       = .false.
  end subroutine af_destroy

  !> Create an array for easy lookup of indices
  subroutine create_index_array(nx, periodic, id_array)
    integer, intent(in)                 :: nx(NDIM)
    logical, intent(in)                 :: periodic(NDIM)
    integer, intent(inout), allocatable :: id_array(DTIMES(:))
    integer                             :: IJK

#if NDIM == 2
    allocate(id_array(0:nx(1)+1, 0:nx(2)+1))
#elif NDIM == 3
    allocate(id_array(0:nx(1)+1, 0:nx(2)+1, 0:nx(3)+1))
#endif

    id_array = af_phys_boundary

#if NDIM == 2
    do j = 1, nx(2)
       do i = 1, nx(1)
          id_array(i, j) = (j-1) * nx(1) + i
       end do
    end do

    if (periodic(1)) then
       id_array(0, :) = id_array(nx(1), :)
       id_array(nx(1)+1, :) = id_array(1, :)
    end if

    if (periodic(2)) then
       id_array(:, 0) = id_array(:, nx(2))
       id_array(:, nx(2)+1) = id_array(:, 1)
    end if
#elif NDIM == 3
    do k = 1, nx(3)
       do j = 1, nx(2)
          do i = 1, nx(1)
             id_array(i, j, k) = (k-1) * nx(2) * nx(1) + (j-1) * nx(1) + i
          end do
       end do
    end do

    if (periodic(1)) then
       id_array(0, :, :) = id_array(nx(1), :, :)
       id_array(nx(1)+1, :, :) = id_array(1, :, :)
    end if

    if (periodic(2)) then
       id_array(:, 0, :) = id_array(:, nx(2), :)
       id_array(:, nx(2)+1, :) = id_array(:, 1, :)
    end if

    if (periodic(3)) then
       id_array(:, :, 0) = id_array(:, :, nx(3))
       id_array(:, :, nx(3)+1) = id_array(:, :, 1)
    end if
#endif
  end subroutine create_index_array

  !> Create a list of leaves and a list of parents for a level
  subroutine set_leaves_parents(boxes, level)
    type(box_t), intent(in)   :: boxes(:) !< List of boxes
    type(lvl_t), intent(inout) :: level !< Level type which contains the indices of boxes
    integer                    :: i, id, i_leaf, i_parent
    integer                    :: n_parents, n_leaves

    n_parents = count(af_has_children(boxes(level%ids)))
    n_leaves = size(level%ids) - n_parents

    if (n_parents /= size(level%parents)) then
       deallocate(level%parents)
       allocate(level%parents(n_parents))
    end if

    if (n_leaves /= size(level%leaves)) then
       deallocate(level%leaves)
       allocate(level%leaves(n_leaves))
    end if

    i_leaf   = 0
    i_parent = 0
    do i = 1, size(level%ids)
       id = level%ids(i)
       if (af_has_children(boxes(id))) then
          i_parent                = i_parent + 1
          level%parents(i_parent) = id
       else
          i_leaf               = i_leaf + 1
          level%leaves(i_leaf) = id
       end if
    end do
  end subroutine set_leaves_parents

  !> Mark box as active and allocate data storage for a box, for its cell- and
  !> face-centered data
  subroutine af_init_box(box, n_cell, n_cc, n_fc)
    type(box_t), intent(inout) :: box !< Box for which we allocate memory
    integer, intent(in)         :: n_cell !< Number of cells per dimension in the box
    integer, intent(in)         :: n_cc   !< Number of cell-centered variables
    integer, intent(in)         :: n_fc   !< Number of face-centered variables

    box%in_use = .true.

    ! Sometimes we re-use a removed box, then we don't have to re-allocate
    if (.not. allocated(box%cc)) then
#if NDIM == 2
       allocate(box%cc(0:n_cell+1, 0:n_cell+1, n_cc))
       allocate(box%fc(n_cell+1,   n_cell+1, NDIM, n_fc))
#elif NDIM == 3
       allocate(box%cc(0:n_cell+1, 0:n_cell+1, 0:n_cell+1, n_cc))
       allocate(box%fc(n_cell+1,   n_cell+1,   n_cell+1, NDIM, n_fc))
#endif
    end if

    ! Initialize to zero
    box%cc = 0
    box%fc = 0
  end subroutine af_init_box

  ! Set the neighbors of id (using their parent)
  subroutine set_neighbs(boxes, id)
    type(box_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    integer                     :: nb, nb_id, IJK

    do KJI_DO(-1, 1)
       if (boxes(id)%neighbor_mat(IJK) == af_no_box) then
          nb_id = find_neighb(boxes, id, [IJK])
          if (nb_id > af_no_box) then
             boxes(id)%neighbor_mat(IJK) = nb_id
#if NDIM == 2
             boxes(nb_id)%neighbor_mat(-i, -j) = id
#elif NDIM == 3
             boxes(nb_id)%neighbor_mat(-i, -j, -k) = id
#endif
          end if
       end if
    end do; CLOSE_DO

    do nb = 1, af_num_neighbors
       if (boxes(id)%neighbors(nb) == af_no_box) then
#if NDIM == 2
          nb_id = boxes(id)%neighbor_mat(af_neighb_dix(1, nb), &
               af_neighb_dix(2, nb))
#elif NDIM == 3
          nb_id = boxes(id)%neighbor_mat(af_neighb_dix(1, nb), &
               af_neighb_dix(2, nb), af_neighb_dix(3, nb))
#endif
          if (nb_id > af_no_box) then
             boxes(id)%neighbors(nb) = nb_id
             boxes(nb_id)%neighbors(af_neighb_rev(nb)) = id
          end if
       end if
    end do
  end subroutine set_neighbs

  !> Get the id of all neighbors of boxes(id), through its parent
  function find_neighb(boxes, id, dix) result(nb_id)
    type(box_t), intent(in) :: boxes(:) !< List with all the boxes
    integer, intent(in)       :: id       !< Box whose neighbor we are looking for
    integer, intent(in)       :: dix(NDIM)
    integer                   :: nb_id, p_id, c_ix, dix_c(NDIM)

    p_id = boxes(id)%parent
    c_ix = af_ix_to_ichild(boxes(id)%ix)

    ! Check if neighbor is in same direction as dix is (low/high). If so, use
    ! neighbor of parent
    where ((dix == -1) .eqv. af_child_low(:, c_ix))
       dix_c = dix
    elsewhere
       dix_c = 0
    end where

#if NDIM == 2
    p_id = boxes(p_id)%neighbor_mat(dix_c(1), dix_c(2))
#elif NDIM == 3
    p_id = boxes(p_id)%neighbor_mat(dix_c(1), dix_c(2), dix_c(3))
#endif

    if (p_id <= af_no_box) then
       nb_id = p_id
    else
       c_ix = af_ix_to_ichild(boxes(id)%ix + dix)
       nb_id = boxes(p_id)%children(c_ix)
    end if
  end function find_neighb

  !> Adjust the refinement of a tree using the user-supplied ref_subr. The
  !> optional argument ref_buffer controls over how many cells neighbors are
  !> affected by refinement flags.
  !>
  !> On input, the tree should be balanced. On output, the tree is still
  !> balanced, and its refinement is updated (with at most one level per call).
  subroutine af_adjust_refinement(tree, ref_subr, ref_info, ref_buffer, &
       modify_refinement)
    type(af_t), intent(inout)      :: tree        !< The tree to adjust
    procedure(af_subr_ref)         :: ref_subr    !< Refinement function
    type(ref_info_t), intent(inout) :: ref_info    !< Information about refinement
    integer, intent(in), optional   :: ref_buffer  !< Buffer width (in cells)
    !> A routine to manually override refinement flags
    procedure(subr_modify_ref), optional :: modify_refinement
    integer                         :: lvl, id, i, c_ids(af_num_children), i_ch
    integer, allocatable            :: ref_flags(:)
    integer                         :: ref_buffer_val

    if (.not. tree%ready) stop "Tree not ready"

    ref_buffer_val = 0          ! Default buffer width (in cells) around refinement
    if (present(ref_buffer)) ref_buffer_val = ref_buffer

    if (ref_buffer_val < 0) &
         error stop "af_adjust_refinement: ref_buffer < 0"
    if (ref_buffer_val > tree%n_cell) &
         error stop "af_adjust_refinement: ref_buffer > tree%n_cell"

    allocate(ref_flags(tree%highest_id))

    ! Set refinement values for all boxes. Only two flags are used below:
    ! af_refine and af_derefine. Other values are ignored.
    call consistent_ref_flags(tree, ref_flags, ref_subr, &
         ref_buffer_val, modify_refinement)

    do lvl = 1, af_max_lvl-1
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)

          if (id > size(ref_flags)) then
             cycle              ! This is a newly added box
          else if (ref_flags(id) == af_refine) then
             ! Add children. First need to get num_children free id's
             call get_free_ids(tree, c_ids)
             call add_children(tree%boxes, id, c_ids, &
                  tree%n_var_cell, tree%n_var_face)
          else if (ref_flags(id) == af_derefine) then
             ! Remove children
             call auto_restrict(tree, id)
             call remove_children(tree, id)
          end if
       end do

       ! Update leaves / parents
       call set_leaves_parents(tree%boxes, tree%lvls(lvl))

       ! Set next level ids to children of this level
       call set_child_ids(tree%lvls(lvl)%parents, &
            tree%lvls(lvl+1)%ids, tree%boxes)

       ! Update connectivity of new children
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          if (ref_flags(id) == af_refine) then
             do i_ch = 1, af_num_children
                call set_neighbs(tree%boxes, tree%boxes(id)%children(i_ch))
             end do
          end if
       end do

       if (size(tree%lvls(lvl+1)%ids) == 0) exit
    end do

    tree%highest_lvl = lvl

    ! We still have to update leaves and parents for the last level, which is
    ! either lvl+1 or af_max_lvl. Note that lvl+1 is empty now, but maybe it was
    ! not not empty before, and that af_max_lvl is skipped in the above loop.
    lvl = min(lvl+1, af_max_lvl)
    call set_leaves_parents(tree%boxes, tree%lvls(lvl))

    ! Set information about the refinement
    call set_ref_info(tree, ref_flags, ref_info)

    call auto_prolong(tree, ref_info)

  end subroutine af_adjust_refinement

  !> Try to automatically restrict to box with index id
  subroutine auto_restrict(tree, id)
    type(af_t), intent(inout) :: tree
    integer, intent(in)        :: id
    integer                    :: iv, i_ch, ch_id

    if (.not. any(tree%has_cc_method(:))) return

    do iv = 1, tree%n_var_cell
       if (tree%has_cc_method(iv)) then
          do i_ch = 1, af_num_children
             ch_id = tree%boxes(id)%children(i_ch)
             call tree%cc_methods(iv)%restrict(tree%boxes(ch_id), &
                  tree%boxes(id), iv)
          end do
       end if
    end do
  end subroutine auto_restrict

  !> Try to automatically prolong to all new boxes
  subroutine auto_prolong(tree, ref_info)
    use m_af_ghostcell, only: af_gc_box
    type(af_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, n, iv, id, p_id

    ! Skip this routine when it won't do anything
    if (.not. any(tree%has_cc_method(:)) .or. ref_info%n_add == 0) then
       return
    end if

    !$omp parallel private(lvl, i, n, iv, id, p_id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          p_id = tree%boxes(id)%parent

          do n = 1, size(tree%cc_method_vars)
             iv = tree%cc_method_vars(n)
             call tree%cc_methods(iv)%prolong(tree%boxes(p_id), &
                  tree%boxes(id), iv)
          end do
       end do
       !$omp end do

       !$omp do
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call af_gc_box(tree, id, [tree%cc_method_vars])
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine auto_prolong

  !> Set information about the refinement for all "normal" levels (>= 1)
  subroutine set_ref_info(tree, ref_flags, ref_info)
    type(af_t), intent(in)         :: tree
    integer, intent(in)             :: ref_flags(:)
    type(ref_info_t), intent(inout) :: ref_info
    integer                         :: id, lvl, n, n_ch
    integer, allocatable            :: ref_count(:), drf_count(:)

    n_ch           = af_num_children
    ref_info%n_add = n_ch * count(ref_flags == af_refine)
    ref_info%n_rm  = n_ch * count(ref_flags == af_derefine)

    ! Use highest_lvl+1 here because this lvl might have been completely removed
    if (allocated(ref_info%lvls)) deallocate(ref_info%lvls)
    allocate(ref_info%lvls(tree%highest_lvl+1))
    allocate(ref_count(tree%highest_lvl+1))
    allocate(drf_count(tree%highest_lvl+1))

    ! Find the number of (de)refined boxes per level
    ref_count = 0
    drf_count = 0

    do id = 1, size(ref_flags)
       lvl = tree%boxes(id)%lvl

       if (ref_flags(id) == af_refine) then
          ref_count(lvl) = ref_count(lvl) + 1
       else if (ref_flags(id) == af_derefine) then
          drf_count(lvl) = drf_count(lvl) + 1
       end if
    end do

    ! Allocate storage per level
    ! There can be no new children at level 1
    allocate(ref_info%lvls(1)%add(0))
    allocate(ref_info%lvls(1)%rm(0))

    do lvl = 2, tree%highest_lvl+1
       n = ref_count(lvl-1) * n_ch
       allocate(ref_info%lvls(lvl)%add(n))
       n = drf_count(lvl-1) * n_ch
       allocate(ref_info%lvls(lvl)%rm(n))
    end do

    ! Set the added and removed id's per level, these are the children of the
    ! (de)refined boxes
    ref_count = 0
    drf_count = 0

    do id = 1, size(ref_flags)
       lvl = tree%boxes(id)%lvl

       if (ref_flags(id) == af_refine) then
          ref_count(lvl) = ref_count(lvl) + 1
          n = n_ch * (ref_count(lvl)-1) + 1
          ref_info%lvls(lvl+1)%add(n:n+n_ch-1) = tree%boxes(id)%children
       else if (ref_flags(id) == af_derefine) then
          drf_count(lvl) = drf_count(lvl) + 1
          n = n_ch * (drf_count(lvl)-1) + 1
          ref_info%lvls(lvl+1)%rm(n:n+n_ch-1) = tree%boxes(id)%children
       end if
    end do
  end subroutine set_ref_info

  !> Get free ids from the boxes(:) array to store new boxes in. These ids are
  !> always consecutive.
  subroutine get_free_ids(tree, ids)
    type(af_t), intent(inout) :: tree
    integer, intent(out)      :: ids(:) !< Array which will be filled with free box ids
    integer                   :: i, highest_id_prev, n_ids

    n_ids = size(ids)

    !> @todo when doing AMR in parallel, perhaps move some of the code outside
    !> the critical construct

    !$omp critical (crit_free_ids)
    if (n_ids <= tree%n_removed_ids) then
       ! Re-use removed boxes
       do i = 1, n_ids
          ids(i) = tree%removed_ids(tree%n_removed_ids-n_ids+i)
       end do
       tree%n_removed_ids = tree%n_removed_ids - n_ids
    else
       ! Add new boxes at the end of the list
       highest_id_prev = tree%highest_id
       tree%highest_id = tree%highest_id + n_ids

       if (tree%highest_id > size(tree%boxes)) then
          print *, "get_free_ids: exceeding memory limit"
          write(*, '(A,E12.2)') " memory_limit (GByte):     ", &
               tree%box_limit * 0.5_dp**30 * &
               af_box_bytes(tree%n_cell, tree%n_var_cell, tree%n_var_face)
          print *, "memory_limit (boxes):     ", tree%box_limit
          print *, "You can increase the memory limit in your call to af_init"
          print *, "by setting mem_limit_gb to a higher value (in GBytes)"
          error stop
       end if

       ids = [(highest_id_prev + i, i=1,n_ids)]
    end if
    !$omp end critical (crit_free_ids)

  end subroutine get_free_ids

  !> Given the refinement function, return consistent refinement flags, that
  !> ensure that the tree is still balanced. Furthermore, it cannot derefine the
  !> base level, and it cannot refine above af_max_lvl. The argument
  !> ref_flags is changed: for boxes that will be refined it holds af_refine,
  !> for boxes that will be derefined it holds af_derefine
  subroutine consistent_ref_flags(tree, ref_flags, ref_subr, &
       ref_buffer, modify_refinement)
    use omp_lib, only: omp_get_max_threads, omp_get_thread_num
    type(af_t), intent(inout) :: tree         !< Tree for which we set refinement flags
    integer, intent(inout)     :: ref_flags(:) !< List of refinement flags for all boxes(:)
    procedure(af_subr_ref)    :: ref_subr     !< User-supplied refinement function.
    integer, intent(in)        :: ref_buffer   !< Buffer width (in cells)
    !> A routine to manually override refinement flags
    procedure(subr_modify_ref), optional :: modify_refinement
    integer              :: lvl, i, i_ch, ch_id, id
    integer              :: p_id
    integer              :: thread_id
    integer, allocatable :: tmp_flags(:, :)
#if NDIM == 2
    integer              :: cell_flags(tree%n_cell, tree%n_cell)
#elif NDIM == 3
    integer              :: cell_flags(tree%n_cell, tree%n_cell, tree%n_cell)
#endif
    integer, parameter   :: unset_flag = -huge(1)

    ! Set refinement flags for each thread individually, because we sometimes
    ! modify the refinement flags of neighbors
    allocate(tmp_flags(size(ref_flags), omp_get_max_threads()))

    tmp_flags(:, :) = unset_flag

    ! Set refinement flags on all leaves and their immediate parents (on other
    ! boxes the flags would not matter)

    !$omp parallel private(lvl, i, id, p_id, cell_flags, thread_id, i_ch, ch_id)
    thread_id = omp_get_thread_num() + 1

    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)

          call ref_subr(tree%boxes(id), cell_flags)
          call cell_to_ref_flags(cell_flags, tree%n_cell, &
               tmp_flags(:, thread_id), tree, id, ref_buffer)

          ! If the parent exists, and this is the first child which is itself
          ! not refined, set refinement flags for the parent
          if (tree%boxes(id)%lvl > 1) then
             p_id = tree%boxes(id)%parent
             do i_ch = 1, af_ix_to_ichild(tree%boxes(id)%ix)-1
                ch_id = tree%boxes(p_id)%children(i_ch)
                if (.not. af_has_children(tree%boxes(ch_id))) exit
             end do

             if (i_ch == af_ix_to_ichild(tree%boxes(id)%ix)) then
                ! This is the first child which is itself not refined
                call ref_subr(tree%boxes(p_id), cell_flags)
                call cell_to_ref_flags(cell_flags, tree%n_cell, &
                     tmp_flags(:, thread_id), tree, p_id, ref_buffer)
             end if
          end if
       end do
       !$omp end do
    end do
    !$omp end parallel

    ! Take the highest value over the threads
    do i = 1, size(ref_flags)
       ref_flags(i) = maxval(tmp_flags(i, :))
       if (ref_flags(i) == unset_flag) ref_flags(i) = af_keep_ref
    end do

    if (maxval(ref_flags) > af_do_ref .or. minval(ref_flags) < af_rm_ref) &
         stop "af_adjust_refinement: invalid refinement flag given"

    ! Cannot refine beyond max level
    do i = 1, size(tree%lvls(af_max_lvl)%ids)
       id = tree%lvls(af_max_lvl)%ids(i)
       if (ref_flags(id) == af_do_ref) ref_flags(id) = af_keep_ref
    end do

    call ensure_two_one_balance(tree, ref_flags)
    call handle_derefinement_flags(tree, ref_flags)

    if (present(modify_refinement)) then
       call modify_refinement(tree, ref_flags)
       call ensure_two_one_balance(tree, ref_flags)
       call handle_derefinement_flags(tree, ref_flags)
    end if

  end subroutine consistent_ref_flags

  !> Adjust refinement flags to ensure 2-1 balance is maintained
  subroutine ensure_two_one_balance(tree, ref_flags)
    type(af_t), intent(inout) :: tree
    integer, intent(inout)    :: ref_flags(:)
    integer                   :: lvl, i, id, nb, nb_id
    integer                   :: p_id, p_nb_id

    ! Ensure 2-1 balance
    do lvl = tree%highest_lvl, 1, -1
       do i = 1, size(tree%lvls(lvl)%leaves) ! We only check leaf tree%boxes
          id = tree%lvls(lvl)%leaves(i)

          if (ref_flags(id) == af_do_ref .or. ref_flags(id) == af_refine) then
             ref_flags(id) = af_refine ! Mark for actual refinement

             ! Ensure we will have the necessary neighbors
             do nb = 1, af_num_neighbors
                nb_id = tree%boxes(id)%neighbors(nb)
                if (nb_id == af_no_box) then
                   ! Mark the parent containing neighbor for refinement
                   p_id = tree%boxes(id)%parent
                   p_nb_id = tree%boxes(p_id)%neighbors(nb)
                   ref_flags(p_nb_id) = af_refine
                end if
             end do

          else if (ref_flags(id) == af_rm_ref) then
             ! Ensure we do not remove a required neighbor
             do nb = 1, af_num_neighbors
                nb_id = tree%boxes(id)%neighbors(nb)
                if (nb_id > af_no_box) then
                   if (af_has_children(tree%boxes(nb_id)) .or. &
                        ref_flags(nb_id) > af_keep_ref) then
                      ref_flags(id) = af_keep_ref
                      exit
                   end if
                end if
             end do
          end if

       end do
    end do
  end subroutine ensure_two_one_balance

  subroutine handle_derefinement_flags(tree, ref_flags)
    type(af_t), intent(inout) :: tree
    integer, intent(inout)    :: ref_flags(:)
    integer                   :: lvl, i, id, c_ids(af_num_children)

    ! Make the (de)refinement flags consistent for blocks with children
    do lvl = tree%highest_lvl-1, 1, -1
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          c_ids = tree%boxes(id)%children

          ! Only consider boxes for which at least one child is a leaf
          if (all(af_has_children(tree%boxes(c_ids)))) cycle

          ! Can only remove children if they are all marked for
          ! derefinement, and the box itself not for refinement.
          if (all(ref_flags(c_ids) == af_rm_ref) .and. &
               ref_flags(id) <= af_keep_ref) then
             ref_flags(id) = af_derefine
          else if (ref_flags(id) /= af_derefine) then
             ref_flags(id) = af_keep_ref
             ! The children cannot be removed. This information is useful when
             ! the modify_refinement() routine is used.
             ref_flags(c_ids) = max(ref_flags(c_ids), af_keep_ref)
          end if
       end do
    end do

  end subroutine handle_derefinement_flags

  !> Given the cell refinement flags of a box, set the refinement flag for that
  !> box and potentially also its neighbors (in case of refinement near a
  !> boundary).
  subroutine cell_to_ref_flags(cell_flags, nc, ref_flags, tree, id, &
       ref_buffer)
    use m_af_utils, only: af_get_loc
    integer, intent(in)     :: nc                     !< n_cell for the box
#if NDIM == 2
    integer, intent(in)     :: cell_flags(nc, nc)     !< Cell refinement flags
#elif NDIM == 3
    integer, intent(in)     :: cell_flags(nc, nc, nc) !< Cell refinement flags
#endif
    integer, intent(inout)  :: ref_flags(:)           !< Box refinement flags for this thread
    type(af_t), intent(in) :: tree                   !< Full tree
    integer, intent(in)     :: id                     !< Which box is considered
    integer, intent(in)     :: ref_buffer             !< Buffer cells around refinement
    integer                 :: ix0(NDIM), ix1(NDIM), IJK, nb_id

    if (minval(cell_flags) < af_rm_ref .or. &
         maxval(cell_flags) > af_do_ref) then
       error stop "Error: invalid cell flags given"
    end if

    ! Check whether the box needs to be refined or keep its refinement
    if (any(cell_flags == af_do_ref)) then
       ref_flags(id) = af_do_ref
    else if (any(cell_flags == af_keep_ref)) then
       ref_flags(id) = max(ref_flags(id), af_keep_ref)
    else    ! All flags are af_rm_ref
       ref_flags(id) = max(ref_flags(id), af_rm_ref)
    end if

    if (ref_buffer <= 0) return ! No need to check neighbors

    ! Check whether neighbors also require refinement, which happens when cells
    ! close to the neighbor are flagged.
    do KJI_DO(-1,1)
       if (all([IJK] == 0)) cycle

       nb_id = tree%boxes(id)%neighbor_mat(IJK)

       ! Skip neighbors that are not there
       if (nb_id <= af_no_box) cycle

       ! Compute index range relevant for neighbor
       ix0 = 1
       ix1 = nc
       where ([IJK] == 1)
          ix0 = nc - ref_buffer + 1
          ix1 = nc
       elsewhere ([IJK] == -1)
          ix0 = 1
          ix1 = ref_buffer
       end where

#if NDIM == 2
       if (any(cell_flags(ix0(1):ix1(1), ix0(2):ix1(2)) == af_do_ref)) then
          ref_flags(nb_id) = af_do_ref
       end if
#elif NDIM == 3
       if (any(cell_flags(ix0(1):ix1(1), ix0(2):ix1(2), &
            ix0(3):ix1(3)) == af_do_ref)) then
          ref_flags(nb_id) = af_do_ref
       end if
#endif
    end do; CLOSE_DO

  end subroutine cell_to_ref_flags

  !> Remove the children of box id
  subroutine remove_children(tree, id)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: id !< Id of box whose children will be removed
    integer                   :: ic, c_id, nb_id, nb_rev, nb, IJK, ix

    !$omp critical (crit_remove_children)
    ix = tree%n_removed_ids
    tree%n_removed_ids = tree%n_removed_ids + af_num_children
    !$omp end critical (crit_remove_children)

    tree%removed_ids(ix+1:ix+af_num_children) = tree%boxes(id)%children

    do ic = 1, af_num_children
       c_id = tree%boxes(id)%children(ic)

       ! Remove from neighbors
       do nb = 1, af_num_neighbors
          nb_id = tree%boxes(c_id)%neighbors(nb)
          if (nb_id > af_no_box) then
             nb_rev = af_neighb_rev(nb)
             tree%boxes(nb_id)%neighbors(nb_rev) = af_no_box
          end if
       end do

       do KJI_DO(-1,1)
          nb_id = tree%boxes(c_id)%neighbor_mat(IJK)
          if (nb_id > af_no_box) then
#if NDIM == 2
             tree%boxes(nb_id)%neighbor_mat(-i, -j) = af_no_box
#elif NDIM == 3
             tree%boxes(nb_id)%neighbor_mat(-i, -j, -k) = af_no_box
#endif
          end if
       end do; CLOSE_DO

       tree%boxes(c_id)%in_use = .false.
    end do

    tree%boxes(id)%children = af_no_box
  end subroutine remove_children

  !> Add children to box id, using the indices in c_ids
  subroutine add_children(boxes, id, c_ids, n_cc, n_fc)
    type(box_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id       !< Id of box that gets children
    integer, intent(in)         :: c_ids(af_num_children) !< Free ids for the children
    integer, intent(in)         :: n_cc                   !< Number of cell-centered variables
    integer, intent(in)         :: n_fc                   !< Number of face-centered variables
    integer                     :: i, nb, child_nb(2**(NDIM-1))
    integer                     :: c_id, c_ix_base(NDIM), dix(NDIM)

    boxes(id)%children = c_ids
    c_ix_base          = 2 * boxes(id)%ix - 1

    do i = 1, af_num_children
       c_id                  = c_ids(i)
       boxes(c_id)%ix        = c_ix_base + af_child_dix(:,i)
       boxes(c_id)%lvl       = boxes(id)%lvl+1
       boxes(c_id)%parent    = id
       boxes(c_id)%tag       = af_init_tag
       boxes(c_id)%children  = af_no_box
       boxes(c_id)%neighbors = af_no_box
       boxes(c_id)%neighbor_mat = af_no_box
       boxes(c_id)%neighbor_mat(DTIMES(0)) = c_id
       boxes(c_id)%n_cell    = boxes(id)%n_cell
       boxes(c_id)%coord_t   = boxes(id)%coord_t
       boxes(c_id)%dr        = 0.5_dp * boxes(id)%dr
       boxes(c_id)%r_min     = boxes(id)%r_min + 0.5_dp * boxes(id)%dr * &
            af_child_dix(:,i) * boxes(id)%n_cell

       call af_init_box(boxes(c_id), boxes(id)%n_cell, n_cc, n_fc)
    end do

    ! Set boundary conditions at children
    do nb = 1, af_num_neighbors
       if (boxes(id)%neighbors(nb) < af_no_box) then
          child_nb = c_ids(af_child_adj_nb(:, nb)) ! Neighboring children
          boxes(child_nb)%neighbors(nb) = boxes(id)%neighbors(nb)
          dix = af_neighb_dix(:, nb)
#if NDIM == 2
          boxes(child_nb)%neighbor_mat(dix(1), dix(2)) = &
               boxes(id)%neighbors(nb)
#elif NDIM == 3
          boxes(child_nb)%neighbor_mat(dix(1), dix(2), dix(3)) = &
               boxes(id)%neighbors(nb)
#endif
       end if
    end do
  end subroutine add_children

  !> Create a list c_ids(:) of all the children of p_ids(:). This is used after
  !> a level has been refined.
  subroutine set_child_ids(p_ids, c_ids, boxes)
    integer, intent(in)                 :: p_ids(:) !< All the parents ids
    integer, allocatable, intent(inout) :: c_ids(:) !< Output: all the children's ids
    type(box_t), intent(in)            :: boxes(:) !< List of all the boxes
    integer                             :: i, i0, i1, n_children

    n_children = af_num_children * size(p_ids)
    if (n_children /= size(c_ids)) then
       deallocate(c_ids)
       allocate(c_ids(n_children))
    end if

    do i = 1, size(p_ids)
       i1 = i * af_num_children
       i0 = i1 - af_num_children + 1
       c_ids(i0:i1) = boxes(p_ids(i))%children
    end do
  end subroutine set_child_ids

  !> Restrict fluxes from children to parents on refinement boundaries.
  subroutine af_consistent_fluxes(tree, f_ixs)
    type(af_t), intent(inout)     :: tree         !< Tree to operate on
    integer, intent(in)           :: f_ixs(:)     !< Indices of the fluxes
    integer                       :: lvl, i, id, nb, nb_id

    if (.not. tree%ready) stop "Tree not ready"
    !$omp parallel private(lvl, i, id, nb, nb_id)
    do lvl = 1, tree%highest_lvl-1
       !$omp do
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          do nb = 1, af_num_neighbors
             nb_id = tree%boxes(id)%neighbors(nb)

             ! If the neighbor exists and has no children, set flux
             if (nb_id > af_no_box) then
                if (.not. af_has_children(tree%boxes(nb_id))) then
                   call flux_from_children(tree%boxes, id, nb, f_ixs)
                end if
             end if
          end do
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine af_consistent_fluxes

  !> The neighbor nb has no children and id does, so set flux on the neighbor
  !> from our children. This ensures flux consistency at refinement boundary.
  subroutine flux_from_children(boxes, id, nb, f_ixs)
    type(box_t), intent(inout) :: boxes(:) !< List of all the boxes
    integer, intent(in)         :: id        !< Id of box for which we set fluxes
    integer, intent(in)         :: nb        !< Direction in which fluxes are set
    integer, intent(in)         :: f_ixs(:)  !< Indices of the fluxes
    integer                     :: nc, nch, c_id, i_ch, i, ic, d
    integer                     :: n_chnb, nb_id, i_nb, ioff(NDIM)
#if NDIM == 2
    integer                     :: n
    real(dp)                    :: w1, w2
#endif


    nc     = boxes(id)%n_cell
    nch    = ishft(nc, -1) ! nc/2
    d      = af_neighb_dim(nb)
    n_chnb = 2**(NDIM-1)
    nb_id  = boxes(id)%neighbors(nb)

    if (af_neighb_low(nb)) then
       i = 1
       i_nb = nc+1
    else
       i = nc+1
       i_nb = 1
    end if

    select case (d)
#if NDIM == 2
    case (1)
       do ic = 1, n_chnb
          ! Get index of child adjacent to neighbor
          i_ch = af_child_adj_nb(ic, nb)
          c_id = boxes(id)%children(i_ch)
          ! Index offset of child w.r.t. parent
          ioff = nch*af_child_dix(:, i_ch)
          boxes(nb_id)%fc(i_nb, ioff(2)+1:ioff(2)+nch, 1, f_ixs) = 0.5_dp * ( &
               boxes(c_id)%fc(i, 1:nc:2, 1, f_ixs) + &
               boxes(c_id)%fc(i, 2:nc:2, 1, f_ixs))
       end do
    case (2)
       if (boxes(nb_id)%coord_t == af_cyl) then
          ! In cylindrical symmetry, we take the weighted average
          do ic = 1, n_chnb
             i_ch = af_child_adj_nb(ic, nb)
             c_id = boxes(id)%children(i_ch)
             ioff = nch*af_child_dix(:, i_ch)

             do n = 1, nch
                call af_cyl_child_weights(boxes(nb_id), ioff(1)+n, w1, w2)
                boxes(nb_id)%fc(ioff(1)+n, i_nb, 2, f_ixs) = 0.5_dp * (&
                     w1 * boxes(c_id)%fc(2*n-1, i, 2, f_ixs) + &
                     w2 * boxes(c_id)%fc(2*n, i, 2, f_ixs))
             end do
          end do
       else
          ! Just take the average of the fine fluxes
          do ic = 1, n_chnb
             i_ch = af_child_adj_nb(ic, nb)
             c_id = boxes(id)%children(i_ch)
             ioff = nch*af_child_dix(:, i_ch)
             boxes(nb_id)%fc(ioff(1)+1:ioff(1)+nch, i_nb, 2, f_ixs) = 0.5_dp * ( &
                  boxes(c_id)%fc(1:nc:2, i, 2, f_ixs) + &
                  boxes(c_id)%fc(2:nc:2, i, 2, f_ixs))
          end do
       end if
#elif NDIM == 3
    case (1)
       do ic = 1, n_chnb
          i_ch = af_child_adj_nb(ic, nb)
          c_id = boxes(id)%children(i_ch)
          ioff = nch*af_child_dix(:, i_ch)
          boxes(nb_id)%fc(i_nb, ioff(2)+1:ioff(2)+nch, &
               ioff(3)+1:ioff(3)+nch, 1, f_ixs) = 0.25_dp * ( &
               boxes(c_id)%fc(i, 1:nc:2, 1:nc:2, 1, f_ixs) + &
               boxes(c_id)%fc(i, 2:nc:2, 1:nc:2, 1, f_ixs) + &
               boxes(c_id)%fc(i, 1:nc:2, 2:nc:2, 1, f_ixs) + &
               boxes(c_id)%fc(i, 2:nc:2, 2:nc:2, 1, f_ixs))
       end do
    case (2)
       do ic = 1, n_chnb
          i_ch = af_child_adj_nb(ic, nb)
          c_id = boxes(id)%children(i_ch)
          ioff = nch*af_child_dix(:, i_ch)
          boxes(nb_id)%fc(ioff(1)+1:ioff(1)+nch, i_nb, &
               ioff(3)+1:ioff(3)+nch, 2, f_ixs) = 0.25_dp * ( &
               boxes(c_id)%fc(1:nc:2, i, 1:nc:2, 2, f_ixs) + &
               boxes(c_id)%fc(2:nc:2, i, 1:nc:2, 2, f_ixs) + &
               boxes(c_id)%fc(1:nc:2, i, 2:nc:2, 2, f_ixs) + &
               boxes(c_id)%fc(2:nc:2, i, 2:nc:2, 2, f_ixs))
       end do
    case (3)
       do ic = 1, n_chnb
          i_ch = af_child_adj_nb(ic, nb)
          c_id = boxes(id)%children(i_ch)
          ioff = nch*af_child_dix(:, i_ch)
          boxes(nb_id)%fc(ioff(1)+1:ioff(1)+nch, &
               ioff(2)+1:ioff(2)+nch, i_nb, 3, f_ixs) = 0.25_dp * ( &
               boxes(c_id)%fc(1:nc:2, 1:nc:2, i, 3, f_ixs) + &
               boxes(c_id)%fc(2:nc:2, 1:nc:2, i, 3, f_ixs) + &
               boxes(c_id)%fc(1:nc:2, 2:nc:2, i, 3, f_ixs) + &
               boxes(c_id)%fc(2:nc:2, 2:nc:2, i, 3, f_ixs))
       end do
#endif
    end select
  end subroutine flux_from_children

end module m_af_core
