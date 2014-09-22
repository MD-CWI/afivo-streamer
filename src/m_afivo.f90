module m_afivo

  implicit none
  public

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: a5_max_levels = 20

  ! Tags that can be set for a block
  integer, parameter :: a5_do_ref = 1
  integer, parameter :: a5_rm_ref = 2
  integer, parameter :: a5_kp_ref = 3
  integer, parameter :: a5_rm_children = 4

  ! Neighbor topology information (todo: add documentation)
  integer, parameter :: a5_no_box = 0

  integer, parameter :: nb_lx = 1
  integer, parameter :: nb_hx = 2
  integer, parameter :: nb_ly = 3
  integer, parameter :: nb_hy = 4

  integer, parameter :: a2_ch_dix(2, 4) = reshape((/0,0,1,0,0,1,1,1/), (/2,4/))
  integer, parameter :: a2_ch_rev(4, 2) = reshape((/2,1,4,3,3,4,1,2/), (/4,2/))
  logical, parameter :: a2_ch_low(4, 2) = reshape((/ .true., .false., .true., &
       .false., .true., .true., .false., .false./), (/4,2/))
  logical, parameter :: a2_nb_low(4) = (/.true., .false., .true., .false./)
  integer, parameter :: a2_nb_rev(4) = (/2, 1, 4, 3/)
  integer, parameter :: a2_nb_dim(4) = (/1, 1, 2, 2/)

  ! Each box contains a tag, for which the following bits are set:
  integer, parameter :: a5_bit_in_use = 1
  integer, parameter :: a5_bit_fresh = 2
  integer, parameter :: a5_bit_just_ref = 3

  type box2_t
     integer               :: lvl
     integer               :: tag
     integer               :: ix(2)
     integer               :: parent
     integer               :: children(4)
     integer               :: neighbors(4)
     real(dp)              :: dr(2)
     real(dp)              :: r_min(2)
     real(dp), allocatable :: cc(:, :, :)
     real(dp), allocatable :: fx(:, :, :)
     real(dp), allocatable :: fy(:, :, :)
  end type box2_t

  type level2_t
     integer, allocatable :: ids(:)
  end type level2_t

  type a2_t
     type(level2_t)            :: levels(a5_max_levels)
     integer                   :: box_cells
     integer                   :: n_cc
     integer                   :: n_fx
     integer                   :: n_fy
     integer                   :: n_boxes
     type(box2_t), allocatable :: boxes(:)
  end type a2_t

  abstract interface
     integer function a2_to_int_f(box)
       import
       type(box2_t), intent(in) :: box
     end function a2_to_int_f

     subroutine a2_subr(box)
       import
       type(box2_t), intent(inout) :: box
     end subroutine a2_subr

     subroutine a2_subr_boxes(boxes, id)
       import
       type(box2_t), intent(inout) :: boxes(:)
       integer, intent(in)         :: id
     end subroutine a2_subr_boxes
  end interface

contains

  subroutine a2_init(self, n_boxes_max, box_cells, n_cc, n_fx, n_fy)
    type(a2_t)          :: self
    integer, intent(in) :: n_boxes_max
    integer, intent(in) :: box_cells, n_cc, n_fx, n_fy
    integer             :: lvl

    self%box_cells = box_cells
    self%n_cc     = n_cc
    self%n_fx     = n_fx
    self%n_fy     = n_fy
    self%n_boxes  = 0
    allocate(self%boxes(n_boxes_max))
    do lvl = 1, a5_max_levels
       allocate(self%levels(lvl)%ids(0))
    end do
  end subroutine a2_init

  integer function a5_n_levels(self)
    type(a2_t), intent(in) :: self
    integer                :: lvl
    do lvl = 1, a5_max_levels-1
       if (size(self%levels(lvl)%ids) == 0) exit
    end do
    a5_n_levels = lvl
  end function a5_n_levels

  subroutine a2_loop_box(self, my_procedure)
    type(a2_t), intent(inout) :: self
    procedure(a2_subr)        :: my_procedure
    integer                   :: lvl, i, id

    do lvl = 1, a5_n_levels(self)
       do i = 1, size(self%levels(lvl)%ids)
          id = self%levels(lvl)%ids(i)
          call my_procedure(self%boxes(id))
       end do
    end do
  end subroutine a2_loop_box

  subroutine a2_loop_boxes(self, my_procedure)
    type(a2_t), intent(inout) :: self
    procedure(a2_subr_boxes)  :: my_procedure
    integer                   :: lvl, i, id

    do lvl = 1, a5_n_levels(self)
       do i = 1, size(self%levels(lvl)%ids)
          id = self%levels(lvl)%ids(i)
          call my_procedure(self%boxes, id)
       end do
    end do
  end subroutine a2_loop_boxes

  subroutine a5_tidy_storage(self, new_size)
    type(a2_t) :: self
    integer, intent(in) :: new_size

    ! TODO
  end subroutine a5_tidy_storage

  subroutine a2_set_neighbors(boxes, ids)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:)
    integer                     :: i, id

    do i = 1, size(ids)
       ! For each "fresh" box, find possible neighbors
       id = ids(i)
       if (btest(boxes(id)%tag, a5_bit_fresh)) then
          call set_nbs_2d(boxes, id)
       end if
    end do
  end subroutine a2_set_neighbors

  subroutine set_nbs_2d(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)          :: id
    integer                      :: nb, nb_id

    do nb = 1, 4                 ! Dimension
       if (boxes(id)%neighbors(nb) == a5_no_box) then
          nb_id = find_nb_2d(boxes, id, nb)
          if (nb_id /= a5_no_box) then
             boxes(id)%neighbors(nb) = nb_id
             boxes(nb_id)%neighbors(a2_nb_rev(nb)) = id
          end if
       end if
    end do
  end subroutine set_nbs_2d

  integer function a2_ix_to_cix(ix)
    integer, intent(in) :: ix(2)
    ! Second index odd: -2, first ix odd: -1
    a2_ix_to_cix = 4 - 2 * iand(ix(2), 1) - iand(ix(1), 1)
  end function a2_ix_to_cix

  ! Get neighbor nb of id, through its parent
  function find_nb_2d(boxes, id, nb) result(nb_id)
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id, nb
    integer                  :: nb_id, p_id, c_ix, d, old_pid

    p_id = boxes(id)%parent
    old_pid = p_id
    c_ix = a2_ix_to_cix(boxes(id)%ix)
    d = a2_nb_dim(nb)

    ! Check if neighbor is in same direction as ix is (low/high). If so,
    ! use neighbor of parent
    if (a2_ch_low(c_ix, d) .eqv. a2_nb_low(nb)) &
         p_id = boxes(p_id)%neighbors(nb)

    ! The child ix of the neighbor is reversed in direction d
    nb_id = boxes(p_id)%children(a2_ch_rev(c_ix, d))
  end function find_nb_2d

  ! After boxes have been added to level 1, this stores them as a "level"
  subroutine a2_set_base(self, ix_list, nb_list, dr, r_min_11)
    type(a2_t), intent(inout) :: self
    integer, intent(in)       :: ix_list(:, :), nb_list(:, :)
    real(dp), intent(in)      :: dr(2), r_min_11(2)
    integer                   :: n_boxes, i, id

    n_boxes       = size(ix_list, 2)
    deallocate(self%levels(1)%ids)
    allocate(self%levels(1)%ids(n_boxes))
    call a2_get_free_ids(self, self%levels(1)%ids)

    do i = 1, n_boxes
       id = self%levels(1)%ids(i)
       self%boxes(id)%ix        = ix_list(:, i)
       self%boxes(id)%lvl       = 1
       self%boxes(id)%parent    = 0
       self%boxes(id)%children  = 0
       self%boxes(id)%neighbors = nb_list(:, i)
       self%boxes(id)%dr        = dr
       self%boxes(id)%r_min     = r_min_11 + (self%boxes(id)%ix - (/1,1/)) * dr

       ! Set "in_use" and "fresh" tag
       self%boxes(id)%tag       = ibset(0, a5_bit_in_use)
       self%boxes(id)%tag       = ibset(self%boxes(id)%tag, a5_bit_fresh)
    end do
  end subroutine a2_set_base

  ! On input, self should be balanced. On output, self is still balanced, and
  ! its refinement is updated (with at most one level per call).
  ! Sets the following bits: a5_bit_fresh, a5_bit_just_ref
  subroutine a2_adjust_refinement(self, ref_func)
    type(a2_t), intent(inout) :: self
    procedure(a2_to_int_f)    :: ref_func
    integer                   :: lvl, id, i, c_ids(4), n_boxes_start
    integer, allocatable      :: ref_vals(:)

    n_boxes_start = self%n_boxes
    allocate(ref_vals(n_boxes_start))
    call a2_set_ref_vals(self%levels, self%boxes, &
         ref_vals, ref_func)

    do lvl = 1, a5_n_levels(self)
       do i = 1, size(self%levels(lvl)%ids)
          id = self%levels(lvl)%ids(i)

          if (id > n_boxes_start) then
             cycle              ! This is a newly added box
          else if (ref_vals(id) == a5_do_ref) then
             ! Add children. First need to get 4 free id's
             call a2_get_free_ids(self, c_ids)
             call a2_add_children(self%boxes, id, c_ids, self%box_cells)
          else if (ref_vals(id) == a5_rm_children) then
             ! Remove children
             call a2_remove_children(self%boxes, id)
          end if
       end do

       ! Set the next level
       if (lvl < a5_max_levels) then
          deallocate(self%levels(lvl+1)%ids)
          call set_child_lvl(self%levels(lvl)%ids, &
               self%levels(lvl+1)%ids, self%boxes)
          ! Set connectivity
          call a2_set_neighbors(self%boxes, self%levels(lvl+1)%ids)
       end if
    end do
  end subroutine a2_adjust_refinement

  subroutine a2_get_free_ids(self, ids)
    type(a2_t), intent(inout) :: self
    integer, intent(out)      :: ids(:)
    integer                   :: i, n_ids, id, bs

    n_ids = size(ids)
    ! Critical part
    ids(1)       = self%n_boxes + 1
    self%n_boxes = self%n_boxes + n_ids
    ! End critical
    do i = 2, n_ids
       ids(i) = ids(1) + i - 1
    end do

    ! Allocate storage
    bs = self%box_cells
    do i = 1, n_ids
       id = ids(i)
       allocate(self%boxes(id)%cc(0:bs+1, 0:bs+1, self%n_cc))
       allocate(self%boxes(id)%fx(bs+1, bs, self%n_fx))
       allocate(self%boxes(id)%fy(bs, bs+1, self%n_fy))
    end do
  end subroutine a2_get_free_ids

  subroutine a2_set_ref_vals(levels, boxes, ref_vals, ref_func)
    type(level2_t), intent(in)  :: levels(:)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(inout)      :: ref_vals(:)
    procedure(a2_to_int_f)      :: ref_func

    integer                     :: n_levels
    integer                     :: lvl, i, id, c_ids(4)
    integer                     :: nb, p_id, nb_id, p_nb_id

    n_levels = size(levels)

    ! Set refinement flags for all boxes using ref_func
    do lvl = 1, n_levels
       do i = 1, size(levels(lvl)%ids)
          id           = levels(lvl)%ids(i)
          ref_vals(id) = ref_func(boxes(id))
       end do
    end do

    ! Cannot derefine lvl 1
    do i = 1, size(levels(1)%ids)
       id = levels(1)%ids(i)
       if (ref_vals(id) == a5_rm_ref) ref_vals(id) = a5_kp_ref
    end do

    ! Cannot refine beyond max level
    do i = 1, size(levels(a5_max_levels)%ids)
       id = levels(a5_max_levels)%ids(i)
       if (ref_vals(id) == a5_do_ref) ref_vals(id) = a5_kp_ref
    end do

    ! Ensure 2-1 balance.
    do lvl = n_levels, 2, -1
       do i = 1, size(levels(lvl)%ids)
          id = levels(lvl)%ids(i)
          ! We only need check leaf boxes
          if (a2_has_children(boxes(id))) cycle

          if (ref_vals(id) == a5_do_ref) then
             ! Ensure we will have the necessary neighbors
             do nb = 1, 4
                nb_id = boxes(id)%neighbors(nb)
                if (nb_id == a5_no_box) then
                   ! Mark the parent containing neighbor for refinement
                   p_id = boxes(id)%parent
                   p_nb_id = boxes(p_id)%neighbors(nb)
                   ref_vals(p_nb_id) = a5_do_ref
                end if
             end do
          else if (ref_vals(id) == a5_rm_ref) then
             ! Ensure we do not remove a required neighbor
             do nb = 1, 4
                nb_id = boxes(id)%neighbors(nb)
                if (nb_id /= a5_no_box) then
                   if (a2_has_children(boxes(nb_id)) .or. &
                        ref_vals(nb_id) == a5_do_ref) then
                      ref_vals(id) = a5_kp_ref
                   end if
                end if
             end do
          end if
       end do
    end do

    ! Make the (de)refinement flags consistent for blocks with children.
    do lvl = n_levels-1, 1, -1
       do i = 1, size(levels(lvl)%ids)
          id = levels(lvl)%ids(i)
          if (.not. a2_has_children(boxes(id))) cycle

          ! Can only remove children if they are all marked for
          ! derefinement, and the box itself not for refinement.
          c_ids = boxes(id)%children
          if (all(ref_vals(c_ids) == a5_rm_ref) .and. &
               ref_vals(id) /= a5_do_ref) then
             ! The children are removed
             ref_vals(id) = a5_rm_children
          else
             ! The children cannot be removed
             ref_vals(id) = a5_kp_ref
             where (ref_vals(c_ids) == a5_rm_ref)
                ref_vals(c_ids) = a5_kp_ref
             end where
          end if
       end do
    end do

  end subroutine a2_set_ref_vals

  subroutine a2_remove_children(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    integer                     :: c_ids(4)
    integer                     :: i, c_id, nb_id, nb_rev, nb

    c_ids              = boxes(id)%children
    boxes(id)%children = a5_no_box

    do i = 1, 4
       c_id            = c_ids(i)
       boxes(c_id)%tag = ibclr(boxes(c_id)%tag, a5_bit_in_use)

       ! Remove this box at the neighbors
       do nb = 1, 4
          nb_id = boxes(c_id)%neighbors(nb)
          if (nb_id /= a5_no_box) then
             nb_rev = a2_nb_rev(nb)
             boxes(nb_id)%neighbors(nb_rev) = a5_no_box
          end if
       end do
    end do
  end subroutine a2_remove_children

  subroutine a2_add_children(boxes, id, c_ids, box_cells)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, box_cells
    integer, intent(in)         :: c_ids(4)
    integer                     :: i, c_id, c_ix_base(2)

    boxes(id)%children = c_ids
    c_ix_base          = 2 * boxes(id)%ix - 1
    boxes(id)%tag      = ibset(boxes(id)%tag, a5_bit_just_ref)

    do i = 1, 4
       c_id                  = c_ids(i)
       boxes(c_id)%ix        = c_ix_base + a2_ch_dix(:,i)
       boxes(c_id)%lvl       = boxes(id)%lvl+1
       boxes(c_id)%dr        = boxes(id)%dr * 0.5_dp
       boxes(c_id)%r_min     = boxes(id)%r_min + &
            a2_ch_dix(:,i) * boxes(c_id)%dr * box_cells
       boxes(c_id)%parent    = id
       boxes(c_id)%children  = a5_no_box
       boxes(c_id)%neighbors = a5_no_box

       ! Set "in_use" and "fresh" tag
       boxes(c_id)%tag       = ibset(0, a5_bit_in_use)
       boxes(c_id)%tag       = ibset(boxes(c_id)%tag, a5_bit_fresh)
    end do
  end subroutine a2_add_children

  subroutine set_child_lvl(p_ids, c_ids, boxes)
    integer, intent(in)                 :: p_ids(:)
    integer, allocatable, intent(inout) :: c_ids(:)
    type(box2_t), intent(in)            :: boxes(:)
    integer                             :: i, ip, ic, n_children

    ! Count 4 times the number of refined parent blocks
    n_children = 4 * count(a2_has_children(boxes(p_ids)))
    allocate(c_ids(n_children))

    ic = 0
    do i = 1, size(p_ids)
       ip = p_ids(i)
       if (a2_has_children(boxes(ip))) then
          c_ids(ic+1:ic+4) = boxes(ip)%children
          ic = ic + 4
       end if
    end do
  end subroutine set_child_lvl

  elemental logical function a2_has_children(box)
    type(box2_t), intent(in) :: box
    a2_has_children = (box%children(1) /= a5_no_box)
  end function a2_has_children

  ! Prolongation using bilinear interpolation. Uses ghost cells and corners.
  subroutine a2_prolong_from(boxes, id, v_ixs)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, v_ixs(:)
    real(dp), parameter         :: one16th=1/16.0_dp
    integer                     :: bs, ic, c_id, i_base, j_base
    integer                     :: i, j, i_c1, i_c2, j_c1, j_c2

    bs = size(boxes(id)%cc, 1) - 2
    do ic = 1, 4
       c_id = boxes(id)%children(ic)

       ! Offset of child w.r.g. parent
       i_base = a2_ch_dix(1, ic) * bs/2
       j_base = a2_ch_dix(2, ic) * bs/2

       ! In these loops, we calculate the closest coarse index (i_c1, j_c1), and
       ! the one-but-closest (i_c2, j_c2). The fine cell lies in between.
       do j = 1, bs
          j_c1 = j_base + ishft(j, -1) + iand(j, 1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          do i = 1, bs
             i_c1 = i_base + ishft(i, -1) + iand(i, 1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1

             boxes(c_id)%cc(i, j, v_ixs) = one16th * ( &
                  9*boxes(id)%cc(i_c1, j_c1, v_ixs) + &
                  3*boxes(id)%cc(i_c2, j_c1, v_ixs) + &
                  3*boxes(id)%cc(i_c1, j_c2, v_ixs) + &
                  boxes(id)%cc(i_c2, j_c2, v_ixs))
          end do
       end do
    end do

  end subroutine a2_prolong_from

  ! subroutine a2_restrict_to(boxes, id)

  ! end subroutine a2_restrict

  subroutine a2_fill_internal_gc(self)
    type(a2_t), intent(inout) :: self
    integer :: lvl, i, id
    do lvl = 1, a5_n_levels(self)
       do i = 1, size(self%levels(lvl)%ids)
          id = self%levels(lvl)%ids(i)
          call a2_internal_gc(self%boxes, id)
       end do
    end do
  end subroutine a2_fill_internal_gc

  subroutine a2_internal_gc(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in) :: id
    integer :: bs, nb_id

    bs = size(boxes(id)%cc, 1) - 2
    nb_id = boxes(id)%neighbors(1)
    if (nb_id > a5_no_box) &
         boxes(id)%cc(0, 1:bs, :) = boxes(nb_id)%cc(bs, 1:bs, :)
    nb_id = boxes(id)%neighbors(2)
    if (nb_id > a5_no_box) &
         boxes(id)%cc(bs+1, 1:bs, :) = boxes(nb_id)%cc(1, 1:bs, :)
    nb_id = boxes(id)%neighbors(3)
    if (nb_id > a5_no_box) &
         boxes(id)%cc(1:bs, 0, :) = boxes(nb_id)%cc(1:bs, bs, :)
    nb_id = boxes(id)%neighbors(4)
    if (nb_id > a5_no_box) &
         boxes(id)%cc(1:bs, bs+1, :) = boxes(nb_id)%cc(1:bs, 1, :)
  end subroutine a2_internal_gc

  subroutine a2_write_tree(self, filename, cc_names, cc_units, n_cycle, time)
    use m_vtk
    type(a2_t), intent(in) :: self
    character(len=*)       :: filename, cc_names(:), cc_units(:)
    integer, intent(in)    :: n_cycle
    real(dp), intent(in)   :: time
    integer                :: lvl, bc, bn, n, n_cells, n_nodes, n_grids
    integer                :: ig, i, j, id, iv, n_ix, c_ix
    integer                :: cell_ix, node_ix
    integer                :: nodes_per_box, cells_per_box
    real(dp), allocatable  :: coords(:), v(:)
    integer, allocatable   :: offsets(:), connects(:), cell_types(:)
    type(vtk_t)            :: vtkf

    bc = self%box_cells
    bn = bc + 1
    nodes_per_box = bn**2
    cells_per_box = bc**2

    n_grids = 0
    do lvl = 1, a5_n_levels(self)
       do n = 1, size(self%levels(lvl)%ids)
          id = self%levels(lvl)%ids(n)
          if (.not. a2_has_children(self%boxes(id))) n_grids = n_grids + 1
       end do
    end do
    n_nodes = nodes_per_box * n_grids
    n_cells = cells_per_box * n_grids

    allocate(coords(2 * n_nodes))
    allocate(v(n_cells))
    allocate(offsets(cells_per_box * n_grids))
    allocate(cell_types(cells_per_box * n_grids))
    allocate(connects(4 * cells_per_box * n_grids))

    cell_types = 8              ! VTK pixel type

    ig = 0
    do lvl = 1, a5_n_levels(self)
       do n = 1, size(self%levels(lvl)%ids)
          id = self%levels(lvl)%ids(n)
          if (a2_has_children(self%boxes(id))) cycle

          ig = ig + 1
          cell_ix = (ig-1) * cells_per_box
          node_ix = (ig-1) * nodes_per_box

          do j = 1, bn
             do i = 1, bn
                n_ix = 2 * (node_ix + (j-1) * bn + i) - 1
                coords(n_ix) = self%boxes(id)%r_min(1) + &
                     (i-1) * self%boxes(id)%dr(1)
                coords(n_ix+1) = self%boxes(id)%r_min(2) + &
                     (j-1) * self%boxes(id)%dr(2)
             end do
          end do

          do j = 1, bc
             do i = 1, bc
                ! In vtk, indexing starts at 0, so subtract 1
                n_ix                      = node_ix + (j-1) * bn + i - 1
                c_ix                      = cell_ix + (j-1) * bc + i
                v(c_ix) = self%boxes(id)%cc(i, j, 1)
                offsets(c_ix)             = 4 * c_ix
                connects(4*c_ix-3:4*c_ix) = (/n_ix, n_ix+1, n_ix+bn, n_ix+bn+1/)
             end do
          end do
       end do
    end do

    call vtk_ini_xml(vtkf, trim(filename), 'UnstructuredGrid')
    call vtk_dat_xml(vtkf, "UnstructuredGrid", .true.)
    call vtk_geo_xml(vtkf, coords, n_nodes, n_cells, 2, n_cycle, time)
    call vtk_con_xml(vtkf, connects, offsets, cell_types, n_cells)
    call vtk_dat_xml(vtkf, "CellData", .true.)
    call vtk_var_r8_xml(vtkf, 'scalars', v, n_cells)
    call vtk_dat_xml(vtkf, "CellData", .false.)
    call vtk_geo_xml_close(vtkf)
    call vtk_dat_xml(vtkf, "UnstructuredGrid", .false.)
    call vtk_end_xml(vtkf)
    print *, "Written ", trim(filename), "n_grids", n_grids
  end subroutine a2_write_tree

  ! subroutine a2_write_tree(self, filename, cc_names, cc_units, n_cycle, time)
  !   use m_write_silo
  !   type(a2_t), intent(in)          :: self
  !   character(len=*)                :: filename, cc_names(:), cc_units(:)
  !   integer, intent(in)             :: n_cycle
  !   real(dp), intent(in)            :: time
  !   character(len=*), parameter     :: grid_name = "gg", var_name  = "vv"
  !   character(len=*), parameter     :: amr_name  = "amr"
  !   character(len=100), allocatable :: grid_list(:), var_list(:, :)
  !   integer                         :: lvl, i, id, ig, iv, bs, n_grids, dbix

  !   bs = self%box_cells
  !   n_grids = 0
  !   do lvl = 1, a5_n_levels(self)
  !      n_grids = n_grids + size(self%levels(lvl)%ids)
  !   end do

  !   allocate(grid_list(n_grids))
  !   allocate(var_list(self%n_cc, n_grids))

  !   call SILO_create_file(filename, dbix)
  !   ig = 0

  !   do lvl = 1, a5_n_levels(self)
  !      do i = 1, size(self%levels(lvl)%ids)
  !         id = self%levels(lvl)%ids(i)
  !         ig = ig + 1
  !         write(grid_list(ig), "(A,I0)") grid_name, ig
  !         call SILO_add_grid(dbix, grid_list(ig), 2, &
  !              (/bs+1, bs+1/), self%boxes(id)%r_min, self%boxes(id)%dr)
  !         print *, id, self%boxes(id)%r_min, self%boxes(id)%dr
  !         do iv = 1, self%n_cc
  !            write(var_list(iv, ig), "(A,I0)") trim(cc_names(iv)) // "_", ig
  !            call SILO_add_var(dbix, var_list(iv, ig), grid_list(ig), &
  !                 pack(self%boxes(id)%cc(1:bs, 1:bs, iv), .true.), (/bs, bs/), &
  !                 trim(cc_units(iv)))
  !         end do
  !      end do
  !   end do

  !   call SILO_set_mmesh_grid(dbix, amr_name, grid_list, n_cycle, time)
  !   do iv = 1, self%n_cc
  !      call SILO_set_mmesh_var(dbix, trim(cc_names(iv)), amr_name, &
  !           var_list(iv, :), n_cycle, time)
  !   end do
  !   call SILO_close_file(dbix)
  !   print *, "Number of grids", ig
  ! end subroutine a2_write_tree

end module m_afivo