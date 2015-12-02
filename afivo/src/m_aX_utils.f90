module m_a$D_utils
  use m_a$D_core

  implicit none
  private

  ! Public subroutines
  public :: a$D_loop_box
  public :: a$D_loop_box_arg
  public :: a$D_loop_boxes
  public :: a$D_loop_boxes_arg
  public :: a$D_box_clear_cc
  public :: a$D_box_add_cc
  public :: a$D_box_sub_cc
  public :: a$D_box_times_cc
  public :: a$D_box_lincomb_cc
  public :: a$D_box_copy_cc_to
  public :: a$D_box_copy_cc
  public :: a$D_boxes_copy_cc
  public :: a$D_tree_copy_cc
  public :: a$D_reduction
  public :: a$D_reduction_loc
  public :: a$D_tree_max_cc
  public :: a$D_tree_min_cc
  public :: a$D_tree_sum_cc
  public :: a$D_box_copy_fc
  public :: a$D_boxes_copy_fc
  public :: a$D_tree_copy_fc

  ! Public functions
  public :: a$D_lvl_dr
  public :: a$D_min_dr
  public :: a$D_r_inside
  public :: a$D_r_center
  public :: a$D_get_child_offset
  public :: a$D_get_loc
  public :: child_that_contains
  public :: a$D_cc_ix
  public :: a$D_r_cc
  public :: a$D_r_loc
#if $D == 2
  public :: a$D_cyl_radius_cc
#endif
  public :: a$D_rr_cc
  public :: a$D_r_node
  public :: a$D_n_cell

contains

  !> Call procedure for each box in tree
  subroutine a$D_loop_box(tree, my_procedure, leaves_only)
    type(a$D_t), intent(inout)     :: tree
    procedure(a$D_subr)           :: my_procedure
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                       :: lvl, i, id

    leaves = .false.; if (present(leaves_only)) leaves = leaves_only
    if (.not. tree%ready) stop "a$D_loop_box: set_base has not been called"

    !$omp parallel private(lvl, i, id)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       if (leaves) then
          !$omp do
          do i = 1, size(tree%lvls(lvl)%leaves)
             id = tree%lvls(lvl)%leaves(i)
             call my_procedure(tree%boxes(id))
          end do
          !$omp end do
       else
          !$omp do
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             call my_procedure(tree%boxes(id))
          end do
          !$omp end do
       end if
    end do
    !$omp end parallel
  end subroutine a$D_loop_box

  !> Call procedure for each box in tree, with argument rarg
  subroutine a$D_loop_box_arg(tree, my_procedure, rarg, leaves_only)
    type(a$D_t), intent(inout)     :: tree
    procedure(a$D_subr_arg)        :: my_procedure
    real(dp), intent(in)          :: rarg(:)
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                       :: lvl, i, id

    if (.not. tree%ready) stop "Tree not ready"
    leaves = .false.; if (present(leaves_only)) leaves = leaves_only

    !$omp parallel private(lvl, i, id)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       if (leaves) then
          !$omp do
          do i = 1, size(tree%lvls(lvl)%leaves)
             id = tree%lvls(lvl)%leaves(i)
             call my_procedure(tree%boxes(id), rarg)
          end do
          !$omp end do
       else
          !$omp do
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             call my_procedure(tree%boxes(id), rarg)
          end do
          !$omp end do
       end if
    end do
    !$omp end parallel
  end subroutine a$D_loop_box_arg

  !> Call procedure for each id in tree, giving the list of boxes
  subroutine a$D_loop_boxes(tree, my_procedure, leaves_only)
    type(a$D_t), intent(inout)     :: tree
    procedure(a$D_subr_boxes)      :: my_procedure
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                       :: lvl, i, id

    if (.not. tree%ready) stop "Tree not ready"
    leaves = .false.; if (present(leaves_only)) leaves = leaves_only

    !$omp parallel private(lvl, i, id)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       if (leaves) then
          !$omp do
          do i = 1, size(tree%lvls(lvl)%leaves)
             id = tree%lvls(lvl)%leaves(i)
             call my_procedure(tree%boxes, id)
          end do
          !$omp end do
       else
          !$omp do
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             call my_procedure(tree%boxes, id)
          end do
          !$omp end do
       end if
    end do
    !$omp end parallel
  end subroutine a$D_loop_boxes

  !> Call procedure for each id in tree, giving the list of boxes
  subroutine a$D_loop_boxes_arg(tree, my_procedure, rarg, leaves_only)
    type(a$D_t), intent(inout)    :: tree
    procedure(a$D_subr_boxes_arg) :: my_procedure
    real(dp), intent(in)         :: rarg(:)
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                      :: lvl, i, id

    if (.not. tree%ready) stop "Tree not ready"
    leaves = .false.; if (present(leaves_only)) leaves = leaves_only

    !$omp parallel private(lvl, i, id)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       if (leaves) then
          !$omp do
          do i = 1, size(tree%lvls(lvl)%leaves)
             id = tree%lvls(lvl)%leaves(i)
             call my_procedure(tree%boxes, id, rarg)
          end do
          !$omp end do
       else
          !$omp do
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             call my_procedure(tree%boxes, id, rarg)
          end do
          !$omp end do
       end if
    end do
    !$omp end parallel
  end subroutine a$D_loop_boxes_arg

  !> Return dr at lvl
  pure function a$D_lvl_dr(tree, lvl) result(dr)
    type(a$D_t), intent(in) :: tree
    integer, intent(in)    :: lvl
    real(dp)               :: dr !< Output: dr at the finest lvl of the tree
    dr = tree%dr_base * 0.5_dp**(lvl-1)
  end function a$D_lvl_dr

  !> Return finest dr that is used in the tree
  pure function a$D_min_dr(tree) result(dr)
    type(a$D_t), intent(in) :: tree
    real(dp)               :: dr !< Output: dr at the finest lvl of the tree
    dr = a$D_lvl_dr(tree, tree%max_lvl)
  end function a$D_min_dr

  !> Returns whether r is inside or within a distance d from box
  pure function a$D_r_inside(box, r, d) result(inside)
    type(box$D_t), intent(in)       :: box
    real(dp), intent(in)           :: r($D)
    real(dp), intent(in), optional :: d
    real(dp)                       :: r_max($D)
    logical                        :: inside

    r_max = box%r_min + box%dr * box%n_cell
    if (present(d)) then
       inside = all(r+d >= box%r_min) .and. all(r-d <= r_max)
    else
       inside = all(r >= box%r_min) .and. all(r <= r_max)
    end if
  end function a$D_r_inside

  !> Return the coordinate of the center of a box
  pure function a$D_r_center(box) result(r_center)
    type(box$D_t), intent(in) :: box
    real(dp)                 :: r_center($D)
    r_center = box%r_min + 0.5_dp * box%n_cell * box%dr
  end function a$D_r_center

  !> Get the offset of a box with respect to its parent (e.g. in 2d, there can
  !> be a child at offset 0,0, one at n_cell/2,0, one at 0,n_cell/2 etc.)
  function a$D_get_child_offset(box, nb) result(ix_offset)
    type(box$D_t), intent(in)           :: box   !< A child box
    integer, intent(in), optional      :: nb     !< Optional: get index on parent neighbor
    integer                            :: ix_offset($D)
    if (box%lvl > 1) then
       ix_offset = iand(box%ix-1, 1) * ishft(box%n_cell, -1) ! * n_cell / 2
       if (present(nb)) ix_offset = ix_offset - a$D_nb_dix(:, nb) * box%n_cell
    else                        ! In the subtree, parents are half the size
       ix_offset = 0
       if (present(nb)) ix_offset = ix_offset - &
            a$D_nb_dix(:, nb) * ishft(box%n_cell, -1) ! n_cell / 2
    endif
  end function a$D_get_child_offset

  !> Get the location of the finest cell containing rr. If max_lvl is present,
  !> do not go to a finer level than max_lvl. If there is no box containing rr,
  !> return a location of -1
  pure function a$D_get_loc(tree, rr, max_lvl) result(loc)
    type(a$D_t), intent(in)       :: tree   !< Tree
    real(dp), intent(in)          :: rr($D) !< Coordinate
    integer, intent(in), optional :: max_lvl !< Maximum level of box
    type(a$D_loc_t)               :: loc    !< Location of cell

    integer :: i, id, i_ch, lvl_max

    lvl_max = tree%lvls_max
    if (present(max_lvl)) lvl_max = max_lvl

    ! Find lvl 1 box that includes rr
    do i = 1, size(tree%lvls(1)%ids)
       id = tree%lvls(1)%ids(i)
       if (a$D_r_inside(tree%boxes(id), rr)) exit
    end do

    ! If not inside any box, return
    if (i > size(tree%lvls(1)%ids)) then
       loc%id = -1
       loc%ix = -1
       return
    end if

    ! Jump into children for as long as possible
    do
       if (tree%boxes(id)%lvl < lvl_max .and. &
            a$D_has_children(tree%boxes(id))) then
          i_ch = child_that_contains(tree%boxes(id), rr)
          id = tree%boxes(id)%children(i_ch)
       else
          exit
       end if
    end do

    loc%id = id
    loc%ix = a$D_cc_ix(tree%boxes(id), rr)
  end function a$D_get_loc

  !> For a box with children that contains rr, find in which child rr lies
  pure function child_that_contains(box, rr) result(i_ch)
    type(box$D_t), intent(in) :: box    !< A box with children
    real(dp), intent(in)      :: rr($D) !< Location inside the box
    integer                   :: i_ch   !< Index of child containing rr
    real(dp)                  :: cntr($D)

    i_ch = 1
    cntr = box%r_min + box%dr * ishft(box%n_cell, -1)

    if (rr(1) > cntr(1)) i_ch = i_ch + 1
    if (rr(2) > cntr(2)) i_ch = i_ch + 2
#if $D==3
    if (rr(3) > cntr(3)) i_ch = i_ch + 4
#endif
  end function child_that_contains

  !> Get the index of the cell that includes point r
  pure function a$D_cc_ix(box, r) result(cc_ix)
    type(box$D_t), intent(in) :: box
    real(dp), intent(in)     :: r($D)
    integer                  :: cc_ix($D)
    cc_ix = ceiling((r - box%r_min) / box%dr)
  end function a$D_cc_ix

  !> Get the location of the cell center with index cc_ix
  pure function a$D_r_cc(box, cc_ix) result(r)
    type(box$D_t), intent(in) :: box
    integer, intent(in)      :: cc_ix($D)
    real(dp)                 :: r($D)
    r = box%r_min + (cc_ix-0.5_dp) * box%dr
  end function a$D_r_cc

  !> Get the location of "loc"
  pure function a$D_r_loc(tree, loc) result(r)
    type(a$D_t), intent(in)     :: tree
    type(a$D_loc_t), intent(in) :: loc
    real(dp)                   :: r($D)
    r = tree%boxes(loc%id)%r_min + &
         (loc%ix-0.5_dp) * tree%boxes(loc%id)%dr
  end function a$D_r_loc

#if $D == 2
  !> Get the radius of the cell center with index cc_ix
  pure function a$D_cyl_radius_cc(box, cc_ix) result(r)
    type(box$D_t), intent(in) :: box
    integer, intent(in)      :: cc_ix($D)
    real(dp)                 :: r
    r = box%r_min(1) + (cc_ix(1)-0.5_dp) * box%dr
  end function a$D_cyl_radius_cc
#endif

  !> Get a general location with real index cc_ix (like a$D_r_cc).
  pure function a$D_rr_cc(box, cc_ix) result(r)
    type(box$D_t), intent(in) :: box
    real(dp), intent(in)     :: cc_ix($D)
    real(dp)                 :: r($D)
    r = box%r_min + (cc_ix-0.5_dp) * box%dr
  end function a$D_rr_cc

  !> Get the location of a node (cell-corner) with index nd_ix
  pure function a$D_r_node(box, nd_ix) result(r)
    type(box$D_t), intent(in) :: box
    integer, intent(in)      :: nd_ix($D)
    real(dp)                 :: r($D)
    r = box%r_min + (nd_ix-1) * box%dr
  end function a$D_r_node

  !> Set cc(..., iv) = 0
  subroutine a$D_box_clear_cc(box, iv)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: iv
#if $D == 2
    box%cc(:,:, iv) = 0
#elif $D == 3
    box%cc(:,:,:, iv) = 0
#endif
  end subroutine a$D_box_clear_cc

  !> Add cc(..., iv_from) to box%cc(..., iv_to)
  subroutine a$D_box_add_cc(box, iv_from, iv_to)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: iv_from, iv_to
#if $D == 2
    box%cc(:,:, iv_to) = box%cc(:,:, iv_to) + box%cc(:,:, iv_from)
#elif $D == 3
    box%cc(:,:,:, iv_to) = box%cc(:,:,:, iv_to) + box%cc(:,:,:, iv_from)
#endif
  end subroutine a$D_box_add_cc

  !> Subtract cc(..., iv_from) from box%cc(..., iv_to)
  subroutine a$D_box_sub_cc(box, iv_from, iv_to)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: iv_from, iv_to
#if $D == 2
    box%cc(:,:, iv_to) = box%cc(:,:, iv_to) - box%cc(:,:, iv_from)
#elif $D == 3
    box%cc(:,:,:, iv_to) = box%cc(:,:,:, iv_to) - box%cc(:,:,:, iv_from)
#endif
  end subroutine a$D_box_sub_cc

  !> Multipy cc(..., iv) with a
  subroutine a$D_box_times_cc(box, a, iv)
    type(box$D_t), intent(inout) :: box
    real(dp), intent(in)        :: a
    integer, intent(in)         :: iv
#if $D == 2
    box%cc(:,:, iv) = a * box%cc(:,:, iv)
#elif $D == 3
    box%cc(:,:,:, iv) = a * box%cc(:,:,:, iv)
#endif
  end subroutine a$D_box_times_cc

  !> Set cc(..., iv_b) = a * cc(..., iv_a) + b * cc(..., iv_b)
  subroutine a$D_box_lincomb_cc(box, a, iv_a, b, iv_b)
    type(box$D_t), intent(inout) :: box
    real(dp), intent(in)        :: a, b
    integer, intent(in)         :: iv_a, iv_b
#if $D == 2
    box%cc(:,:, iv_b) = a * box%cc(:,:, iv_a) + b * box%cc(:,:, iv_b)
#elif $D == 3
    box%cc(:,:,:, iv_b) = a * box%cc(:,:,:, iv_a) + b * box%cc(:,:,:, iv_b)
#endif
  end subroutine a$D_box_lincomb_cc

  !> Copy cc(..., iv_from) from box_in to cc(..., iv_to) on box_out
  subroutine a$D_box_copy_cc_to(box_from, iv_from, box_to, iv_to)
    type(box$D_t), intent(in)    :: box_from
    type(box$D_t), intent(inout) :: box_to
    integer, intent(in)         :: iv_from, iv_to
#if $D == 2
    box_to%cc(:,:, iv_to) = box_from%cc(:,:, iv_from)
#elif $D == 3
    box_to%cc(:,:,:, iv_to) = box_from%cc(:,:,:, iv_from)
#endif
  end subroutine a$D_box_copy_cc_to

  !> Copy cc(..., iv_from) to box%cc(..., iv_to)
  subroutine a$D_box_copy_cc(box, iv_from, iv_to)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: iv_from, iv_to
#if $D == 2
    box%cc(:,:, iv_to) = box%cc(:,:, iv_from)
#elif $D == 3
    box%cc(:,:,:, iv_to) = box%cc(:,:,:, iv_from)
#endif
  end subroutine a$D_box_copy_cc

  !> Copy cc(..., iv_from) to box%cc(..., iv_to) for all ids
  subroutine a$D_boxes_copy_cc(boxes, ids, iv_from, iv_to)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:), iv_from, iv_to
    integer                     :: i

    !$omp parallel do
    do i = 1, size(ids)
       call a$D_box_copy_cc(boxes(ids(i)), iv_from, iv_to)
    end do
    !$omp end parallel do
  end subroutine a$D_boxes_copy_cc

  !> Copy cc(..., iv_from) to box%cc(..., iv_to) for full tree
  subroutine a$D_tree_copy_cc(tree, iv_from, iv_to)
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)       :: iv_from, iv_to
    integer                   :: lvl

    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       call a$D_boxes_copy_cc(tree%boxes, tree%lvls(lvl)%ids, iv_from, iv_to)
    end do
  end subroutine a$D_tree_copy_cc

    !> A general scalar reduction method
  subroutine a$D_reduction(tree, box_func, reduction, init_val, out_val)
    type(a$D_t), intent(in) :: tree    !< Tree to do the reduction on
    real(dp), intent(in)   :: init_val !< Initial value for the reduction
    real(dp), intent(out)  :: out_val  !< Result of the reduction
    real(dp)               :: tmp, my_val
    integer                :: i, id, lvl

    interface
       real(dp) function box_func(box)
         import
         type(box$D_t), intent(in) :: box
       end function box_func

       real(dp) function reduction(a, b)
         import
         real(dp), intent(in) :: a, b
       end function reduction
    end interface

    if (.not. tree%ready) stop "Tree not ready"
    out_val = init_val
    my_val  = init_val

    !$omp parallel private(lvl, i, id, tmp) firstprivate(my_val)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          tmp = box_func(tree%boxes(id))
          my_val = reduction(tmp, my_val)
       end do
       !$omp end do
    end do

    !$omp critical
    out_val = reduction(my_val, out_val)
    !$omp end critical
    !$omp end parallel
  end subroutine a$D_reduction

  !> A general scalar reduction method, that returns the location of the
  !> minimum/maximum value found
  !> TODO: test
  subroutine a$D_reduction_loc(tree, box_subr, reduction, &
       init_val, out_val, out_loc)
    type(a$D_t), intent(in)      :: tree     !< Tree to do the reduction on
    real(dp), intent(in)        :: init_val !< Initial value for the reduction
    real(dp), intent(out)       :: out_val  !< Result of the reduction
    type(a$D_loc_t), intent(out) :: out_loc  !< Location
    real(dp)                    :: tmp, new_val, my_val
    integer                     :: i, id, lvl, tmp_ix($D)
    type(a$D_loc_t)             :: my_loc

    interface
       subroutine box_subr(box, val, ix)
         import
         type(box$D_t), intent(in) :: box
         real(dp), intent(out)    :: val
         integer, intent(out)     :: ix($D)
       end subroutine box_subr

       real(dp) function reduction(a, b)
         import
         real(dp), intent(in) :: a, b
       end function reduction
    end interface

    if (.not. tree%ready) stop "Tree not ready"
    out_val   = init_val
    my_val    = init_val
    my_loc%id = -1
    my_loc%ix = -1

    !$omp parallel private(lvl, i, id, tmp, tmp_ix, new_val) &
    !$omp firstprivate(my_val, my_loc)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call box_subr(tree%boxes(id), tmp, tmp_ix)
          new_val = reduction(tmp, my_val)
          if (abs(new_val - my_val) > 0) then
             my_loc%id = id
             my_loc%ix = tmp_ix
             my_val = tmp
          end if
       end do
       !$omp end do
    end do

    !$omp critical
    new_val = reduction(my_val, out_val)
    if (abs(new_val - out_val) > 0) then
       out_loc%id = my_loc%id
       out_loc%ix = my_loc%ix
       out_val = my_val
    end if
    !$omp end critical
    !$omp end parallel
  end subroutine a$D_reduction_loc

  !> Find maximum value of cc(..., iv). Only loop over leaves, and ghost cells
  !> are not used.
  subroutine a$D_tree_max_cc(tree, iv, cc_max)
    type(a$D_t), intent(in) :: tree
    integer, intent(in)    :: iv
    real(dp), intent(out)  :: cc_max
    real(dp)               :: tmp, my_max
    integer                :: i, id, lvl, nc

    if (.not. tree%ready) stop "Tree not ready"
    my_max = -huge(1.0_dp)

    !$omp parallel reduction(max: my_max) private(lvl, i, id, nc, tmp)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          nc = tree%boxes(id)%n_cell
#if $D == 2
          tmp = maxval(tree%boxes(id)%cc(1:nc, 1:nc, iv))
#elif $D == 3
          tmp = maxval(tree%boxes(id)%cc(1:nc, 1:nc, 1:nc, iv))
#endif
          if (tmp > my_max) my_max = tmp
       end do
       !$omp end do
    end do
    !$omp end parallel

    cc_max = my_max
  end subroutine a$D_tree_max_cc

  !> Find minimum value of cc(..., iv). Only loop over leaves, and ghost cells
  !> are not used.
  subroutine a$D_tree_min_cc(tree, iv, cc_min)
    type(a$D_t), intent(in) :: tree
    integer, intent(in)    :: iv
    real(dp), intent(out)  :: cc_min
    real(dp)               :: tmp, my_min
    integer                :: i, id, lvl, nc

    if (.not. tree%ready) stop "Tree not ready"
    my_min = huge(1.0_dp)

    !$omp parallel reduction(min: my_min) private(lvl, i, id, nc, tmp)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          nc = tree%boxes(id)%n_cell
#if $D == 2
          tmp = minval(tree%boxes(id)%cc(1:nc, 1:nc, iv))
#elif $D == 3
          tmp = minval(tree%boxes(id)%cc(1:nc, 1:nc, 1:nc, iv))
#endif
          if (tmp < my_min) my_min = tmp
       end do
       !$omp end do
    end do
    !$omp end parallel

    cc_min = my_min
  end subroutine a$D_tree_min_cc

  !> Find weighted sum of cc(..., iv). Only loop over leaves, and ghost cells
  !> are not used.
  subroutine a$D_tree_sum_cc(tree, iv, cc_sum)
    type(a$D_t), intent(in) :: tree
    integer, intent(in)    :: iv
    real(dp), intent(out)  :: cc_sum
    real(dp)               :: tmp, my_sum, fac
    integer                :: i, id, lvl, nc

    if (.not. tree%ready) stop "Tree not ready"
    my_sum = 0

    !$omp parallel reduction(+: my_sum) private(lvl, i, id, nc, tmp, fac)
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       fac = a$D_lvl_dr(tree, lvl)**$D

       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          nc = tree%boxes(id)%n_cell
#if $D == 2
          if (tree%coord_t == a5_cyl) then
             tmp = sum_2pr_box(tree%boxes(id), iv)
          else
             tmp = sum(tree%boxes(id)%cc(1:nc, 1:nc, iv))
          end if
#elif $D == 3
          tmp = sum(tree%boxes(id)%cc(1:nc, 1:nc, 1:nc, iv))
#endif
          my_sum = my_sum + fac * tmp
       end do
       !$omp end do
    end do
    !$omp end parallel

    cc_sum = my_sum

#if $D == 2
  contains

    ! Sum of 2 * pi * r * values
    pure function sum_2pr_box(box, iv) result(res)
      type(box2_t), intent(in) :: box
      integer, intent(in)      :: iv
      real(dp), parameter      :: twopi = 2 * acos(-1.0_dp)
      real(dp)                 :: res
      integer                  :: i, j, nc

      res = 0
      nc  = box%n_cell

      do j = 1, nc
         do i = 1, nc
            res = res + box%cc(i, j, iv) * a2_cyl_radius_cc(box, [i, j])
         end do
      end do
      res = res * twopi
    end function sum_2pr_box
#endif
  end subroutine a$D_tree_sum_cc

  !> Copy fx/fy/fz(..., iv_from) to fx/fy/fz(..., iv_to)
  subroutine a$D_box_copy_fc(box, iv_from, iv_to)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: iv_from, iv_to
#if $D == 2
    box%fx(:,:, iv_to) = box%fx(:,:, iv_from)
    box%fy(:,:, iv_to) = box%fy(:,:, iv_from)
#elif $D == 3
    box%fx(:,:,:, iv_to) = box%fx(:,:,:, iv_from)
    box%fy(:,:,:, iv_to) = box%fy(:,:,:, iv_from)
    box%fz(:,:,:, iv_to) = box%fz(:,:,:, iv_from)
#endif
  end subroutine a$D_box_copy_fc

  !> Copy fx/fy/fz(..., iv_from) to fx/fy/fz(..., iv_to) for all ids
  subroutine a$D_boxes_copy_fc(boxes, ids, iv_from, iv_to)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:), iv_from, iv_to
    integer                     :: i

    !$omp parallel do
    do i = 1, size(ids)
       call a$D_box_copy_fc(boxes(ids(i)), iv_from, iv_to)
    end do
    !$omp end parallel do
  end subroutine a$D_boxes_copy_fc

  !> Copy fx/fy/fz(..., iv_from) to fx/fy/fz(..., iv_to) for full tree
  subroutine a$D_tree_copy_fc(tree, iv_from, iv_to)
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)       :: iv_from, iv_to
    integer                   :: lvl

    if (.not. tree%ready) stop "Tree not ready"
    do lvl = lbound(tree%lvls, 1), tree%max_lvl
       call a$D_boxes_copy_fc(tree%boxes, tree%lvls(lvl)%ids, iv_from, iv_to)
    end do
  end subroutine a$D_tree_copy_fc

    !> Return n_cell at lvl. For all lvls >= 1, n_cell has the same value, but
  !> for lvls <= 0, n_cell changes.
  !> TODO: remove this in future
  pure function a$D_n_cell(tree, lvl) result(n_cell)
    type(a$D_t), intent(in) :: tree
    integer, intent(in)    :: lvl
    integer                :: n_cell !< Output: n_cell at lvl

    if (lvl >= 1) then
       n_cell = tree%n_cell
    else
       n_cell = tree%n_cell / (2**(1-lvl))
    end if
  end function a$D_n_cell

end module m_a$D_utils
