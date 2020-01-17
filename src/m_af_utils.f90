#include "cpp_macros.h"
!> This module contains all kinds of different 'helper' routines for Afivo. If
!> the number of routines for a particular topic becomes large, they should
!> probably be put in a separate module.
module m_af_utils
  use m_af_types

  implicit none
  private

  ! Public subroutines
  public :: af_loop_box
  public :: af_loop_box_arg
  public :: af_loop_boxes
  public :: af_loop_boxes_arg
  public :: af_loop_tree
  public :: af_loop_tree_arg
  public :: af_tree_clear_cc
  public :: af_box_clear_cc
  public :: af_tree_clear_ghostcells
  public :: af_box_clear_ghostcells
  public :: af_box_add_cc
  public :: af_box_sub_cc
  public :: af_tree_times_cc
  public :: af_tree_apply
  public :: af_box_times_cc
  public :: af_box_lincomb_cc
  public :: af_box_copy_cc_to
  public :: af_box_copy_cc
  public :: af_boxes_copy_cc
  public :: af_tree_copy_cc
  public :: af_reduction
  public :: af_reduction_vec
  public :: af_reduction_loc
  public :: af_tree_max_cc
  public :: af_tree_maxabs_cc
  public :: af_tree_min_cc
  public :: af_tree_max_fc
  public :: af_tree_min_fc
  public :: af_tree_sum_cc
  public :: af_box_copy_fc
  public :: af_boxes_copy_fc
  public :: af_tree_copy_fc

  ! Public functions
  public :: af_get_id_at
  public :: af_get_loc
  public :: af_r_loc
  public :: af_r_inside
  public :: af_n_cell

contains

  !> Call procedure for each box in tree
  subroutine af_loop_box(tree, my_procedure, leaves_only)
    type(af_t), intent(inout)     :: tree
    procedure(af_subr)           :: my_procedure
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                       :: lvl, i, id

    leaves = .false.; if (present(leaves_only)) leaves = leaves_only
    if (.not. tree%ready) stop "af_loop_box: set_base has not been called"

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
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
  end subroutine af_loop_box

  !> Call procedure for each box in tree, with argument rarg
  subroutine af_loop_box_arg(tree, my_procedure, rarg, leaves_only)
    type(af_t), intent(inout)     :: tree
    procedure(af_subr_arg)        :: my_procedure
    real(dp), intent(in)          :: rarg(:)
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                       :: lvl, i, id

    if (.not. tree%ready) stop "Tree not ready"
    leaves = .false.; if (present(leaves_only)) leaves = leaves_only

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
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
  end subroutine af_loop_box_arg

  !> Call procedure for each id in tree, giving the list of boxes
  subroutine af_loop_boxes(tree, my_procedure, leaves_only)
    type(af_t), intent(inout)     :: tree
    procedure(af_subr_boxes)      :: my_procedure
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                       :: lvl, i, id

    if (.not. tree%ready) stop "Tree not ready"
    leaves = .false.; if (present(leaves_only)) leaves = leaves_only

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
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
  end subroutine af_loop_boxes

  !> Call procedure for each id in tree, giving the list of boxes
  subroutine af_loop_boxes_arg(tree, my_procedure, rarg, leaves_only)
    type(af_t), intent(inout)    :: tree
    procedure(af_subr_boxes_arg) :: my_procedure
    real(dp), intent(in)         :: rarg(:)
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                      :: lvl, i, id

    if (.not. tree%ready) stop "Tree not ready"
    leaves = .false.; if (present(leaves_only)) leaves = leaves_only

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
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
  end subroutine af_loop_boxes_arg

  !> Call procedure for each id in tree, passing the tree as first argument
  subroutine af_loop_tree(tree, my_procedure, leaves_only)
    type(af_t), intent(inout)     :: tree
    procedure(af_subr_tree)       :: my_procedure
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                       :: lvl, i

    if (.not. tree%ready) stop "Tree not ready"
    leaves = .false.; if (present(leaves_only)) leaves = leaves_only

    !$omp parallel private(lvl, i)
    do lvl = 1, tree%highest_lvl
       if (leaves) then
          !$omp do
          do i = 1, size(tree%lvls(lvl)%leaves)
             call my_procedure(tree, tree%lvls(lvl)%leaves(i))
          end do
          !$omp end do
       else
          !$omp do
          do i = 1, size(tree%lvls(lvl)%ids)
             call my_procedure(tree, tree%lvls(lvl)%ids(i))
          end do
          !$omp end do
       end if
    end do
    !$omp end parallel
  end subroutine af_loop_tree

  !> Call procedure for each id in tree, passing the tree as first argument
  subroutine af_loop_tree_arg(tree, my_procedure, rarg, leaves_only)
    type(af_t), intent(inout)     :: tree
    procedure(af_subr_tree_arg)   :: my_procedure
    real(dp), intent(in)          :: rarg(:)
    logical, intent(in), optional :: leaves_only
    logical                       :: leaves
    integer                       :: lvl, i

    if (.not. tree%ready) stop "Tree not ready"
    leaves = .false.; if (present(leaves_only)) leaves = leaves_only

    !$omp parallel private(lvl, i)
    do lvl = 1, tree%highest_lvl
       if (leaves) then
          !$omp do
          do i = 1, size(tree%lvls(lvl)%leaves)
             call my_procedure(tree, tree%lvls(lvl)%leaves(i), rarg)
          end do
          !$omp end do
       else
          !$omp do
          do i = 1, size(tree%lvls(lvl)%ids)
             call my_procedure(tree, tree%lvls(lvl)%leaves(i), rarg)
          end do
          !$omp end do
       end if
    end do
    !$omp end parallel
  end subroutine af_loop_tree_arg

  !> Returns whether r is inside or within a distance d from box
  pure function af_r_inside(box, r, d) result(inside)
    type(box_t), intent(in)       :: box
    real(dp), intent(in)           :: r(NDIM)
    real(dp), intent(in), optional :: d
    real(dp)                       :: r_max(NDIM)
    logical                        :: inside

    r_max = box%r_min + box%dr * box%n_cell
    if (present(d)) then
       inside = all(r+d >= box%r_min) .and. all(r-d <= r_max)
    else
       inside = all(r >= box%r_min) .and. all(r <= r_max)
    end if
  end function af_r_inside

  !> Get the id of the finest box containing rr. If highest_lvl is present, do
  !> not go to a finer level than highest_lvl. If there is no box containing rr,
  !> return a location of -1
  pure function af_get_id_at(tree, rr, highest_lvl, guess) result(id)
    type(af_t), intent(in)       :: tree        !< Full grid
    real(dp), intent(in)          :: rr(NDIM)      !< Coordinate
    integer, intent(in), optional :: highest_lvl !< Maximum level of box
    integer, intent(in), optional :: guess       !< Guess of box id, cannot be used with highest_lvl
    integer                       :: id          !< Id of finest box containing rr

    integer                       :: i, id_tmp, i_ch, lvl_max

    lvl_max = af_max_lvl
    if (present(highest_lvl)) lvl_max = highest_lvl

    id = -1

    if (present(guess)) then
       ! Check whether the guess is valid
       if (guess > 0 .and. guess < tree%highest_id) then
          if (tree%boxes(guess)%in_use .and. &
               tree%boxes(guess)%lvl <= lvl_max .and. &
               af_r_inside(tree%boxes(guess), rr)) then
             id = guess
          end if
       end if
    end if

    if (id == -1) then
       ! If there was no (valid) guess, find lvl 1 box that includes rr
       do i = 1, size(tree%lvls(1)%ids)
          id_tmp = tree%lvls(1)%ids(i)
          if (af_r_inside(tree%boxes(id_tmp), rr)) then
             id = id_tmp
             exit
          end if
       end do
    end if

    if (id > 0) then
       ! Jump into children for as long as possible
       do
          if (tree%boxes(id)%lvl >= lvl_max .or. &
               .not. af_has_children(tree%boxes(id))) exit
          i_ch = child_that_contains(tree%boxes(id), rr)
          id = tree%boxes(id)%children(i_ch)
       end do
    end if

  end function af_get_id_at

  !> Get the location of the finest cell containing rr. If highest_lvl is present,
  !> do not go to a finer level than highest_lvl. If there is no box containing rr,
  !> return a location of -1
  pure function af_get_loc(tree, rr, highest_lvl, guess) result(loc)
    type(af_t), intent(in)       :: tree        !< Full grid
    real(dp), intent(in)          :: rr(NDIM)      !< Coordinate
    integer, intent(in), optional :: highest_lvl !< Maximum level of box
    !> Guess of box id, cannot be used with highest_lvl
    integer, intent(in), optional :: guess
    type(af_loc_t)               :: loc         !< Location of cell

    loc%id = af_get_id_at(tree, rr, highest_lvl, guess)

    if (loc%id == -1) then
       loc%ix = -1
    else
       loc%ix = af_cc_ix(tree%boxes(loc%id), rr)

       ! Fix indices for points exactly on the boundaries of a box (which could
       ! get a ghost cell index)
       where (loc%ix < 1) loc%ix = 1
       where (loc%ix > tree%n_cell) loc%ix = tree%n_cell
    end if
  end function af_get_loc

  !> For a box with children that contains rr, find in which child rr lies
  pure function child_that_contains(box, rr) result(i_ch)
    type(box_t), intent(in) :: box    !< A box with children
    real(dp), intent(in)      :: rr(NDIM) !< Location inside the box
    integer                   :: i_ch   !< Index of child containing rr
    real(dp)                  :: center(NDIM)

    i_ch   = 1
    center = box%r_min + box%dr * ishft(box%n_cell, -1)

    if (rr(1) > center(1)) i_ch = i_ch + 1
    if (rr(2) > center(2)) i_ch = i_ch + 2
#if NDIM==3
    if (rr(3) > center(3)) i_ch = i_ch + 4
#endif
  end function child_that_contains

  !> Get the coordinate of the cell-center at loc
  pure function af_r_loc(tree, loc) result(r)
    type(af_t), intent(in)     :: tree !< Full grid
    type(af_loc_t), intent(in) :: loc !< Location object
    real(dp)                   :: r(NDIM) !< Coordinate at cell center
    r = tree%boxes(loc%id)%r_min + &
         (loc%ix-0.5_dp) * tree%boxes(loc%id)%dr
  end function af_r_loc

  subroutine af_tree_clear_cc(tree, iv)
    type(af_t), intent(inout) :: tree
    integer, intent(in)        :: iv !< Variable to clear
    integer                    :: lvl, i, id

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call af_box_clear_cc(tree%boxes(id), iv)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine af_tree_clear_cc

  !> Set cc(..., iv) = 0
  subroutine af_box_clear_cc(box, iv)
    type(box_t), intent(inout) :: box
    integer, intent(in)         :: iv !< Variable to clear
#if NDIM == 2
    box%cc(:,:, iv) = 0
#elif NDIM == 3
    box%cc(:,:,:, iv) = 0
#endif
  end subroutine af_box_clear_cc

  subroutine af_tree_clear_ghostcells(tree, iv)
    type(af_t), intent(inout) :: tree
    integer, intent(in)        :: iv !< Variable to clear
    integer                    :: lvl, i, id

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call af_box_clear_ghostcells(tree%boxes(id), iv)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine af_tree_clear_ghostcells

  !> Set cc(..., iv) = 0
  subroutine af_box_clear_ghostcells(box, iv)
    type(box_t), intent(inout) :: box
    integer, intent(in)          :: iv !< Variable to clear
    integer                      :: nc

    nc = box%n_cell

#if NDIM == 2
    box%cc(0, :, iv)       = 0
    box%cc(nc+1, :, iv)    = 0
    box%cc(:, 0, iv)       = 0
    box%cc(:, nc+1, iv)    = 0
#elif NDIM == 3
    box%cc(0, :, :, iv)    = 0
    box%cc(nc+1, :, :, iv) = 0
    box%cc(:, 0, :, iv)    = 0
    box%cc(:, nc+1, :, iv) = 0
    box%cc(:, :, 0, iv)    = 0
    box%cc(:, :, nc+1, iv) = 0
#endif
  end subroutine af_box_clear_ghostcells

  !> Add cc(..., iv_from) to box%cc(..., iv_to)
  subroutine af_box_add_cc(box, iv_from, iv_to)
    type(box_t), intent(inout) :: box
    integer, intent(in)         :: iv_from, iv_to
#if NDIM == 2
    box%cc(:,:, iv_to) = box%cc(:,:, iv_to) + box%cc(:,:, iv_from)
#elif NDIM == 3
    box%cc(:,:,:, iv_to) = box%cc(:,:,:, iv_to) + box%cc(:,:,:, iv_from)
#endif
  end subroutine af_box_add_cc

  !> Subtract cc(..., iv_from) from box%cc(..., iv_to)
  subroutine af_box_sub_cc(box, iv_from, iv_to)
    type(box_t), intent(inout) :: box
    integer, intent(in)         :: iv_from, iv_to
#if NDIM == 2
    box%cc(:,:, iv_to) = box%cc(:,:, iv_to) - box%cc(:,:, iv_from)
#elif NDIM == 3
    box%cc(:,:,:, iv_to) = box%cc(:,:,:, iv_to) - box%cc(:,:,:, iv_from)
#endif
  end subroutine af_box_sub_cc

  !> Perform cc(..., iv_a) = cc(..., iv_a) 'op' cc(..., iv_b), where 'op' can be
  !> '+', '*' or '/'
  subroutine af_tree_apply(tree, iv_a, iv_b, op, eps)
    type(af_t), intent(inout)     :: tree
    integer, intent(in)            :: iv_a, iv_b
    character(len=*), intent(in)   :: op
    real(dp), intent(in), optional :: eps !< Min value for division
    integer                        :: lvl, i, id
    real(dp)                       :: use_eps

    use_eps = sqrt(tiny(1.0_dp))
    if (present(eps)) use_eps = eps

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          select case (op)
          case ('+')
             tree%boxes(id)%cc(DTIMES(:), iv_a) = &
                  tree%boxes(id)%cc(DTIMES(:), iv_a) + &
                  tree%boxes(id)%cc(DTIMES(:), iv_b)
          case ('*')
             tree%boxes(id)%cc(DTIMES(:), iv_a) = &
                  tree%boxes(id)%cc(DTIMES(:), iv_a) * &
                  tree%boxes(id)%cc(DTIMES(:), iv_b)
          case ('/')
             tree%boxes(id)%cc(DTIMES(:), iv_a) = &
                  tree%boxes(id)%cc(DTIMES(:), iv_a) / &
                  max(tree%boxes(id)%cc(DTIMES(:), iv_b), eps)
          case default
             error stop "af_tree_apply: unknown operand"
          end select
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine af_tree_apply

  subroutine af_tree_times_cc(tree, ivs, facs)
    type(af_t), intent(inout) :: tree
    integer, intent(in)        :: ivs(:)
    real(dp), intent(in)       :: facs(:)
    integer                    :: lvl, i, id, n

    if (size(ivs) /= size(facs)) &
         error stop "af_times_cc: invalid array size"

    !$omp parallel private(lvl, i, id, n)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          do n = 1, size(ivs)
             call af_box_times_cc(tree%boxes(id), facs(n), ivs(n))
          end do
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine af_tree_times_cc

  !> Multipy cc(..., iv) with a
  subroutine af_box_times_cc(box, a, iv)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)        :: a
    integer, intent(in)         :: iv
#if NDIM == 2
    box%cc(:,:, iv) = a * box%cc(:,:, iv)
#elif NDIM == 3
    box%cc(:,:,:, iv) = a * box%cc(:,:,:, iv)
#endif
  end subroutine af_box_times_cc

  !> Set cc(..., iv_b) = a * cc(..., iv_a) + b * cc(..., iv_b)
  subroutine af_box_lincomb_cc(box, a, iv_a, b, iv_b)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)        :: a, b
    integer, intent(in)         :: iv_a, iv_b
#if NDIM == 2
    box%cc(:,:, iv_b) = a * box%cc(:,:, iv_a) + b * box%cc(:,:, iv_b)
#elif NDIM == 3
    box%cc(:,:,:, iv_b) = a * box%cc(:,:,:, iv_a) + b * box%cc(:,:,:, iv_b)
#endif
  end subroutine af_box_lincomb_cc

  !> Copy cc(..., iv_from) from box_in to cc(..., iv_to) on box_out
  subroutine af_box_copy_cc_to(box_from, iv_from, box_to, iv_to)
    type(box_t), intent(in)    :: box_from
    type(box_t), intent(inout) :: box_to
    integer, intent(in)         :: iv_from, iv_to
#if NDIM == 2
    box_to%cc(:,:, iv_to) = box_from%cc(:,:, iv_from)
#elif NDIM == 3
    box_to%cc(:,:,:, iv_to) = box_from%cc(:,:,:, iv_from)
#endif
  end subroutine af_box_copy_cc_to

  !> Copy cc(..., iv_from) to box%cc(..., iv_to)
  subroutine af_box_copy_cc(box, iv_from, iv_to)
    type(box_t), intent(inout) :: box
    integer, intent(in)         :: iv_from, iv_to
#if NDIM == 2
    box%cc(:,:, iv_to) = box%cc(:,:, iv_from)
#elif NDIM == 3
    box%cc(:,:,:, iv_to) = box%cc(:,:,:, iv_from)
#endif
  end subroutine af_box_copy_cc

  !> Copy cc(..., iv_from) to box%cc(..., iv_to) for all ids
  subroutine af_boxes_copy_cc(boxes, ids, iv_from, iv_to)
    type(box_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:), iv_from, iv_to
    integer                     :: i

    !$omp parallel do
    do i = 1, size(ids)
       call af_box_copy_cc(boxes(ids(i)), iv_from, iv_to)
    end do
    !$omp end parallel do
  end subroutine af_boxes_copy_cc

  !> Copy cc(..., iv_from) to box%cc(..., iv_to) for full tree
  subroutine af_tree_copy_cc(tree, iv_from, iv_to)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: iv_from, iv_to
    integer                   :: lvl

    do lvl = 1, tree%highest_lvl
       call af_boxes_copy_cc(tree%boxes, tree%lvls(lvl)%ids, iv_from, iv_to)
    end do
  end subroutine af_tree_copy_cc

  !> A general scalar reduction method
  subroutine af_reduction(tree, box_func, reduction, init_val, out_val)
    type(af_t), intent(in) :: tree    !< Tree to do the reduction on
    real(dp), intent(in)   :: init_val !< Initial value for the reduction
    real(dp), intent(out)  :: out_val  !< Result of the reduction
    real(dp)               :: tmp, my_val
    integer                :: i, id, lvl

    interface
       !> Function that returns a scalar
       real(dp) function box_func(box)
         import
         type(box_t), intent(in) :: box
       end function box_func

       !> Reduction method (e.g., min, max, sum)
       real(dp) function reduction(a, b)
         import
         real(dp), intent(in) :: a, b
       end function reduction
    end interface

    if (.not. tree%ready) stop "Tree not ready"
    out_val = init_val
    my_val  = init_val

    !$omp parallel private(lvl, i, id, tmp) firstprivate(my_val)
    do lvl = 1, tree%highest_lvl
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
  end subroutine af_reduction

  !> A general vector reduction method
  subroutine af_reduction_vec(tree, box_func, reduction, init_val, &
       out_val, n_vals)
    type(af_t), intent(in) :: tree             !< Tree to do the reduction on
    integer, intent(in)     :: n_vals           !< Size of vector
    real(dp), intent(in)    :: init_val(n_vals) !< Initial value for the reduction
    real(dp), intent(out)   :: out_val(n_vals)  !< Result of the reduction
    real(dp)                :: tmp(n_vals), my_val(n_vals)
    integer                 :: i, id, lvl

    interface
       !> Function that returns a vector
       function box_func(box, n_vals) result(vec)
         import
         type(box_t), intent(in) :: box
         integer, intent(in)       :: n_vals
         real(dp)                  :: vec(n_vals)
       end function box_func

       !> Reduction method (e.g., min, max, sum)
       function reduction(vec_1, vec_2, n_vals) result(vec)
         import
         integer, intent(in)  :: n_vals
         real(dp), intent(in) :: vec_1(n_vals), vec_2(n_vals)
         real(dp)             :: vec(n_vals)
       end function reduction
    end interface

    if (.not. tree%ready) stop "Tree not ready"
    out_val = init_val
    my_val  = init_val

    !$omp parallel private(lvl, i, id, tmp) firstprivate(my_val)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          tmp = box_func(tree%boxes(id), n_vals)
          my_val = reduction(tmp, my_val, n_vals)
       end do
       !$omp end do
    end do

    !$omp critical
    out_val = reduction(my_val, out_val, n_vals)
    !$omp end critical
    !$omp end parallel
  end subroutine af_reduction_vec

  !> A general scalar reduction method, that returns the location of the
  !> minimum/maximum value found
  subroutine af_reduction_loc(tree, iv, box_subr, reduction, &
       init_val, out_val, out_loc)
    type(af_t), intent(in)      :: tree     !< Tree to do the reduction on
    integer, intent(in)          :: iv       !< Variable to operate on (can be ignored)
    real(dp), intent(in)         :: init_val !< Initial value for the reduction
    real(dp), intent(out)        :: out_val  !< Result of the reduction
    type(af_loc_t), intent(out) :: out_loc  !< Location
    real(dp)                     :: tmp, new_val, my_val
    integer                      :: i, id, lvl, tmp_ix(NDIM)
    type(af_loc_t)              :: my_loc

    interface
       !> Subroutine that returns a scalar and a cell index
       subroutine box_subr(box, iv, val, ix)
         import
         type(box_t), intent(in) :: box
         integer, intent(in)       :: iv
         real(dp), intent(out)     :: val
         integer, intent(out)      :: ix(NDIM)
       end subroutine box_subr

       !> Reduction method (e.g., min, max, sum)
       real(dp) function reduction(a, b)
         import
         real(dp), intent(in)      :: a, b
       end function reduction
    end interface

    if (.not. tree%ready) stop "Tree not ready"
    out_val   = init_val
    my_val    = init_val
    my_loc%id = -1
    my_loc%ix = -1

    !$omp parallel private(lvl, i, id, tmp, tmp_ix, new_val) &
    !$omp firstprivate(my_val, my_loc)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          call box_subr(tree%boxes(id), iv, tmp, tmp_ix)
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
  end subroutine af_reduction_loc

  !> Find maximum value of cc(..., iv). By default, only loop over leaves, and
  !> ghost cells are not used.
  subroutine af_tree_max_cc(tree, iv, cc_max, loc)
    type(af_t), intent(in)      :: tree   !< Full grid
    integer, intent(in)          :: iv     !< Index of variable
    real(dp), intent(out)        :: cc_max !< Maximum value
    !> Location of maximum
    type(af_loc_t), intent(out), optional :: loc
    type(af_loc_t)                        :: tmp_loc

    call af_reduction_loc(tree, iv, box_max_cc, reduce_max, &
         -huge(1.0_dp)/10, cc_max, tmp_loc)
    if (present(loc)) loc = tmp_loc
  end subroutine af_tree_max_cc

  !> Find maximum value of abs(cc(..., iv)). By default, only loop over leaves,
  !> and ghost cells are not used.
  subroutine af_tree_maxabs_cc(tree, iv, cc_max, loc)
    type(af_t), intent(in)      :: tree   !< Full grid
    integer, intent(in)          :: iv     !< Index of variable
    real(dp), intent(out)        :: cc_max !< Maximum value
    !> Location of maximum
    type(af_loc_t), intent(out), optional :: loc
    type(af_loc_t)                        :: tmp_loc

    call af_reduction_loc(tree, iv, box_maxabs_cc, reduce_max, &
         -huge(1.0_dp)/10, cc_max, tmp_loc)
    if (present(loc)) loc = tmp_loc
  end subroutine af_tree_maxabs_cc

  !> Find minimum value of cc(..., iv). By default, only loop over leaves, and
  !> ghost cells are not used.
  subroutine af_tree_min_cc(tree, iv, cc_min, loc)
    type(af_t), intent(in)      :: tree   !< Full grid
    integer, intent(in)          :: iv     !< Index of variable
    real(dp), intent(out)        :: cc_min !< Minimum value
    !> Location of minimum
    type(af_loc_t), intent(out), optional :: loc
    type(af_loc_t)                        :: tmp_loc

    call af_reduction_loc(tree, iv, box_min_cc, reduce_min, &
         huge(1.0_dp)/10, cc_min, tmp_loc)

    if (present(loc)) loc = tmp_loc
  end subroutine af_tree_min_cc

  !> Find maximum value of fc(..., dim, iv). By default, only loop over leaves.
  subroutine af_tree_max_fc(tree, dim, iv, fc_max, loc)
    type(af_t), intent(in)      :: tree   !< Full grid
    integer, intent(in)          :: dim    !< Flux dimension
    integer, intent(in)          :: iv     !< Index of face variable
    real(dp), intent(out)        :: fc_max !< Maximum value
    !> Location of maximum
    type(af_loc_t), intent(out), optional :: loc
    type(af_loc_t)                        :: tmp_loc
    integer                                :: dim_iv

    ! Encode dim and iv in a single variable
    dim_iv = (dim-1) * tree%n_var_face + iv - 1

    call af_reduction_loc(tree, dim_iv, box_max_fc, reduce_max, &
         -huge(1.0_dp)/10, fc_max, tmp_loc)
    if (present(loc)) loc = tmp_loc
  end subroutine af_tree_max_fc

  !> Find maximum value of fc(..., dim, iv). By default, only loop over leaves.
  subroutine af_tree_min_fc(tree, dim, iv, fc_min, loc)
    type(af_t), intent(in)      :: tree   !< Full grid
    integer, intent(in)          :: dim    !< Flux dimension
    integer, intent(in)          :: iv     !< Index of face variable
    real(dp), intent(out)        :: fc_min !< Minimum value
    !> Location of minimum
    type(af_loc_t), intent(out), optional :: loc
    type(af_loc_t)                        :: tmp_loc
    integer                                :: dim_iv

    ! Encode dim and iv in a single variable
    dim_iv = (dim-1) * tree%n_var_face + iv - 1

    call af_reduction_loc(tree, dim_iv, box_min_fc, reduce_min, &
         huge(1.0_dp)/10, fc_min, tmp_loc)
    if (present(loc)) loc = tmp_loc
  end subroutine af_tree_min_fc

  subroutine box_max_cc(box, iv, val, ix)
    type(box_t), intent(in) :: box
    integer, intent(in)       :: iv
    real(dp), intent(out)     :: val
    integer, intent(out)      :: ix(NDIM)
    integer                   :: nc

    nc = box%n_cell
#if NDIM == 2
    ix = maxloc(box%cc(1:nc, 1:nc, iv))
    val = box%cc(ix(1), ix(2), iv)
#elif NDIM == 3
    ix = maxloc(box%cc(1:nc, 1:nc, 1:nc, iv))
    val = box%cc(ix(1), ix(2), ix(3), iv)
#endif
  end subroutine box_max_cc

  subroutine box_maxabs_cc(box, iv, val, ix)
    type(box_t), intent(in) :: box
    integer, intent(in)       :: iv
    real(dp), intent(out)     :: val
    integer, intent(out)      :: ix(NDIM)
    integer                   :: nc

    nc = box%n_cell
#if NDIM == 2
    ix = maxloc(abs(box%cc(1:nc, 1:nc, iv)))
    val = abs(box%cc(ix(1), ix(2), iv))
#elif NDIM == 3
    ix = maxloc(abs(box%cc(1:nc, 1:nc, 1:nc, iv)))
    val = abs(box%cc(ix(1), ix(2), ix(3), iv))
#endif
  end subroutine box_maxabs_cc

  subroutine box_min_cc(box, iv, val, ix)
    type(box_t), intent(in) :: box
    integer, intent(in)       :: iv
    real(dp), intent(out)     :: val
    integer, intent(out)      :: ix(NDIM)
    integer                   :: nc

    nc = box%n_cell
#if NDIM == 2
    ix = minloc(box%cc(1:nc, 1:nc, iv))
    val = box%cc(ix(1), ix(2), iv)
#elif NDIM == 3
    ix = minloc(box%cc(1:nc, 1:nc, 1:nc, iv))
    val = box%cc(ix(1), ix(2), ix(3), iv)
#endif
  end subroutine box_min_cc

  subroutine box_max_fc(box, dim_iv, val, ix)
    type(box_t), intent(in) :: box
    integer, intent(in)       :: dim_iv
    real(dp), intent(out)     :: val
    integer, intent(out)      :: ix(NDIM)
    integer                   :: dim, iv, n_fc, nc, dix(NDIM)

    nc     = box%n_cell
    dix(:) = 0

#if NDIM == 2
    n_fc = size(box%fc, 4)

    ! Decode dim and iv
    dim      = dim_iv / n_fc + 1
    iv       = dim_iv - (dim-1) * n_fc + 1
    ! Also include fluxes at 'upper' boundary
    dix(dim) = 1
    ix       = maxloc(box%fc(1:nc+dix(1), 1:nc+dix(2), dim, iv))
    val      = box%fc(ix(1), ix(2), dim, iv)
#elif NDIM == 3
    n_fc     = size(box%fc, 5)
    ! Decode dim and iv
    dim      = dim_iv / n_fc + 1
    dix(dim) = 1
    iv       = dim_iv - (dim-1) * n_fc + 1
    ix       = maxloc(box%fc(1:nc+dix(1), 1:nc+dix(2), 1:nc+dix(3), dim, iv))
    val      = box%fc(ix(1), ix(2), ix(3), dim, iv)
#endif
  end subroutine box_max_fc

  subroutine box_min_fc(box, dim_iv, val, ix)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: dim_iv
    real(dp), intent(out)   :: val
    integer, intent(out)    :: ix(NDIM)
    integer                 :: dim, iv, n_fc, nc, dix(NDIM)

    nc     = box%n_cell
    dix(:) = 0

#if NDIM == 2
    n_fc = size(box%fc, 4)

    ! Decode dim and iv
    dim      = dim_iv / n_fc + 1
    iv       = dim_iv - (dim-1) * n_fc + 1
    ! Also include fluxes at 'upper' boundary
    dix(dim) = 1
    ix       = minloc(box%fc(1:nc+dix(1), 1:nc+dix(2), dim, iv))
    val      = box%fc(ix(1), ix(2), dim, iv)
#elif NDIM == 3
    n_fc     = size(box%fc, 5)
    ! Decode dim and iv
    dim      = dim_iv / n_fc + 1
    dix(dim) = 1
    iv       = dim_iv - (dim-1) * n_fc + 1
    ix       = minloc(box%fc(1:nc+dix(1), 1:nc+dix(2), 1:nc+dix(3), dim, iv))
    val      = box%fc(ix(1), ix(2), ix(3), dim, iv)
#endif
  end subroutine box_min_fc

  real(dp) function reduce_max(a, b)
    real(dp), intent(in) :: a, b
    reduce_max = max(a, b)
  end function reduce_max

  real(dp) function reduce_min(a, b)
    real(dp), intent(in) :: a, b
    reduce_min = min(a, b)
  end function reduce_min

  !> Find weighted sum of cc(..., iv). Only loop over leaves, and ghost cells
  !> are not used.
  subroutine af_tree_sum_cc(tree, iv, cc_sum)
    type(af_t), intent(in) :: tree !< Full grid
    integer, intent(in)    :: iv !< Index of variable
    real(dp), intent(out)  :: cc_sum !< Volume-integrated sum of variable
    real(dp)               :: tmp, my_sum, fac
    integer                :: i, id, lvl, nc

    if (.not. tree%ready) stop "Tree not ready"
    my_sum = 0

    !$omp parallel reduction(+: my_sum) private(lvl, i, id, nc, tmp, fac)
    do lvl = 1, tree%highest_lvl
       fac = product(af_lvl_dr(tree, lvl))

       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          nc = tree%boxes(id)%n_cell
#if NDIM == 2
          if (tree%coord_t == af_cyl) then
             tmp = sum_2pr_box(tree%boxes(id), iv)
          else
             tmp = sum(tree%boxes(id)%cc(1:nc, 1:nc, iv))
          end if
#elif NDIM == 3
          tmp = sum(tree%boxes(id)%cc(1:nc, 1:nc, 1:nc, iv))
#endif
          my_sum = my_sum + fac * tmp
       end do
       !$omp end do
    end do
    !$omp end parallel

    cc_sum = my_sum

#if NDIM == 2
  contains

    ! Sum of 2 * pi * r * values
    pure function sum_2pr_box(box, iv) result(res)
      type(box_t), intent(in) :: box
      integer, intent(in)      :: iv
      real(dp), parameter      :: twopi = 2 * acos(-1.0_dp)
      real(dp)                 :: res
      integer                  :: i, j, nc

      res = 0
      nc  = box%n_cell

      do j = 1, nc
         do i = 1, nc
            res = res + box%cc(i, j, iv) * af_cyl_radius_cc(box, i)
         end do
      end do
      res = res * twopi
    end function sum_2pr_box
#endif
  end subroutine af_tree_sum_cc

  !> Copy fx/fy/fz(..., iv_from) to fx/fy/fz(..., iv_to)
  subroutine af_box_copy_fc(box, iv_from, iv_to)
    type(box_t), intent(inout) :: box !< Operate on this box
    integer, intent(in)         :: iv_from !< From this variable
    integer, intent(in)         :: iv_to !< To this variable
#if NDIM == 2
    box%fc(:,:,:, iv_to) = box%fc(:,:,:, iv_from)
#elif NDIM == 3
    box%fc(:,:,:,:, iv_to) = box%fc(:,:,:,:, iv_from)
#endif
  end subroutine af_box_copy_fc

  !> Copy fx/fy/fz(..., iv_from) to fx/fy/fz(..., iv_to) for all ids
  subroutine af_boxes_copy_fc(boxes, ids, iv_from, iv_to)
    type(box_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: ids(:) !< Operate on these boxes
    integer, intent(in)         :: iv_from !< From this variable
    integer, intent(in)         :: iv_to !< To this variable
    integer                     :: i

    !$omp parallel do
    do i = 1, size(ids)
       call af_box_copy_fc(boxes(ids(i)), iv_from, iv_to)
    end do
    !$omp end parallel do
  end subroutine af_boxes_copy_fc

  !> Copy fx/fy/fz(..., iv_from) to fx/fy/fz(..., iv_to) for full tree
  subroutine af_tree_copy_fc(tree, iv_from, iv_to)
    type(af_t), intent(inout) :: tree !< Full grid
    integer, intent(in)       :: iv_from !< From this variable
    integer, intent(in)       :: iv_to !< To this variable
    integer                   :: lvl

    if (.not. tree%ready) stop "Tree not ready"
    do lvl = 1, tree%highest_lvl
       call af_boxes_copy_fc(tree%boxes, tree%lvls(lvl)%ids, iv_from, iv_to)
    end do
  end subroutine af_tree_copy_fc

  !> Return n_cell at lvl. For all lvls >= 1, n_cell has the same value, but
  !> for lvls <= 0, n_cell changes.
  !> @todo remove this in future
  pure function af_n_cell(tree, lvl) result(n_cell)
    type(af_t), intent(in) :: tree !< Full grid
    integer, intent(in)    :: lvl !< Refinement level
    integer                :: n_cell !< Output: n_cell at lvl

    if (lvl >= 1) then
       n_cell = tree%n_cell
    else
       n_cell = tree%n_cell / (2**(1-lvl))
    end if
  end function af_n_cell

end module m_af_utils
