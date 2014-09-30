module m_mg2d
  use m_afivo

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  type mg2_subt_t
     integer                   :: min_lvl
     type(lvl2_t), allocatable :: lvls(:)
     type(box2_t), allocatable :: boxes(:)
  end type mg2_subt_t

  type mg2_t
     private
     integer  :: i_phi
     integer  :: i_phi_old
     integer  :: i_rhs
     integer  :: i_res
     integer  :: n_cycle_down
     integer  :: n_cycle_up
     integer :: n_cycle_base
     procedure(a2_subr_gc), pointer, nopass    :: sides_bc
     procedure(a2_subr_gc), pointer, nopass    :: corners_bc
     procedure(mg2d_box_op), pointer, nopass :: box_op
  end type mg2_t

  abstract interface
     subroutine mg2d_box_op(box, i_in, i_out, mg)
       import
       type(box2_t), intent(inout) :: box
       integer, intent(in)         :: i_in, i_out
       type(mg2_t), intent(in)     :: mg
     end subroutine mg2d_box_op
  end interface

  ! Public types
  public :: mg2d_box_op
  public :: mg2_t
  public :: mg2_subt_t

  ! Public methods
  public :: mg2d_set
  public :: mg2d_create_subtree
  public :: mg2d_restrict_trees
  public :: mg2d_fas_fmg
  public :: mg2d_fas_vcycle
  public :: mg2d_laplacian_box

contains

  subroutine mg2d_set(mg, i_phi, i_phi_old, i_rhs, i_res, &
       n_cycle_down, n_cycle_up, n_cycle_base, &
       sides_bc, corners_bc, my_operator)
    type(mg2_t), intent(out) :: mg
    integer, intent(in)      :: i_phi, i_phi_old, i_rhs, i_res
    integer, intent(in)      :: n_cycle_down, n_cycle_up, n_cycle_base
    procedure(a2_subr_gc)    :: sides_bc, corners_bc
    procedure(mg2d_box_op)   :: my_operator

    mg%i_phi        = i_phi
    mg%i_phi_old    = i_phi_old
    mg%i_rhs        = i_rhs
    mg%i_res        = i_res
    mg%n_cycle_down = n_cycle_down
    mg%n_cycle_up   = n_cycle_up
    mg%n_cycle_base = n_cycle_base
    mg%sides_bc     => sides_bc
    mg%corners_bc   => corners_bc
    mg%box_op       => my_operator
  end subroutine mg2d_set

  subroutine mg2d_create_subtree(tree, subt)
    type(a2_t), intent(in)             :: tree
    type(mg2_subt_t), intent(inout) :: subt

    integer :: i, id, n, lvl, n_lvls, boxes_per_lvl, offset

    ! Determine number of lvls for subtree
    n = tree%n_cell
    n_lvls = 0
    do
       n = ishft(n, -1)
       if (btest(n, 0)) exit     ! Exit if n is odd
       n_lvls = n_lvls + 1
    end do

    if (n_lvls == 0) stop "mg2d_create_subtree: the subtree has no levels"
    subt%min_lvl = 1 - n_lvls

    ! Will always have this many boxes per level
    boxes_per_lvl = size(tree%lvls(1)%ids)

    ! Allocate tree
    allocate(subt%lvls(-n_lvls+1:0))
    allocate(subt%boxes(n_lvls*boxes_per_lvl))

    do lvl = 0, -n_lvls+1, -1
       allocate(subt%lvls(lvl)%ids(boxes_per_lvl))

       offset             = (n_lvls+lvl-1) * boxes_per_lvl
       subt%lvls(lvl)%ids = tree%lvls(1)%ids + offset

       do n = 1, boxes_per_lvl
          id                     = tree%lvls(1)%ids(n)
          i                      = subt%lvls(lvl)%ids(n)
          subt%boxes(i)%lvl      = lvl
          subt%boxes(i)%ix       = tree%boxes(id)%ix
          subt%boxes(i)%tag      = ibset(0, a5_bit_in_use)
          subt%boxes(i)%dr       = tree%dr_base * 2**(1-lvl)
          subt%boxes(i)%r_min    = tree%boxes(id)%r_min
          subt%boxes(i)%n_cell   = tree%n_cell / 2**(1-lvl)

          ! Parent and children are not set for the subtree
          subt%boxes(i)%parent   = a5_no_box
          subt%boxes(i)%children = a5_no_box

          ! But neighbors should be set
          subt%boxes(i)%neighbors = tree%boxes(id)%neighbors
          where (subt%boxes(i)%neighbors > a5_no_box)
             subt%boxes(i)%neighbors = &
                  subt%boxes(i)%neighbors + offset
          end where

          call alloc_box(subt%boxes(i), subt%boxes(i)%n_cell, &
               tree%n_var_cell, tree%n_var_face)
       end do
    end do

  end subroutine mg2d_create_subtree

  ! Restrict phi and rhs
  subroutine mg2d_restrict_trees(tree, subt, mg, use_subtree)
    type(a2_t), intent(inout)       :: tree
    type(mg2_subt_t), intent(inout) :: subt
    logical, intent(in)             :: use_subtree
    type(mg2_t), intent(in)         :: mg
    integer                         :: lvl, i, id, jd

    ! Restrict phi and rhs on tree
    do lvl = tree%n_lvls-1, 1, -1
       call a2_restrict_to_boxes(tree%boxes, tree%lvls(lvl)%parents, &
            [mg%i_rhs, mg%i_phi])
    end do

    if (use_subtree) then
       ! Move over to subtree
       do i = 1, size(tree%lvls(1)%ids)
          id = tree%lvls(1)%ids(i)
          jd = subt%lvls(0)%ids(i)
          call a2_restrict_box(tree%boxes(id), &
               subt%boxes(jd), [0, 0], [mg%i_rhs, mg%i_phi])
       end do

       ! Do rest of subtree
       do lvl = 0, subt%min_lvl+1, -1
          do i = 1, size(subt%lvls(lvl)%ids)
             id = subt%lvls(lvl)%ids(i)
             jd = subt%lvls(lvl-1)%ids(i)
             call a2_restrict_box(subt%boxes(id), &
                  subt%boxes(jd), [0, 0], [mg%i_rhs, mg%i_phi])
          end do
       end do
    end if
  end subroutine mg2d_restrict_trees

  ! Need valid ghost cells on input, has valid gc on output
  subroutine mg2d_fas_fmg(tree, subt, mg, use_subtree)
    type(a2_t), intent(inout)          :: tree
    type(mg2_subt_t), intent(inout) :: subt
    type(mg2_t), intent(in)            :: mg
    logical, intent(in) :: use_subtree
    integer                            :: i, id, jd, lvl

    if (use_subtree) then
       do lvl = subt%min_lvl, 0
          ! Store phi in phi_old
          call a2_boxes_copy_cc(subt%boxes, subt%lvls(lvl)%ids, &
               mg%i_phi, mg%i_phi_old)

          if (lvl > subt%min_lvl) then
             ! Correct solution at this lvl using lvl-1 data
             ! phi = phi + prolong(phi_coarse - phi_old_coarse)
             do i = 1, size(subt%lvls(lvl)%ids)
                id = subt%lvls(lvl)%ids(i)
                jd = subt%lvls(lvl-1)%ids(i)
                call correct_child_box(subt%boxes(jd), subt%boxes(id), &
                     [0,0], mg%i_phi, mg%i_phi_old)
             end do

             ! Update ghost cells
             call mg2d_fill_gc(subt%boxes, subt%lvls(lvl)%ids, [mg%i_phi], &
                  mg%sides_bc, mg%corners_bc)
          endif

          ! Perform V-cycle on subtree
          call mg2d_fas_vcycle_subtree(subt, mg, lvl)
       end do
    end if

    do lvl = 1, tree%n_lvls
       ! Store phi in phi_old
       call a2_boxes_copy_cc(tree%boxes, tree%lvls(lvl)%ids, &
            mg%i_phi, mg%i_phi_old)

       if (lvl > 1) then
          ! Correct solution at this lvl using lvl-1 data
          ! phi = phi + prolong(phi_coarse - phi_old_coarse)
          do i = 1, size(tree%lvls(lvl-1)%parents)
             id = tree%lvls(lvl-1)%parents(i)
             call correct_children(tree%boxes, id, mg%i_phi, mg%i_phi_old)
          end do

          ! Update ghost cells
          call mg2d_fill_gc(tree%boxes, tree%lvls(lvl)%ids, [mg%i_phi], &
               mg%sides_bc, mg%corners_bc)

       else if (lvl == 1 .and. use_subtree) then
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             jd = subt%lvls(lvl-1)%ids(i)
             call correct_child_box(subt%boxes(jd), tree%boxes(id), &
                  [0,0], mg%i_phi, mg%i_phi_old)
          end do
       endif

       ! Perform V-cycle
       call mg2d_fas_vcycle(tree, subt, mg, use_subtree, lvl)
    end do

  end subroutine mg2d_fas_fmg

  ! On entrance, need valid ghost cell data. On exit, leave valid ghost cell
  ! data
  recursive subroutine mg2d_fas_vcycle(tree, subt, mg, use_subtree, lvl)
    type(a2_t), intent(inout)          :: tree
    type(mg2_subt_t), intent(inout) :: subt
    type(mg2_t), intent(in)            :: mg
    logical, intent(in)                :: use_subtree
    integer, intent(in)                :: lvl
    integer                            :: i, id, jd

    if (use_subtree .and. lvl < 1) then
       ! We use a subtree to get an approximation for lvl 1 efficiently
       call mg2d_fas_vcycle_subtree(subt, mg, lvl)

    else if (.not. use_subtree .and. lvl == 1) then
       ! Perform base level relaxation because we don't have a subtree
       call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_base)

    else ! Normal part of the V-cycle

       ! Downwards relaxation
       call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_down)

       ! Calculate residual at current lvl
       call residual_boxes(tree%boxes, tree%lvls(lvl)%ids, mg)

       ! Restrict phi and res
       if (use_subtree .and. lvl == 1) then
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             jd = subt%lvls(lvl-1)%ids(i)
             call a2_restrict_box(tree%boxes(id), subt%boxes(jd), &
                  [0,0], [mg%i_phi, mg%i_res])
          end do
       else
          call a2_restrict_to_boxes(tree%boxes, tree%lvls(lvl-1)%parents, &
               [mg%i_phi, mg%i_res])
       end if

       ! Have to update ghost cells for phi_c (todo: not everywhere?)
       if (use_subtree .and. lvl == 1) then
          call mg2d_fill_gc(subt%boxes, subt%lvls(lvl-1)%ids, [mg%i_phi], &
               mg%sides_bc, mg%corners_bc)
       else
          call mg2d_fill_gc(tree%boxes, tree%lvls(lvl-1)%ids, [mg%i_phi], &
               mg%sides_bc, mg%corners_bc)
       end if

       ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined
       if (use_subtree .and. lvl == 1) then
          do i = 1, size(subt%lvls(lvl-1)%ids)
             id = subt%lvls(lvl-1)%ids(i)
             call mg%box_op(subt%boxes(id), mg%i_phi, mg%i_rhs, mg)
             call a2_box_add_cc(subt%boxes(id), mg%i_res, mg%i_rhs)
          end do
       else
          do i = 1, size(tree%lvls(lvl-1)%parents)
             id = tree%lvls(lvl-1)%parents(i)
             call mg%box_op(tree%boxes(id), mg%i_phi, mg%i_rhs, mg)
             call a2_box_add_cc(tree%boxes(id), mg%i_res, mg%i_rhs)
          end do
       end if

       ! Store current coarse phi in phi_old
       if (use_subtree .and. lvl == 1) then
          call a2_boxes_copy_cc(subt%boxes, subt%lvls(lvl-1)%ids, &
               mg%i_phi, mg%i_phi_old)
       else
          call a2_boxes_copy_cc(tree%boxes, tree%lvls(lvl-1)%ids, &
               mg%i_phi, mg%i_phi_old)
       end if

       call mg2d_fas_vcycle(tree, subt, mg, use_subtree, lvl-1)

       ! Correct solution at this lvl using lvl-1 data
       ! phi = phi + prolong(phi_coarse - phi_old_coarse)
       if (use_subtree .and. lvl == 1) then
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             jd = subt%lvls(lvl-1)%ids(i)
             call correct_child_box(subt%boxes(jd), tree%boxes(id), &
                  [0,0], mg%i_phi, mg%i_phi_old)
          end do
       else
          do i = 1, size(tree%lvls(lvl-1)%parents)
             id = tree%lvls(lvl-1)%parents(i)
             call correct_children(tree%boxes, id, mg%i_phi, mg%i_phi_old)
          end do
       end if

       ! Have to fill ghost cells again (todo: not everywhere?)
       call mg2d_fill_gc(tree%boxes, tree%lvls(lvl)%ids, [mg%i_phi], &
            mg%sides_bc, mg%corners_bc)

       ! Upwards relaxation
       call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_up)

    end if

    if (lvl > 0) then
       ! Calculate updated residual at current lvl for output
       call residual_boxes(tree%boxes, tree%lvls(lvl)%ids, mg)
    end if
  end subroutine mg2d_fas_vcycle

  recursive subroutine mg2d_fas_vcycle_subtree(subt, mg, lvl)
    type(mg2_subt_t), intent(inout) :: subt
    type(mg2_t), intent(in)         :: mg
    integer, intent(in)             :: lvl
    integer                         :: i, id, jd

    if (lvl == subt%min_lvl) then
       ! Perform base level relaxation on the subtree
       call gsrb_boxes(subt%boxes, subt%lvls(lvl)%ids, mg, mg%n_cycle_base)
    else
       ! Downwards relaxation
       call gsrb_boxes(subt%boxes, subt%lvls(lvl)%ids, mg, mg%n_cycle_down)

       ! Calculate residual
       call residual_boxes(subt%boxes, subt%lvls(lvl)%ids, mg)

       ! Restrict phi and res
       do i = 1, size(subt%lvls(lvl)%ids)
          id = subt%lvls(lvl)%ids(i)
          jd = subt%lvls(lvl-1)%ids(i)
          call a2_restrict_box(subt%boxes(id), subt%boxes(jd), &
               [0,0], [mg%i_phi, mg%i_res])
       end do

       ! Have to update ghost cells for phi_c (todo: not everywhere?)
       call mg2d_fill_gc(subt%boxes, subt%lvls(lvl-1)%ids, [mg%i_phi], &
            mg%sides_bc, mg%corners_bc)

       ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined
       do i = 1, size(subt%lvls(lvl-1)%ids)
          id = subt%lvls(lvl-1)%ids(i)
          call mg%box_op(subt%boxes(id), mg%i_phi, mg%i_rhs, mg)
          call a2_box_add_cc(subt%boxes(id), mg%i_res, mg%i_rhs)
       end do

       ! Store current coarse phi in phi_old
       call a2_boxes_copy_cc(subt%boxes, subt%lvls(lvl-1)%ids, &
            mg%i_phi, mg%i_phi_old)

       call mg2d_fas_vcycle_subtree(subt, mg, lvl-1)

       ! Correct solution at this lvl using lvl-1 data
       ! phi = phi + prolong(phi_coarse - phi_old_coarse)
       do i = 1, size(subt%lvls(lvl)%ids)
          id = subt%lvls(lvl)%ids(i)
          jd = subt%lvls(lvl-1)%ids(i)
          call correct_child_box(subt%boxes(jd), subt%boxes(id), &
               [0,0], mg%i_phi, mg%i_phi_old)
       end do

       ! Have to fill ghost cells again (todo: not everywhere?)
       call mg2d_fill_gc(subt%boxes, subt%lvls(lvl)%ids, [mg%i_phi], &
            mg%sides_bc, mg%corners_bc)

       ! Upwards relaxation
       call gsrb_boxes(subt%boxes, subt%lvls(lvl)%ids, mg, mg%n_cycle_up)
    end if
  end subroutine mg2d_fas_vcycle_subtree

  subroutine mg2d_fill_gc(boxes, ids, ivs, sides_bc, corners_bc)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:), ivs(:)
    procedure(a2_subr_gc)       :: sides_bc, corners_bc
    integer                     :: i

    do i = 1, size(ids)
       call a2_gc_box_sides(boxes, ids(i), ivs, &
            a2_sides_extrap, sides_bc)
    end do
    do i = 1, size(ids)
       call a2_gc_box_corners(boxes, ids(i), ivs, &
            a2_corners_extrap, corners_bc)
    end do
  end subroutine mg2d_fill_gc

  ! Sets phi = phi + prolong(phi_coarse - phi_old_coarse)
  subroutine correct_children(boxes, id, i_phi, i_phi_old)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, i_phi, i_phi_old
    integer                     :: nc, i_c, c_id, ix_offset(2)

    nc = boxes(id)%n_cell
    do i_c = 1, 4
       c_id = boxes(id)%children(i_c)
       ! Offset of child w.r.t. parent
       ix_offset = a2_ch_dix(:, i_c) * ishft(nc, -1)

       call correct_child_box(boxes(id), boxes(c_id), &
            ix_offset, i_phi, i_phi_old)
    end do
  end subroutine correct_children

  subroutine correct_child_box(box_p, box_c, ix_offset, i_phi, i_phi_old)
    type(box2_t), intent(inout) :: box_c
    type(box2_t), intent(in)    :: box_p
    integer, intent(in)         :: i_phi, i_phi_old, ix_offset(2)
    real(dp), parameter         :: f1=1/16.0_dp, f3=3/16.0_dp, f9=9/16.0_dp
    integer                     :: nc, i, j, i_c1, i_c2, j_c1, j_c2

    nc = box_c%n_cell
    ! In these loops, we calculate the closest coarse index (_c1), and the
    ! one-but-closest (_c2). The fine cell lies in between.
    do j = 1, nc
       j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

          box_c%cc(i, j, i_phi) = box_c%cc(i, j, i_phi) &
               + f9 * (box_p%cc(i_c1, j_c1, i_phi) &
               - box_p%cc(i_c1, j_c1, i_phi_old)) &
               + f3 * (box_p%cc(i_c2, j_c1, i_phi) &
               - box_p%cc(i_c2, j_c1, i_phi_old) &
               + box_p%cc(i_c1, j_c2, i_phi) &
               - box_p%cc(i_c1, j_c2, i_phi_old)) &
               + f1 * (box_p%cc(i_c2, j_c2, i_phi) &
               - box_p%cc(i_c2, j_c2, i_phi_old))
       end do
    end do
  end subroutine correct_child_box

  subroutine gsrb_boxes(boxes, ids, mg, n_cycle)
    type(box2_t), intent(inout) :: boxes(:)
    type(mg2_t), intent(in)     :: mg
    integer, intent(in)         :: ids(:), n_cycle
    integer                     :: n, i

    do n = 1, 2 * n_cycle
       do i = 1, size(ids)
          call gsrb_box(boxes(ids(i)), mg%i_phi, mg%i_rhs, n)
       end do

       ! Communicate updated boundary cells
       call mg2d_fill_gc(boxes, ids, [mg%i_phi], mg%sides_bc, mg%corners_bc)
    end do
  end subroutine gsrb_boxes

  subroutine gsrb_box(box, i_phi, i_rhs, redblack_cntr)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: i_phi, i_rhs, redblack_cntr
    integer                     :: i, j, nc, di(2)
    real(dp)                    :: dxdy

    dxdy = box%dr**2
    nc   = box%n_cell

    ! The parity of redblack_cntr determines which cells we use
    di(1) = iand(redblack_cntr, 1)
    di(2) = ieor(di(1), 1)

    do j = 1, nc, 2
       do i = 1+di(1), nc, 2
          box%cc(i, j, i_phi) = 0.25_dp * ( &
               box%cc(i+1, j, i_phi) + box%cc(i-1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi) - &
               dxdy * box%cc(i, j, i_rhs))
       end do
    end do

    do j = 2, nc, 2
       do i = 1+di(2), nc, 2
          box%cc(i, j, i_phi) = 0.25_dp * ( &
               box%cc(i+1, j, i_phi) + box%cc(i-1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi) - &
               dxdy * box%cc(i, j, i_rhs))
       end do
    end do
  end subroutine gsrb_box

  subroutine mg2d_laplacian_box(box, i_in, i_out, mg)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: i_in, i_out
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, j, nc
    real(dp)                    :: inv_dr_sq(2)

    nc = box%n_cell
    inv_dr_sq = 1 / box%dr**2

    do j = 1, nc
       do i = 1, nc
          box%cc(i, j, i_out) = inv_dr_sq(1) * (box%cc(i-1, j, i_in) &
               - 2 * box%cc(i, j, i_in) + box%cc(i+1, j, i_in)) + &
               inv_dr_sq(2) * (box%cc(i, j-1, i_in) &
               - 2 * box%cc(i, j, i_in) + box%cc(i, j+1, i_in))
       end do
    end do
  end subroutine mg2d_laplacian_box

  subroutine residual_boxes(boxes, ids, mg)
    type(box2_t), intent(inout) :: boxes(:)
    type(mg2_t), intent(in)     :: mg
    integer, intent(in)         :: ids(:)
    integer                     :: i

    do i = 1, size(ids)
       call residual_box(boxes(ids(i)), mg)
    end do
  end subroutine residual_boxes

  subroutine residual_box(box, mg)
    type(box2_t), intent(inout) :: box
    type(mg2_t), intent(in)     :: mg
    integer                     :: nc

    call mg%box_op(box, mg%i_phi, mg%i_res, mg)
    nc = box%n_cell
    box%cc(1:nc, 1:nc, mg%i_res) = box%cc(1:nc, 1:nc, mg%i_rhs) &
         - box%cc(1:nc, 1:nc, mg%i_res)
  end subroutine residual_box

end module m_mg2d