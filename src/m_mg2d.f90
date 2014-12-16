module m_mg2d
  use m_afivo_2d

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  type mg2_t
     private
     integer :: i_phi
     integer :: i_phi_old
     integer :: i_rhs
     integer :: i_res
     integer :: n_cycle_down
     integer :: n_cycle_up
     integer :: n_cycle_base
     procedure(a2_subr_gc), pointer, nopass    :: sides_bc
     procedure(a2_subr_gc), pointer, nopass    :: corners_bc
     procedure(mg2d_box_op), pointer, nopass   :: box_op
     procedure(mg2d_box_gsrb), pointer, nopass :: box_gsrb
  end type mg2_t

  abstract interface
     subroutine mg2d_box_op(box, i_in, i_out, mg)
       import
       type(box2_t), intent(inout) :: box
       integer, intent(in)         :: i_in, i_out
       type(mg2_t), intent(in)     :: mg
     end subroutine mg2d_box_op

     subroutine mg2d_box_gsrb(box, i_phi, i_rhs, redblack_cntr, mg)
       import
       type(box2_t), intent(inout) :: box
       integer, intent(in)         :: i_phi, i_rhs, redblack_cntr
       type(mg2_t), intent(in)     :: mg
     end subroutine mg2d_box_gsrb
  end interface

  ! Public types
  public :: mg2d_box_op
  public :: mg2_t

  ! Public methods
  public :: mg2d_set
  public :: mg2d_create_subtree
  public :: mg2d_restrict_trees
  public :: mg2d_fas_fmg
  public :: mg2d_fas_vcycle
  public :: mg2d_lpl_box
  public :: mg2d_lpl_cyl_box
  public :: mg2d_gsrb_lpl_box
  public :: mg2d_gsrb_lpl_cyl_box

contains

  subroutine mg2d_set(mg, i_phi, i_phi_old, i_rhs, i_res, &
       n_cycle_down, n_cycle_up, n_cycle_base, &
       sides_bc, corners_bc, my_operator, my_gsrb)
    type(mg2_t), intent(out) :: mg
    integer, intent(in)      :: i_phi, i_phi_old, i_rhs, i_res
    integer, intent(in)      :: n_cycle_down, n_cycle_up, n_cycle_base
    procedure(a2_subr_gc)    :: sides_bc, corners_bc
    procedure(mg2d_box_op)   :: my_operator
    procedure(mg2d_box_gsrb)   :: my_gsrb

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
    mg%box_gsrb     => my_gsrb
  end subroutine mg2d_set

  subroutine mg2d_create_subtree(tree)
    type(a2_t), intent(inout) :: tree

    integer :: id, c_id, n, lvl, n_lvls
    integer :: min_lvl, boxes_per_lvl, offset
    type(lvl_t), allocatable :: tmp_lvls(:)

    ! Determine number of lvls for subtree
    n = tree%n_cell
    n_lvls = 0
    do
       n = ishft(n, -1)
       if (btest(n, 0)) exit     ! Exit if n is odd
       n_lvls = n_lvls + 1
    end do

    if (n_lvls == 0) stop "mg2d_create_subtree: the subtree has no levels"

    min_lvl = 1 - n_lvls

    ! Will always have this many boxes per level
    boxes_per_lvl = size(tree%lvls(1)%ids)

    ! Allocate subtree levels
    allocate(tmp_lvls(min_lvl:ubound(tree%lvls, 1)))
    tmp_lvls(1:) = tree%lvls
    deallocate(tree%lvls)
    call move_alloc(tmp_lvls, tree%lvls)

    ! Create coarser levels which are copies of lvl 1
    do lvl = 0, min_lvl, -1
       allocate(tree%lvls(lvl)%ids(boxes_per_lvl))
       allocate(tree%lvls(lvl)%parents(boxes_per_lvl))
       allocate(tree%lvls(lvl)%leaves(0))

       call a2_get_free_ids(tree, tree%lvls(lvl)%ids)
       tree%lvls(lvl)%parents = tree%lvls(lvl)%ids
       offset = tree%lvls(lvl)%ids(1) - tree%lvls(lvl+1)%ids(1)

       do n = 1, boxes_per_lvl
          c_id                        = tree%lvls(lvl+1)%ids(n)
          id                          = tree%lvls(lvl)%ids(n)
          tree%boxes(id)%lvl          = lvl
          tree%boxes(id)%ix           = tree%boxes(c_id)%ix
          tree%boxes(id)%tag          = ibset(0, a5_bit_in_use)
          tree%boxes(id)%dr           = tree%boxes(c_id)%dr * 2
          tree%boxes(id)%r_min        = tree%boxes(c_id)%r_min
          tree%boxes(id)%n_cell       = tree%boxes(c_id)%n_cell / 2

          tree%boxes(id)%parent       = a5_no_box
          tree%boxes(id)%children(1)  = c_id
          tree%boxes(id)%children(2:) = a5_no_box

          ! Connectivity stays the same
          tree%boxes(id)%neighbors = tree%boxes(c_id)%neighbors
          where (tree%boxes(id)%neighbors > a5_no_box)
             tree%boxes(id)%neighbors = &
                  tree%boxes(id)%neighbors + offset
          end where

          call alloc_box(tree%boxes(id), tree%boxes(id)%n_cell, &
               tree%n_var_cell, tree%n_var_face)
       end do
    end do
  end subroutine mg2d_create_subtree

  ! Restrict cell centered variables ivs(:) on the trees
  subroutine mg2d_restrict_trees(tree, iv, mg)
    type(a2_t), intent(inout)       :: tree
    integer, intent(in)             :: iv
    type(mg2_t), intent(in)         :: mg
    integer                         :: lvl

    ! Restrict phi and rhs on tree
    do lvl = tree%n_lvls-1, lbound(tree%lvls, 1), -1
       call a2_restrict_to_boxes(tree%boxes, tree%lvls(lvl)%parents, iv)
       ! TODO: doesn't make sense to set gc for rhs
       call mg2d_fill_gc(tree%boxes, tree%lvls(lvl)%ids, iv, &
               mg%sides_bc, mg%corners_bc)
    end do
  end subroutine mg2d_restrict_trees

  ! Need valid ghost cells on input, has valid gc on output
  subroutine mg2d_fas_fmg(tree, mg)
    type(a2_t), intent(inout)       :: tree
    type(mg2_t), intent(in)         :: mg
    integer                         :: lvl, min_lvl

    min_lvl = lbound(tree%lvls, 1)

    do lvl = min_lvl, tree%n_lvls
       ! Store phi in phi_old
       call a2_boxes_copy_cc(tree%boxes, tree%lvls(lvl)%ids, &
            mg%i_phi, mg%i_phi_old)

       if (lvl > min_lvl) then
          ! Correct solution at this lvl using lvl-1 data
          ! phi = phi + prolong(phi_coarse - phi_old_coarse)
          call correct_children(tree%boxes, tree%lvls(lvl-1)%parents, mg)

          ! Update ghost cells
          call mg2d_fill_gc(tree%boxes, tree%lvls(lvl)%ids, mg%i_phi, &
               mg%sides_bc, mg%corners_bc)
       end if

       ! Perform V-cycle
       call mg2d_fas_vcycle(tree, mg, lvl)
    end do

  end subroutine mg2d_fas_fmg

  ! On entrance, need valid ghost cell data. On exit, leave valid ghost cell
  ! data
  subroutine mg2d_fas_vcycle(tree, mg, max_lvl)
    type(a2_t), intent(inout)          :: tree
    type(mg2_t), intent(in)            :: mg
    integer, intent(in)                :: max_lvl
    integer                            :: i, id, lvl, min_lvl

    min_lvl = lbound(tree%lvls, 1)

    do lvl = max_lvl,  min_lvl+1, -1
       ! Downwards relaxation
       call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_down)

       ! Calculate residual at current lvl
       call residual_boxes(tree%boxes, tree%lvls(lvl)%ids, mg)

       call a2_restrict_to_boxes(tree%boxes, tree%lvls(lvl-1)%parents, mg%i_phi)
       call a2_restrict_to_boxes(tree%boxes, tree%lvls(lvl-1)%parents, mg%i_res)

       call mg2d_fill_gc(tree%boxes, tree%lvls(lvl-1)%ids, mg%i_phi, &
            mg%sides_bc, mg%corners_bc)

       ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined, and
       ! store current coarse phi in phi_old
       !$omp do private(id)
       do i = 1, size(tree%lvls(lvl-1)%parents)
          id = tree%lvls(lvl-1)%parents(i)
          call mg%box_op(tree%boxes(id), mg%i_phi, mg%i_rhs, mg)
          call a2_box_add_cc(tree%boxes(id), mg%i_res, mg%i_rhs)
          call a2_box_copy_cc(tree%boxes(id), mg%i_phi, mg%i_phi_old)
       end do
       !$omp end do
    end do

    lvl = min_lvl
    call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_base)

    ! Do the upwards part of the v-cycle in the tree
    do lvl = min_lvl+1, max_lvl
       ! Correct solution at this lvl using lvl-1 data
       ! phi = phi + prolong(phi_coarse - phi_old_coarse)
       call correct_children(tree%boxes, tree%lvls(lvl-1)%parents, mg)

       ! Have to fill ghost cells again (todo: not everywhere?)
       call mg2d_fill_gc(tree%boxes, tree%lvls(lvl)%ids, mg%i_phi, &
            mg%sides_bc, mg%corners_bc)

       ! Upwards relaxation
       call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_up)
    end do
  end subroutine mg2d_fas_vcycle

  subroutine mg2d_fill_gc(boxes, ids, iv, sides_bc, corners_bc)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:), iv
    procedure(a2_subr_gc)       :: sides_bc, corners_bc
    integer                     :: i

    !$omp do
    do i = 1, size(ids)
       call a2_gc_box_sides(boxes, ids(i), iv, &
            a2_sides_extrap, sides_bc)
    end do
    !$omp end do
    !$omp do
    do i = 1, size(ids)
       call a2_gc_box_corners(boxes, ids(i), iv, &
            a2_corners_extrap, corners_bc)
    end do
    !$omp end do
  end subroutine mg2d_fill_gc

  ! Sets phi = phi + prolong(phi_coarse - phi_old_coarse)
  subroutine correct_children(boxes, ids, mg)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:)
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, nc, i_c, c_id, ix_offset(2)

    !$omp do private(nc, i_c, c_id, ix_offset)
    do i = 1, size(ids)

       nc = boxes(ids(i))%n_cell
       do i_c = 1, 4
          c_id = boxes(ids(i))%children(i_c)
          if (c_id == a5_no_box) cycle

          ! Offset of child w.r.t. parent
          ix_offset = a2_ch_dix(:, i_c) * ishft(nc, -1)

          call correct_child_box(boxes(ids(i)), boxes(c_id), &
               ix_offset, mg%i_phi, mg%i_phi_old)
       end do
    end do
    !$omp end do
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
    use omp_lib
    type(box2_t), intent(inout) :: boxes(:)
    type(mg2_t), intent(in)     :: mg
    integer, intent(in)         :: ids(:), n_cycle
    integer                     :: n, i

    do n = 1, 2 * n_cycle
       !$omp do
       do i = 1, size(ids)
          call mg%box_gsrb(boxes(ids(i)), mg%i_phi, mg%i_rhs, n, mg)
       end do
       !$omp end do

       !$omp do
       do i = 1, size(ids)
          call a2_gc_box_sides(boxes, ids(i), mg%i_phi, &
               a2_sides_extrap, mg%sides_bc)
       end do
       !$omp end do

       !$omp do
       do i = 1, size(ids)
          call a2_gc_box_corners(boxes, ids(i), mg%i_phi, &
               a2_corners_extrap, mg%corners_bc)
       end do
       !$omp end do
    end do

  end subroutine gsrb_boxes

  subroutine mg2d_gsrb_lpl_cyl_box(box, i_phi, i_rhs, redblack_cntr, mg)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: i_phi, i_rhs, redblack_cntr
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, j, nc, di(2), ioff
    real(dp)                    :: dxdy

    dxdy = box%dr**2
    nc   = box%n_cell
    ioff = (box%ix(1)-1) * nc

    ! The parity of redblack_cntr determines which cells we use
    di(1) = iand(redblack_cntr, 1)
    di(2) = ieor(di(1), 1)

    do j = 1, nc, 2
       do i = 1+di(1), nc, 2
          box%cc(i, j, i_phi) = 0.25_dp * ( &
               (i+ioff)/(i+ioff-0.5_dp) * box%cc(i+1, j, i_phi) + &
               (i+ioff-1)/(i+ioff-0.5_dp) * box%cc(i-1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi) - &
               dxdy * box%cc(i, j, i_rhs))
       end do
    end do

    do j = 2, nc, 2
       do i = 1+di(2), nc, 2
          box%cc(i, j, i_phi) = 0.25_dp * ( &
               (i+ioff)/(i+ioff-0.5_dp) * box%cc(i+1, j, i_phi) + &
               (i+ioff-1)/(i+ioff-0.5_dp) * box%cc(i-1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi) - &
               dxdy * box%cc(i, j, i_rhs))
       end do
    end do
  end subroutine mg2d_gsrb_lpl_cyl_box

  subroutine mg2d_gsrb_lpl_box(box, i_phi, i_rhs, redblack_cntr, mg)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: i_phi, i_rhs, redblack_cntr
    type(mg2_t), intent(in)     :: mg
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
  end subroutine mg2d_gsrb_lpl_box

  subroutine mg2d_lpl_box(box, i_in, i_out, mg)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: i_in, i_out
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, j, nc
    real(dp)                    :: inv_dr_sq

    nc = box%n_cell
    inv_dr_sq = 1 / box%dr**2

    do j = 1, nc
       do i = 1, nc
          box%cc(i, j, i_out) = (box%cc(i-1, j, i_in) + box%cc(i+1, j, i_in) + &
               box%cc(i, j-1, i_in) + box%cc(i, j+1, i_in) - &
               4 * box%cc(i, j, i_in)) * inv_dr_sq
       end do
    end do
  end subroutine mg2d_lpl_box

  subroutine mg2d_lpl_cyl_box(box, i_in, i_out, mg)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: i_in, i_out
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, j, nc, ioff
    real(dp)                    :: inv_dr_sq

    nc = box%n_cell
    inv_dr_sq = 1 / box%dr**2
    ioff = (box%ix(1)-1) * nc

    do j = 1, nc
       do i = 1, nc
          box%cc(i, j, i_out) = &
               ((i+ioff-1)/(i+ioff-0.5_dp) * box%cc(i-1, j, i_in) + &
               (i+ioff)/(i+ioff-0.5_dp) * box%cc(i+1, j, i_in) + &
               box%cc(i, j-1, i_in) + box%cc(i, j+1, i_in) - &
               4 * box%cc(i, j, i_in)) * inv_dr_sq
       end do
    end do
  end subroutine mg2d_lpl_cyl_box

  subroutine residual_boxes(boxes, ids, mg)
    type(box2_t), intent(inout) :: boxes(:)
    type(mg2_t), intent(in)     :: mg
    integer, intent(in)         :: ids(:)
    integer                     :: i

    !$omp do
    do i = 1, size(ids)
       call residual_box(boxes(ids(i)), mg)
    end do
    !$omp end do
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
