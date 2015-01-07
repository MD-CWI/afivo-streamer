! Multigrid code for Xd simulations. The following replacements take place on
! this code:
! 1. $ D -> 2 or 3 (dimension of code)
! 2. preprocess file with cpp
! 3. cat -s (merge multiple blank lines)
module m_mg_$Dd
  use m_afivo_$Dd

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  type mg$D_t
     private
     integer :: i_phi
     integer :: i_phi_old
     integer :: i_rhs
     integer :: i_res
     integer :: n_cycle_down
     integer :: n_cycle_up
     integer :: n_cycle_base
     procedure(a$D_subr_gc), pointer, nopass    :: sides_bc
     procedure(a$D_subr_gc), pointer, nopass    :: corners_bc
     procedure(mg$Dd_box_op), pointer, nopass   :: box_op
     procedure(mg$Dd_box_gsrb), pointer, nopass :: box_gsrb
  end type mg$D_t

  abstract interface
     subroutine mg$Dd_box_op(box, i_in, i_out)
       import
       type(box$D_t), intent(inout) :: box
       integer, intent(in)         :: i_in, i_out
     end subroutine mg$Dd_box_op

     subroutine mg$Dd_box_gsrb(box, i_phi, i_rhs, redblack_cntr)
       import
       type(box$D_t), intent(inout) :: box
       integer, intent(in)         :: i_phi, i_rhs, redblack_cntr
     end subroutine mg$Dd_box_gsrb
  end interface

  ! Public types
  public :: mg$Dd_box_op
  public :: mg$D_t

  ! Public methods
  public :: mg$Dd_set
  public :: mg$Dd_fas_fmg
  public :: mg$Dd_fas_vcycle
  public :: mg$Dd_lpl_box
  public :: mg$Dd_gsrb_lpl_box

contains

  subroutine mg$Dd_set(mg, i_phi, i_phi_old, i_rhs, i_res, &
       n_cycle_down, n_cycle_up, n_cycle_base, &
       sides_bc, corners_bc, my_operator, my_gsrb)
    type(mg$D_t), intent(out) :: mg
    integer, intent(in)      :: i_phi, i_phi_old, i_rhs, i_res
    integer, intent(in)      :: n_cycle_down, n_cycle_up, n_cycle_base
    procedure(a$D_subr_gc)    :: sides_bc, corners_bc
    procedure(mg$Dd_box_op)   :: my_operator
    procedure(mg$Dd_box_gsrb)   :: my_gsrb

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
  end subroutine mg$Dd_set

  ! Need valid ghost cells on input, has valid gc on output
  subroutine mg$Dd_fas_fmg(tree, mg)
    type(a$D_t), intent(inout)       :: tree
    type(mg$D_t), intent(in)         :: mg
    integer                         :: lvl, min_lvl

    min_lvl = lbound(tree%lvls, 1)

    do lvl = min_lvl, tree%n_lvls
       ! Store phi in phi_old
       call a$D_boxes_copy_cc(tree%boxes, tree%lvls(lvl)%ids, &
            mg%i_phi, mg%i_phi_old)

       if (lvl > min_lvl) then
          ! Correct solution at this lvl using lvl-1 data
          ! phi = phi + prolong(phi_coarse - phi_old_coarse)
          call correct_children(tree%boxes, tree%lvls(lvl-1)%parents, mg)

          ! Update ghost cells
          call mg$Dd_fill_gc(tree%boxes, tree%lvls(lvl)%ids, mg%i_phi, &
               mg%sides_bc, mg%corners_bc)
       end if

       ! Perform V-cycle
       call mg$Dd_fas_vcycle(tree, mg, lvl)
    end do

  end subroutine mg$Dd_fas_fmg

  ! On entrance, need valid ghost cell data. On exit, leave valid ghost cell
  ! data
  subroutine mg$Dd_fas_vcycle(tree, mg, max_lvl)
    type(a$D_t), intent(inout)          :: tree
    type(mg$D_t), intent(in)            :: mg
    integer, intent(in)                :: max_lvl
    integer                            :: i, id, lvl, min_lvl

    min_lvl = lbound(tree%lvls, 1)

    do lvl = max_lvl,  min_lvl+1, -1
       ! Downwards relaxation
       call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_down)

       ! Calculate residual at current lvl
       call residual_boxes(tree%boxes, tree%lvls(lvl)%ids, mg)

       call a$D_restrict_to_boxes(tree%boxes, tree%lvls(lvl-1)%parents, mg%i_phi)
       call a$D_restrict_to_boxes(tree%boxes, tree%lvls(lvl-1)%parents, mg%i_res)

       call mg$Dd_fill_gc(tree%boxes, tree%lvls(lvl-1)%ids, mg%i_phi, &
            mg%sides_bc, mg%corners_bc)

       ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined, and
       ! store current coarse phi in phi_old
       !$omp do private(id)
       do i = 1, size(tree%lvls(lvl-1)%parents)
          id = tree%lvls(lvl-1)%parents(i)
          call mg%box_op(tree%boxes(id), mg%i_phi, mg%i_rhs)
          call a$D_box_add_cc(tree%boxes(id), mg%i_res, mg%i_rhs)
          call a$D_box_copy_cc(tree%boxes(id), mg%i_phi, mg%i_phi_old)
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
       call mg$Dd_fill_gc(tree%boxes, tree%lvls(lvl)%ids, mg%i_phi, &
            mg%sides_bc, mg%corners_bc)

       ! Upwards relaxation
       call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_up)
    end do
  end subroutine mg$Dd_fas_vcycle

  subroutine mg$Dd_fill_gc(boxes, ids, iv, sides_bc, corners_bc)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:), iv
    procedure(a$D_subr_gc)       :: sides_bc, corners_bc
    integer                     :: i

    !$omp do
    do i = 1, size(ids)
       call a$D_gc_box_sides(boxes, ids(i), iv, &
            a$D_sides_extrap, sides_bc)
    end do
    !$omp end do
    !$omp do
    do i = 1, size(ids)
       call a$D_gc_box_corners(boxes, ids(i), iv, &
            a$D_corners_extrap, corners_bc)
    end do
    !$omp end do
  end subroutine mg$Dd_fill_gc

  ! Sets phi = phi + prolong(phi_coarse - phi_old_coarse)
  subroutine correct_children(boxes, ids, mg)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:)
    type(mg$D_t), intent(in)     :: mg
    integer                     :: i, nc, i_c, c_id, ix_offset($D)

    !$omp do private(nc, i_c, c_id, ix_offset)
    do i = 1, size(ids)
       nc = boxes(ids(i))%n_cell

       do i_c = 1, a$D_num_children
          c_id = boxes(ids(i))%children(i_c)
          if (c_id == a5_no_box) cycle

          ! Offset of child w.r.t. parent
          ix_offset = a$D_ch_dix(:, i_c) * ishft(nc, -1)

          call correct_child_box(boxes(ids(i)), boxes(c_id), &
               ix_offset, mg%i_phi, mg%i_phi_old)
       end do
    end do
    !$omp end do
  end subroutine correct_children

  subroutine correct_child_box(box_p, box_c, ix_offset, i_phi, i_phi_old)
    type(box$D_t), intent(inout) :: box_c
    type(box$D_t), intent(in)    :: box_p
    integer, intent(in)         :: i_phi, i_phi_old, ix_offset($D)
    integer                     :: nc, i, j, i_c1, i_c2, j_c1, j_c2
#if $D == 2
    real(dp), parameter         :: f1=1/16.0_dp, f3=3/16.0_dp, f9=9/16.0_dp
#elif $D == 3
    integer                     :: k, k_c1, k_c2
    real(dp), parameter         :: f1=1/64.0_dp, f3=3/64.0_dp
    real(dp), parameter         :: f9=9/64.0_dp, f27=27/64.0_dp
#endif

    nc = box_c%n_cell
    ! In these loops, we calculate the closest coarse index (_c1), and the
    ! one-but-closest (_c2). The fine cell lies in between.
#if $D == 2
    do j = 1, nc
       j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

          box_c%cc(i, j, i_phi) = box_c%cc(i, j, i_phi) &
               + f9 * (box_p%cc(i_c1, j_c1, i_phi) &
               -      box_p%cc(i_c1, j_c1, i_phi_old)) &
               + f3 * (box_p%cc(i_c2, j_c1, i_phi) &
               -       box_p%cc(i_c2, j_c1, i_phi_old) &
               +       box_p%cc(i_c1, j_c2, i_phi) &
               -       box_p%cc(i_c1, j_c2, i_phi_old)) &
               + f1 * (box_p%cc(i_c2, j_c2, i_phi) &
               -       box_p%cc(i_c2, j_c2, i_phi_old))
       end do
    end do
#elif $D == 3
    do k = 1, nc
       k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
       k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)          ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)          ! even: +1, odd: -1

             box_c%cc(i, j, k, i_phi) = box_c%cc(i, j, k, i_phi) &
                  + f27 * (box_p%cc(i_c1, j_c1, k_c1, i_phi) &
                  -        box_p%cc(i_c1, j_c1, k_c1, i_phi_old)) &
                  + f9 * (box_p%cc(i_c2, j_c1, k_c1, i_phi) &
                  -       box_p%cc(i_c2, j_c1, k_c1, i_phi_old)) &
                  + f9 * (box_p%cc(i_c1, j_c2, k_c1, i_phi) &
                  -       box_p%cc(i_c1, j_c2, k_c1, i_phi_old)) &
                  + f9 * (box_p%cc(i_c1, j_c1, k_c2, i_phi) &
                  -       box_p%cc(i_c1, j_c1, k_c2, i_phi_old)) &
                  + f3 * (box_p%cc(i_c1, j_c2, k_c2, i_phi) &
                  -       box_p%cc(i_c1, j_c2, k_c2, i_phi_old)) &
                  + f3 * (box_p%cc(i_c2, j_c1, k_c2, i_phi) &
                  -       box_p%cc(i_c2, j_c1, k_c2, i_phi_old)) &
                  + f3 * (box_p%cc(i_c2, j_c2, k_c1, i_phi) &
                  -       box_p%cc(i_c2, j_c2, k_c1, i_phi_old)) &
                  + f1 * (box_p%cc(i_c2, j_c2, k_c2, i_phi) &
                  -       box_p%cc(i_c2, j_c2, k_c2, i_phi_old))
          end do
       end do
    end do
#endif
  end subroutine correct_child_box

  subroutine gsrb_boxes(boxes, ids, mg, n_cycle)
    type(box$D_t), intent(inout) :: boxes(:)
    type(mg$D_t), intent(in)     :: mg
    integer, intent(in)         :: ids(:), n_cycle
    integer                     :: n, i

    do n = 1, 2 * n_cycle
       !$omp do
       do i = 1, size(ids)
          call mg%box_gsrb(boxes(ids(i)), mg%i_phi, mg%i_rhs, n)
       end do
       !$omp end do

       !$omp do
       do i = 1, size(ids)
          call a$D_gc_box_sides(boxes, ids(i), mg%i_phi, &
               a$D_sides_extrap, mg%sides_bc)
       end do
       !$omp end do

       !$omp do
       do i = 1, size(ids)
          call a$D_gc_box_corners(boxes, ids(i), mg%i_phi, &
               a$D_corners_extrap, mg%corners_bc)
       end do
       !$omp end do
    end do
  end subroutine gsrb_boxes

  subroutine mg$Dd_gsrb_lpl_box(box, i_phi, i_rhs, redblack_cntr)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: i_phi, i_rhs, redblack_cntr
    integer                     :: i, i0, j, nc
    real(dp)                    :: dx2
#if $D == 3
    integer                     :: k
    real(dp), parameter         :: sixth = 1/6.0_dp
#endif

    dx2 = box%dr**2
    nc  = box%n_cell

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
#if $D == 2
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          box%cc(i, j, i_phi) = 0.25_dp * ( &
               box%cc(i+1, j, i_phi) + box%cc(i-1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi) - &
               dx2 * box%cc(i, j, i_rhs))
       end do
    end do
#elif $D == 3
    do k = 1, nc
       do j = 1, nc
          i0 = 2 - iand(ieor(redblack_cntr, k+j), 1)
          do i = i0, nc, 2
             box%cc(i, j, k, i_phi) = sixth * ( &
                  box%cc(i+1, j, k, i_phi) + box%cc(i-1, j, k, i_phi) + &
                  box%cc(i, j+1, k, i_phi) + box%cc(i, j-1, k, i_phi) + &
                  box%cc(i, j, k+1, i_phi) + box%cc(i, j, k-1, i_phi) - &
                  dx2 * box%cc(i, j, k, i_rhs))
          end do
       end do
    end do
#endif
  end subroutine mg$Dd_gsrb_lpl_box

  subroutine mg$Dd_lpl_box(box, i_in, i_out)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: i_in, i_out
    integer                     :: i, j, nc
    real(dp)                    :: inv_dr_sq
#if $D == 3
    integer                     :: k
#endif

    nc = box%n_cell
    inv_dr_sq = 1 / box%dr**2

#if $D == 2
    do j = 1, nc
       do i = 1, nc
          box%cc(i, j, i_out) = inv_dr_sq * (box%cc(i-1, j, i_in) + &
               box%cc(i+1, j, i_in) + box%cc(i, j-1, i_in) + &
               box%cc(i, j+1, i_in) - 4 * box%cc(i, j, i_in))
       end do
    end do
#elif $D == 3
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             box%cc(i, j, k, i_out) = inv_dr_sq * (box%cc(i-1, j, k, i_in) + &
                  box%cc(i+1, j, k, i_in) + box%cc(i, j-1, k, i_in) + &
                  box%cc(i, j+1, k, i_in) + box%cc(i, j, k-1, i_in) + &
                  box%cc(i, j, k+1, i_in) - 6 * box%cc(i, j, k, i_in))
          end do
       end do
    end do
#endif
  end subroutine mg$Dd_lpl_box

  subroutine residual_boxes(boxes, ids, mg)
    type(box$D_t), intent(inout) :: boxes(:)
    type(mg$D_t), intent(in)     :: mg
    integer, intent(in)         :: ids(:)
    integer                     :: i

    !$omp do
    do i = 1, size(ids)
       call residual_box(boxes(ids(i)), mg)
    end do
    !$omp end do
  end subroutine residual_boxes

  subroutine residual_box(box, mg)
    type(box$D_t), intent(inout) :: box
    type(mg$D_t), intent(in)     :: mg
    integer                     :: nc

    call mg%box_op(box, mg%i_phi, mg%i_res)
    nc = box%n_cell
#if $D == 2
    box%cc(1:nc, 1:nc, mg%i_res) = box%cc(1:nc, 1:nc, mg%i_rhs) &
         - box%cc(1:nc, 1:nc, mg%i_res)
#elif $D == 3
    box%cc(1:nc, 1:nc, 1:nc, mg%i_res) = box%cc(1:nc, 1:nc, 1:nc, mg%i_rhs) &
         - box%cc(1:nc, 1:nc, 1:nc, mg%i_res)
#endif
  end subroutine residual_box

end module m_mg_$Dd
