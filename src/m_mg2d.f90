module m_mg2d
  use m_afivo

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  type mg2_subtree_t
     integer :: n_lvls
     
  end type mg2_subtree_t

  type mg2_t
     private
     integer  :: i_phi
     integer  :: i_phi_old
     integer  :: i_rhs
     integer  :: i_res
     integer  :: n_cycle_down
     integer  :: n_cycle_up
     real(dp) :: gmres_reduction = 1.0e-10_dp
     type(mg2_subtree_t) :: subtree
     procedure(a2_subr_gc), pointer, nopass :: sides_bc, corners_bc
  end type mg2_t

  ! Public types
  public :: mg2_t

  ! Public methods
  public :: mg2d_set
  public :: mg2d_restrict_rhs
  public :: mg2d_fas_fmg
  public :: mg2d_fas_vcycle

contains

  subroutine mg2d_set(mg, i_phi, i_phi_old, i_rhs, i_res, &
       n_cycle_down, n_cycle_up, sides_bc, corners_bc)
    type(mg2_t), intent(out) :: mg
    integer, intent(in)      :: i_phi, i_phi_old, i_rhs, i_res
    integer, intent(in)      :: n_cycle_down, n_cycle_up
    procedure(a2_subr_gc)    :: sides_bc, corners_bc
    mg%i_phi        = i_phi
    mg%i_phi_old    = i_phi_old
    mg%i_rhs        = i_rhs
    mg%i_res        = i_res
    mg%n_cycle_down = n_cycle_down
    mg%n_cycle_up   = n_cycle_up
    mg%sides_bc     => sides_bc
    mg%corners_bc   => corners_bc
  end subroutine mg2d_set

  subroutine mg2d_restrict_rhs(tree, mg)
    type(a2_t), intent(inout) :: tree
    type(mg2_t), intent(inout) :: mg
    integer :: lvl, i, id

    do lvl = tree%n_lvls-1, 1, -1
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call a2_restrict_to(tree%boxes, id, [mg%i_rhs])
       end do
    end do
  end subroutine mg2d_restrict_rhs

  ! Need valid ghost cells on input, has valid gc on output
  subroutine mg2d_fas_fmg(tree, mg)
    type(a2_t), intent(inout) :: tree
    type(mg2_t), intent(in)   :: mg
    integer                   :: i, id, lvl

    do lvl = 1, tree%n_lvls
       ! Store phi in phi_old
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          tree%boxes(id)%cc(:, :, mg%i_phi_old) = &
               tree%boxes(id)%cc(:, :, mg%i_phi)
       end do

       if (lvl > 1) then
          ! Correct solution at this lvl using lvl-1 data
          ! phi = phi + prolong(phi_coarse - phi_old_coarse)
          do i = 1, size(tree%lvls(lvl-1)%parents)
             id = tree%lvls(lvl-1)%parents(i)
             call correct_from_box(tree%boxes, id, mg%i_phi, mg%i_phi_old)
          end do

          ! Update ghost cells
          call mg2d_fill_gc(tree%boxes, tree%lvls(lvl)%ids, [mg%i_phi], &
               mg%sides_bc, mg%corners_bc)
       endif

       call mg2d_fas_vcycle(tree, mg, lvl) ! Perform V-cycle
    end do
  end subroutine mg2d_fas_fmg

  ! On entrance, need valid ghost cell data. On exit, leave valid ghost cell data
  recursive subroutine mg2d_fas_vcycle(tree, mg, lvl)
    type(a2_t), intent(inout) :: tree
    type(mg2_t), intent(in)   :: mg
    integer, intent(in)       :: lvl
    integer                   :: i, id, n

    if (lvl == 1) then
       do n = 1, 2 * (mg%n_cycle_down + mg%n_cycle_up)
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             call gsrb_box(tree%boxes(id), mg%i_phi, mg%i_rhs, n)
          end do

          ! Communicate updated boundary cells
          call mg2d_fill_gc(tree%boxes, tree%lvls(lvl)%ids, [mg%i_phi], &
               mg%sides_bc, mg%corners_bc)
       end do
    else
       ! Downwards relaxation. Half the cycle are "red", half are "black".
       do n = 1, 2 * mg%n_cycle_down
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             call gsrb_box(tree%boxes(id), mg%i_phi, mg%i_rhs, n)
          end do

          ! Communicate updated boundary cells
          call mg2d_fill_gc(tree%boxes, tree%lvls(lvl)%ids, [mg%i_phi], &
               mg%sides_bc, mg%corners_bc)
       end do

       ! Calculate residual at current lvl
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call residual_box(tree%boxes(id), mg%i_res, mg%i_phi, mg%i_rhs)
       end do

       ! Restrict phi
       do i = 1, size(tree%lvls(lvl-1)%parents)
          id = tree%lvls(lvl-1)%parents(i)
          call a2_restrict_to(tree%boxes, id, [mg%i_phi])
       end do

       ! Have to update ghost cells for phi_c (todo: not everywhere?)
       call mg2d_fill_gc(tree%boxes, tree%lvls(lvl-1)%ids, [mg%i_phi], &
            mg%sides_bc, mg%corners_bc)

       ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined
       do i = 1, size(tree%lvls(lvl-1)%parents)
          id = tree%lvls(lvl-1)%parents(i)
          call laplacian_box(tree%boxes(id), mg%i_rhs, mg%i_phi)
          call a2_restrict_to_iadd(tree%boxes, id, mg%i_res, mg%i_rhs)
       end do

       ! Store current coarse phi in phi_old
       do i = 1, size(tree%lvls(lvl-1)%ids)
          id = tree%lvls(lvl-1)%ids(i)
          tree%boxes(id)%cc(:, :, mg%i_phi_old) = &
               tree%boxes(id)%cc(:, :, mg%i_phi)
       end do

       call mg2d_fas_vcycle(tree, mg, lvl-1)

       ! Correct solution at this lvl using lvl-1 data
       ! phi = phi + prolong(phi_coarse - phi_old_coarse)
       do i = 1, size(tree%lvls(lvl-1)%parents)
          id = tree%lvls(lvl-1)%parents(i)
          call correct_from_box(tree%boxes, id, mg%i_phi, mg%i_phi_old)
       end do

       ! Have to fill ghost cells again (todo: not everywhere?)
       call mg2d_fill_gc(tree%boxes, tree%lvls(lvl)%ids, [mg%i_phi], &
            mg%sides_bc, mg%corners_bc)

       ! Upwards relaxation
       do n = 1, 2 * mg%n_cycle_down
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             call gsrb_box(tree%boxes(id), mg%i_phi, mg%i_rhs, n)
          end do

          ! Communicate updated boundary cells
          call mg2d_fill_gc(tree%boxes, tree%lvls(lvl)%ids, [mg%i_phi], &
               mg%sides_bc, mg%corners_bc)
       end do
    end if

    ! Calculate residual at current lvl for output
    do i = 1, size(tree%lvls(lvl)%ids)
       id = tree%lvls(lvl)%ids(i)
       call residual_box(tree%boxes(id), mg%i_res, mg%i_phi, mg%i_rhs)
    end do
  end subroutine mg2d_fas_vcycle

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
  subroutine correct_from_box(boxes, id, i_phi, i_phi_old)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, i_phi, i_phi_old
    real(dp), parameter         :: f1=1/16.0_dp, f3=3/16.0_dp, f9=9/16.0_dp
    integer                     :: nc, i_c, c_id, ix_offset(2)
    integer                     :: i, j, i_c1, i_c2, j_c1, j_c2

    nc = boxes(id)%cfg%n_cell
    do i_c = 1, 4
       c_id = boxes(id)%children(i_c)

       ! Offset of child w.r.t. parent
       ix_offset = a2_ch_dix(:, i_c) * ishft(nc, -1)

       ! In these loops, we calculate the closest coarse index (_c1), and the
       ! one-but-closest (_c2). The fine cell lies in between.
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

             boxes(c_id)%cc(i, j, i_phi) = boxes(c_id)%cc(i, j, i_phi) &
                  + f9 * (boxes(id)%cc(i_c1, j_c1, i_phi) &
                  - boxes(id)%cc(i_c1, j_c1, i_phi_old)) &
                  + f3 * (boxes(id)%cc(i_c2, j_c1, i_phi) &
                  - boxes(id)%cc(i_c2, j_c1, i_phi_old) &
                  + boxes(id)%cc(i_c1, j_c2, i_phi) &
                  - boxes(id)%cc(i_c1, j_c2, i_phi_old)) &
                  + f1 * (boxes(id)%cc(i_c2, j_c2, i_phi) &
                  - boxes(id)%cc(i_c2, j_c2, i_phi_old))
          end do
       end do
    end do
  end subroutine correct_from_box

  subroutine gsrb_box(box, i_phi, i_rhs, redblack_cntr)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: i_phi, i_rhs, redblack_cntr
    integer                     :: i, j, nc, di(2)
    real(dp)                    :: dxdy

    dxdy = product(a2_dr(box))
    nc   = box%cfg%n_cell

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

  subroutine laplacian_box(box, i_lpl, i_phi)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: i_lpl, i_phi
    integer                     :: i, j, nc
    real(dp)                    :: inv_dr_sq(2)

    nc = box%cfg%n_cell
    inv_dr_sq = 1 / a2_dr(box)**2

    do j = 1, nc
       do i = 1, nc
          box%cc(i, j, i_lpl) = inv_dr_sq(1) * (box%cc(i-1, j, i_phi) &
               - 2 * box%cc(i, j, i_phi) + box%cc(i+1, j, i_phi)) + &
               inv_dr_sq(2) * (box%cc(i, j-1, i_phi) &
               - 2 * box%cc(i, j, i_phi) + box%cc(i, j+1, i_phi))
       end do
    end do
  end subroutine laplacian_box

  subroutine residual_box(box, i_res, i_phi, i_rhs)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: i_res, i_phi, i_rhs
    integer                     :: i, j, nc
    real(dp)                    :: inv_dr_sq(2)

    nc = box%cfg%n_cell
    inv_dr_sq = 1 / a2_dr(box)**2

    do j = 1, nc
       do i = 1, nc
          box%cc(i, j, i_res) = box%cc(i, j, i_rhs) - ( &
               inv_dr_sq(1) * (box%cc(i-1, j, i_phi) &
               - 2 * box%cc(i, j, i_phi) + box%cc(i+1, j, i_phi)) + &
               inv_dr_sq(2) * (box%cc(i, j-1, i_phi) &
               - 2 * box%cc(i, j, i_phi) + box%cc(i, j+1, i_phi)))
       end do
    end do
  end subroutine residual_box

end module m_mg2d