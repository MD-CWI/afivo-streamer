!> Multigrid code for $D-dimensional problems
!> \author Jannis Teunissen
!> \copyright GPLv3

! The following replacements take place on this file (m_mgb_Xd.f90) to generate
! 2D and 3D versions:
! 1. $D -> 2 or 3 (dimension of code)
! 2. preprocess file with cpp
! 3. cat -s (merge multiple blank lines)
module m_mgb_$Dd
  use m_afivo_$Dd

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  !> Type to store multigrid options in
  type, public :: mg$D_t
     integer :: i_phi           !< Variable holding solution
     integer :: i_phi_old       !< Internal variable (buffer for prev. solution)
     integer :: i_rhs           !< Variable holding right-hand side
     integer :: i_res           !< Variable holding residual
     integer :: i_lsf           !< Variable holding level set function
     integer :: n_cycle_down    !< Number of relaxation cycles in downward sweep
     integer :: n_cycle_up      !< Number of relaxation cycles in upward sweep
     integer :: n_cycle_base    !< Number of relaxation cycles at bottom level
     !> Routine to call for filling ghost cells near physical boundaries
     procedure(a$D_subr_gc), pointer, nopass    :: sides_bc
     !> Subroutine that performs the (non)linear operator
     procedure(mg$Dd_box_op), pointer, nopass   :: box_op
     !> Subroutine that performs Gauss-Seidel relaxation on a box
     procedure(mg$Dd_box_gsrb), pointer, nopass :: box_gsrb
  end type mg$D_t

  abstract interface
     !> Subroutine that performs A * cc(..., i_in) = cc(..., i_out)
     subroutine mg$Dd_box_op(box, i_in, i_out, i_lsf)
       import
       type(box$D_t), intent(inout) :: box !< The box to operate on
       integer, intent(in)         :: i_in !< Index of input variable
       integer, intent(in)         :: i_out !< Index of output variable
       integer, intent(in)         :: i_lsf !< Index of level set function
     end subroutine mg$Dd_box_op

     !> Subroutine that performs Gauss-Seidel relaxation
     subroutine mg$Dd_box_gsrb(box, i_phi, i_rhs, i_lsf, redblack_cntr)
       import
       type(box$D_t), intent(inout) :: box !< The box to operate on
       integer, intent(in)         :: i_phi !< Index of solution variable
       integer, intent(in)         :: i_rhs !< Index of right-hand side variable
       integer, intent(in)         :: i_lsf !< Index of level set function
       integer, intent(in)         :: redblack_cntr !< Iteration counter
     end subroutine mg$Dd_box_gsrb
  end interface

  public :: mg$Dd_box_op
  public :: mg$Dd_box_gsrb
  public :: mg$Dd_set
  public :: mg$Dd_fas_fmg
  public :: mg$Dd_fas_vcycle
  public :: mg$Dd_lpl_box
  public :: mg$Dd_gsrb_lpl_box

contains

  !> Store multigrid options in a mg type
  subroutine mg$Dd_set(mg, i_phi, i_phi_old, i_rhs, i_res, i_lsf, &
       n_cycle_down, n_cycle_up, n_cycle_base, &
       sides_bc, my_operator, my_gsrb)
    type(mg$D_t), intent(out) :: mg           !< Store multigrid options in this variable
    integer, intent(in)       :: i_phi        !< Variable holding solution
    integer, intent(in)       :: i_phi_old    !< Internal variable (holding prev. solution)
    integer, intent(in)       :: i_rhs        !< Variable holding right-hand side
    integer, intent(in)       :: i_res        !< Variable holding residual
    integer, intent(in)       :: i_lsf        !< Variable holding residual
    integer, intent(in)       :: n_cycle_down !< Number of relaxation cycles in downward sweep
    integer, intent(in)       :: n_cycle_up   !< Number of relaxation cycles in upward sweep
    integer, intent(in)       :: n_cycle_base !< Number of relaxation cycles at bottom level
    !> Routine to call for filling ghost cells near physical boundaries
    procedure(a$D_subr_gc)    :: sides_bc
    !> Subroutine that performs the (non)linear operator
    procedure(mg$Dd_box_op)   :: my_operator
    !> Subroutine that performs Gauss-Seidel relaxation on a box
    procedure(mg$Dd_box_gsrb) :: my_gsrb

    mg%i_phi        = i_phi
    mg%i_phi_old    = i_phi_old
    mg%i_rhs        = i_rhs
    mg%i_res        = i_res
    mg%i_lsf        = i_lsf
    mg%n_cycle_down = n_cycle_down
    mg%n_cycle_up   = n_cycle_up
    mg%n_cycle_base = n_cycle_base
    mg%sides_bc     => sides_bc
    mg%box_op       => my_operator
    mg%box_gsrb     => my_gsrb
  end subroutine mg$Dd_set

  !> Perform FAS-FMG cycle (full approximation scheme, full multigrid). Note
  !> that this routine needs valid ghost cells (for i_phi) on input, and gives
  !> back valid ghost cells on output
  subroutine mg$Dd_fas_fmg(tree, mg)
    type(a$D_t), intent(inout)       :: tree !< Tree to do multigrid on
    type(mg$D_t), intent(in)         :: mg   !< Multigrid options
    integer                         :: lvl, min_lvl

    min_lvl = lbound(tree%lvls, 1)

    do lvl = min_lvl, tree%max_lvl
       ! Store phi in phi_old
       call a$D_boxes_copy_cc(tree%boxes, tree%lvls(lvl)%ids, &
            mg%i_phi, mg%i_phi_old)

       if (lvl > min_lvl) then
          ! Correct solution at this lvl using lvl-1 data
          ! phi = phi + prolong(phi_coarse - phi_old_coarse)
          call correct_children(tree%boxes, tree%lvls(lvl-1)%parents, mg)

          ! Update ghost cells
          call fill_gc(tree%boxes, tree%lvls(lvl)%ids, mg%i_phi, &
               mg%sides_bc)
       end if

       ! Perform V-cycle
       call mg$Dd_fas_vcycle(tree, mg, lvl)
    end do

  end subroutine mg$Dd_fas_fmg

  !> Perform FAS V-cycle (full approximation scheme). Note that this routine
  !> needs valid ghost cells (for i_phi) on input, and gives back valid ghost
  !> cells on output
  subroutine mg$Dd_fas_vcycle(tree, mg, max_lvl)
    type(a$D_t), intent(inout) :: tree !< Tree to do multigrid on
    type(mg$D_t), intent(in)   :: mg   !< Multigrid options
    integer, intent(in)        :: max_lvl !< Maximum level for V-cycle
    integer                    :: i, id, lvl, min_lvl

    min_lvl = lbound(tree%lvls, 1)

    do lvl = max_lvl,  min_lvl+1, -1
       ! Downwards relaxation
       call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_down)

       ! Calculate residual at current lvl
       call residual_boxes(tree%boxes, tree%lvls(lvl)%ids, mg)

       call a$D_restrict_to_boxes(tree%boxes, tree%lvls(lvl-1)%parents, mg%i_phi)
       call a$D_restrict_to_boxes(tree%boxes, tree%lvls(lvl-1)%parents, mg%i_res)

       call fill_gc(tree%boxes, tree%lvls(lvl-1)%ids, mg%i_phi, &
            mg%sides_bc)

       ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined, and
       ! store current coarse phi in phi_old
       !$omp do private(id)
       do i = 1, size(tree%lvls(lvl-1)%parents)
          id = tree%lvls(lvl-1)%parents(i)
          call mg%box_op(tree%boxes(id), mg%i_phi, mg%i_rhs, mg%i_lsf)
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
       call fill_gc(tree%boxes, tree%lvls(lvl)%ids, mg%i_phi, &
            mg%sides_bc)

       ! Upwards relaxation
       call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_up)
    end do
  end subroutine mg$Dd_fas_vcycle

  subroutine fill_gc(boxes, ids, iv, sides_bc)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:), iv
    procedure(a$D_subr_gc)       :: sides_bc
    integer                     :: i

    !$omp do
    do i = 1, size(ids)
       call a$D_gc_box_sides(boxes, ids(i), iv, &
            a$D_sides_extrap, sides_bc)
    end do
    !$omp end do
  end subroutine fill_gc

  ! Sets phi = phi + prolong(phi_coarse - phi_old_coarse)
  subroutine correct_children(boxes, ids, mg)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:)
    type(mg$D_t), intent(in)     :: mg
    integer                     :: i, id, nc, i_c, c_id, ix_offset($D)

    !$omp do private(nc, i_c, c_id, ix_offset)
    do i = 1, size(ids)
       id = ids(i)
       nc = boxes(id)%n_cell

       ! Store the correction in i_phi_old
#if $D == 2
       boxes(id)%cc(:, :, mg%i_phi_old) = boxes(id)%cc(:, :, mg%i_phi) - &
            boxes(id)%cc(:, :, mg%i_phi_old)
#elif $D == 3
       boxes(id)%cc(:, :, :, mg%i_phi_old) = boxes(id)%cc(:, :, :, mg%i_phi) - &
            boxes(id)%cc(:, :, :, mg%i_phi_old)
#endif

       do i_c = 1, a$D_num_children
          c_id = boxes(id)%children(i_c)
          if (c_id == a5_no_box) cycle

          ! Offset of child w.r.t. parent
          ix_offset = a$D_ch_dix(:, i_c) * ishft(nc, -1)

          call correct_child_box(boxes(id), boxes(c_id), &
               ix_offset, mg%i_phi, mg%i_phi_old, mg%i_lsf)
       end do
    end do
    !$omp end do
  end subroutine correct_children

  subroutine correct_child_box(box_p, box_c, ix_offset, i_phi, i_corr, i_lsf)
    type(box$D_t), intent(inout) :: box_c
    type(box$D_t), intent(in)    :: box_p
    integer, intent(in)         :: i_phi, i_corr, i_lsf, ix_offset($D)
    integer                     :: nc, i, j, i_c1, i_c2, j_c1, j_c2
#if $D == 3
    integer                     :: k, k_c1, k_c2
#endif
    real(dp) :: lsf, val(3), dist(3), tmp

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

          call lsf_dist_val(box_c%cc(i, j, i_lsf), &
               box_p%cc(i_c1, j_c1, [i_corr, i_lsf]), 0.0_dp, dist(1), val(1))
          call lsf_dist_val(box_c%cc(i, j, i_lsf), &
               box_p%cc(i_c2, j_c1, [i_corr, i_lsf]), 0.0_dp, dist(2), val(2))
          call lsf_dist_val(box_c%cc(i, j, i_lsf), &
               box_p%cc(i_c1, j_c2, [i_corr, i_lsf]), 0.0_dp, dist(3), val(3))

          tmp = 1 / (dist(1) * dist(2) + dist(1) * dist(3) + 2 * dist(2) * dist(3))
          box_c%cc(i, j, i_phi) = box_c%cc(i, j, i_phi) + &
               2 * dist(2) * dist(3) * tmp * val(1) + &
               dist(1) * dist(3) * tmp * val(2) + &
               dist(1) * dist(2) * tmp * val(3)
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

             box_c%cc(i, j, k, i_phi) = box_c%cc(i, j, k, i_phi) + 0.25_dp * ( &
                  box_p%cc(i_c1, j_c1, k_c1, i_corr) &
                  + box_p%cc(i_c2, j_c1, k_c1, i_corr) &
                  + box_p%cc(i_c1, j_c2, k_c1, i_corr) &
                  + box_p%cc(i_c1, j_c1, k_c2, i_corr))
          end do
       end do
    end do
#endif
  end subroutine correct_child_box

  subroutine lsf_dist_val(lsf_a, v_b, b_value, dist, val)
    real(dp), intent(in)  :: lsf_a, v_b(2), b_value
    real(dp), intent(out) :: dist, val
    real(dp)              :: tmp

    ! Determine whether there is a boundary
    if (lsf_a * v_b(2) < 0) then
       dist = lsf_a / (lsf_a - v_b(2))
       val  = b_value
    else
       dist = 1
       val  = v_b(1)
    end if
  end subroutine lsf_dist_val

  subroutine gsrb_boxes(boxes, ids, mg, n_cycle)
    type(box$D_t), intent(inout) :: boxes(:)
    type(mg$D_t), intent(in)     :: mg
    integer, intent(in)         :: ids(:), n_cycle
    integer                     :: n, i

    do n = 1, 2 * n_cycle
       !$omp do
       do i = 1, size(ids)
          call mg%box_gsrb(boxes(ids(i)), mg%i_phi, mg%i_rhs, mg%i_lsf, n)
       end do
       !$omp end do

       !$omp do
       do i = 1, size(ids)
          call a$D_gc_box_sides(boxes, ids(i), mg%i_phi, &
               a$D_sides_extrap, mg%sides_bc)
       end do
       !$omp end do
    end do
  end subroutine gsrb_boxes

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator
  subroutine mg$Dd_gsrb_lpl_box(box, i_phi, i_rhs, i_lsf, redblack_cntr)
    type(box$D_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: i_phi !< Index of solution variable
    integer, intent(in)         :: i_rhs !< Index of right-hand side
    integer, intent(in)         :: i_lsf
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    integer                     :: i, i0, j, nc
    real(dp)                    :: dx2, dist(4), val(4)
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
          call lsf_dist_val(box%cc(i, j, i_lsf), &
               box%cc(i-1, j, [i_phi, i_lsf]), 0.0_dp, dist(1), val(1))
          call lsf_dist_val(box%cc(i, j, i_lsf), &
               box%cc(i+1, j, [i_phi, i_lsf]), 0.0_dp, dist(2), val(2))
          call lsf_dist_val(box%cc(i, j, i_lsf), &
               box%cc(i, j-1, [i_phi, i_lsf]), 0.0_dp, dist(3), val(3))
          call lsf_dist_val(box%cc(i, j, i_lsf), &
               box%cc(i, j+1, [i_phi, i_lsf]), 0.0_dp, dist(4), val(4))

          box%cc(i, j, i_phi) = 0.5_dp / &
               (dist(1) * dist(2) + dist(3) * dist(4)) * ( &
               (dist(2) * val(1) + dist(1) * val(2)) * &
               dist(3) * dist(4) / (0.5 * (dist(1) + dist(2))) + &
               (dist(4) * val(3) + dist(3) * val(4)) * &
               dist(1) * dist(2) / (0.5 * (dist(3) + dist(4))) - &
               product(dist) * dx2 * box%cc(i, j, i_rhs))
       end do
    end do
#elif $D == 3
    stop
#endif
  end subroutine mg$Dd_gsrb_lpl_box

  !> Perform Laplacian operator on a box
  subroutine mg$Dd_lpl_box(box, i_in, i_out, i_lsf)
    type(box$D_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: i_in !< Index of variable to take Laplacian of
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    integer, intent(in)         :: i_lsf
    integer                     :: i, j, nc
    real(dp)                    :: inv_dr_sq, dist(4), val(4), f0
#if $D == 3
    integer                     :: k
#endif

    nc = box%n_cell
    inv_dr_sq = 1 / box%dr**2

#if $D == 2
    do j = 1, nc
       do i = 1, nc
          call lsf_dist_val(box%cc(i, j, i_lsf), &
               box%cc(i-1, j, [i_in, i_lsf]), 0.0_dp, dist(1), val(1))
          call lsf_dist_val(box%cc(i, j, i_lsf), &
               box%cc(i+1, j, [i_in, i_lsf]), 0.0_dp, dist(2), val(2))
          call lsf_dist_val(box%cc(i, j, i_lsf), &
               box%cc(i, j-1, [i_in, i_lsf]), 0.0_dp, dist(3), val(3))
          call lsf_dist_val(box%cc(i, j, i_lsf), &
               box%cc(i, j+1, [i_in, i_lsf]), 0.0_dp, dist(4), val(4))

          f0 = box%cc(i, j, i_in)
          box%cc(i, j, i_out) = inv_dr_sq * ( &
               (dist(2) * val(1) + dist(1) * val(2) - (dist(1)+dist(2)) * f0) / &
               (0.5_dp * (dist(1) + dist(2)) * dist(1) * dist(2)) + &
               (dist(4) * val(3) + dist(3) * val(4) - (dist(3)+dist(4)) * f0) / &
               (0.5_dp * (dist(3) + dist(4)) * dist(3) * dist(4)))
       end do
    end do
#elif $D == 3
    stop
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

    call mg%box_op(box, mg%i_phi, mg%i_res, mg%i_lsf)
    nc = box%n_cell
#if $D == 2
    box%cc(1:nc, 1:nc, mg%i_res) = box%cc(1:nc, 1:nc, mg%i_rhs) &
         - box%cc(1:nc, 1:nc, mg%i_res)
#elif $D == 3
    box%cc(1:nc, 1:nc, 1:nc, mg%i_res) = box%cc(1:nc, 1:nc, 1:nc, mg%i_rhs) &
         - box%cc(1:nc, 1:nc, 1:nc, mg%i_res)
#endif
  end subroutine residual_box

end module m_mgb_$Dd
