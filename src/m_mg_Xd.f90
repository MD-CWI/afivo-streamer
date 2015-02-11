!> Multigrid code for $D-dimensional problems
!> \author Jannis Teunissen
!> \copyright GPLv3

! The following replacements take place on this file (m_mg_Xd.f90) to generate
! 2D and 3D versions:
! 1. $D -> 2 or 3 (dimension of code)
! 2. preprocess file with cpp
! 3. cat -s (merge multiple blank lines)
module m_mg_$Dd
  use m_afivo_$Dd

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  !> Type to store multigrid options in
  type, public :: mg$D_t
     integer :: i_phi        = -1 !< Variable holding solution
     integer :: i_tmp        = -1 !< Internal variable (holding prev. solution)
     integer :: i_rhs        = -1 !< Variable holding right-hand side
     integer :: i_res        = -1 !< Variable holding residual

     integer :: i_eps        = -1 !< Optional variable (diel. permittivity)
     integer :: i_lsf        = -1 !< Optional variable for level set function

     integer :: n_cycle_down = -1 !< Number of relaxation cycles in downward sweep
     integer :: n_cycle_up   = -1 !< Number of relaxation cycles in upward sweep
     integer :: n_cycle_base = -1 !< Number of relaxation cycles at bottom level
     logical :: initialized  = .false.

     !> Routine to call for filling ghost cells near physical boundaries
     procedure(a$D_subr_gc), pointer, nopass   :: sides_bc => null()

     !> Routine to call for filling ghost cells near refinement boundaries
     procedure(a$D_subr_gc), pointer, nopass   :: sides_rb => null()

     !> Subroutine that performs the (non)linear operator
     procedure(mg$D_box_op), pointer, nopass   :: box_op => null()

     !> Subroutine that performs Gauss-Seidel relaxation on a box
     procedure(mg$D_box_gsrb), pointer, nopass :: box_gsrb => null()

     !> Subroutine that corrects the children of a box
     procedure(mg$D_box_corr), pointer, nopass :: box_corr => null()

     !> Todo: restriction method
  end type mg$D_t

  abstract interface
     !> Subroutine that performs A * cc(..., i_in) = cc(..., i_out)
     subroutine mg$D_box_op(box, i_out, mg)
       import
       type(box$D_t), intent(inout) :: box !< The box to operate on
       type(mg$D_t), intent(in)     :: mg
       integer, intent(in)         :: i_out !< Index of output variable
     end subroutine mg$D_box_op

     !> Subroutine that performs Gauss-Seidel relaxation
     subroutine mg$D_box_gsrb(box, redblack_cntr, mg)
       import
       type(box$D_t), intent(inout) :: box !< The box to operate on
       type(mg$D_t), intent(in)     :: mg
       integer, intent(in)         :: redblack_cntr !< Iteration counter
     end subroutine mg$D_box_gsrb

     subroutine mg$D_box_corr(box_p, box_c, mg)
       import
       type(box$D_t), intent(inout) :: box_c
       type(box$D_t), intent(in)    :: box_p
       type(mg$D_t), intent(in)     :: mg
     end subroutine mg$D_box_corr
  end interface

  public :: mg$D_init_mg
  public :: mg$D_fas_fmg
  public :: mg$D_fas_vcycle
  public :: mg$D_set_curvature
  public :: mg$D_box_op
  public :: mg$D_box_gsrb
  public :: mg$D_box_corr
  public :: mg$D_box_lpl
  public :: mg$D_box_gsrb_lpl
  public :: mg$D_box_corr_lpl

contains

  !> Check multigrid options or set them to default
  subroutine mg$D_init_mg(mg)
    type(mg$D_t), intent(inout) :: mg           !< Multigrid options

    if (mg%i_phi < 0)                  stop "mg$D_init_mg: i_phi not set"
    if (mg%i_tmp < 0)                  stop "mg$D_init_mg: i_tmp not set"
    if (mg%i_rhs < 0)                  stop "mg$D_init_mg: i_rhs not set"
    if (mg%i_res < 0)                  stop "mg$D_init_mg: i_res not set"

    if (.not. associated(mg%sides_bc)) stop "mg$D_init_mg: sides_bc not set"

    ! Check whether these are set, otherwise use default
    if (mg%n_cycle_down < 0)           mg%n_cycle_down = 2
    if (mg%n_cycle_up < 0)             mg%n_cycle_up = 2
    if (mg%n_cycle_base < 0)           mg%n_cycle_base = 2

    ! Check whether methods are set, otherwise use default (for laplacian)
    if (.not. associated(mg%sides_rb)) mg%sides_rb => a$D_sides_extrap
    if (.not. associated(mg%box_op))   mg%box_op => mg$D_box_lpl
    if (.not. associated(mg%box_gsrb)) mg%box_gsrb => mg$D_box_gsrb_lpl
    if (.not. associated(mg%box_corr)) mg%box_corr => mg$D_box_corr_lpl

    mg%initialized = .true.
  end subroutine mg$D_init_mg

  subroutine check_mg(mg)
    type(mg$D_t), intent(in) :: mg           !< Multigrid options
    if (.not. mg%initialized) stop "check_mg: you haven't called mg$D_init_mg"
  end subroutine check_mg

  !> Perform FAS-FMG cycle (full approximation scheme, full multigrid). Note
  !> that this routine needs valid ghost cells (for i_phi) on input, and gives
  !> back valid ghost cells on output
  subroutine mg$D_fas_fmg(tree, mg)
    type(a$D_t), intent(inout)       :: tree !< Tree to do multigrid on
    type(mg$D_t), intent(in)         :: mg   !< Multigrid options
    integer                         :: lvl, min_lvl

    call check_mg(mg)           ! Check whether mg options are set
    min_lvl = lbound(tree%lvls, 1)

    do lvl = min_lvl, tree%max_lvl
       ! Store phi_old in tmp
       call a$D_boxes_copy_cc(tree%boxes, tree%lvls(lvl)%ids, &
            mg%i_phi, mg%i_tmp)

       if (lvl > min_lvl) then
          ! Correct solution at this lvl using lvl-1 data
          ! phi = phi + prolong(phi_coarse - phi_old_coarse)
          call correct_children(tree%boxes, tree%lvls(lvl-1)%parents, mg)

          ! Update ghost cells
          call fill_gc_phi(tree%boxes, tree%lvls(lvl)%ids, mg)
       end if

       ! Perform V-cycle
       call mg$D_fas_vcycle(tree, mg, lvl)
    end do
  end subroutine mg$D_fas_fmg

  !> Perform FAS V-cycle (full approximation scheme). Note that this routine
  !> needs valid ghost cells (for i_phi) on input, and gives back valid ghost
  !> cells on output
  subroutine mg$D_fas_vcycle(tree, mg, max_lvl)
    type(a$D_t), intent(inout) :: tree !< Tree to do multigrid on
    type(mg$D_t), intent(in)   :: mg   !< Multigrid options
    integer, intent(in)        :: max_lvl !< Maximum level for V-cycle
    integer                    :: i, id, lvl, min_lvl

    call check_mg(mg)           ! Check whether mg options are set
    min_lvl = lbound(tree%lvls, 1)

    do lvl = max_lvl,  min_lvl+1, -1
       ! Downwards relaxation
       call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_down)

       ! Calculate residual at current lvl
       call residual_boxes(tree%boxes, tree%lvls(lvl)%ids, mg)

       call a$D_restrict_to_boxes(tree%boxes, tree%lvls(lvl-1)%parents, mg%i_phi)
       call a$D_restrict_to_boxes(tree%boxes, tree%lvls(lvl-1)%parents, mg%i_res)

       call fill_gc_phi(tree%boxes, tree%lvls(lvl-1)%ids, mg)

       ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined, and
       ! store current coarse phi in tmp.

       !$omp do private(id)
       do i = 1, size(tree%lvls(lvl-1)%parents)
          id = tree%lvls(lvl-1)%parents(i)
          call mg%box_op(tree%boxes(id), mg%i_rhs, mg)
          call a$D_box_add_cc(tree%boxes(id), mg%i_res, mg%i_rhs)
          call a$D_box_copy_cc(tree%boxes(id), mg%i_phi, mg%i_tmp)
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
       call fill_gc_phi(tree%boxes, tree%lvls(lvl)%ids, mg)

       ! Upwards relaxation
       call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_up)
    end do
    !$omp barrier
  end subroutine mg$D_fas_vcycle

  !> Set variable i_crv to an estimate of the curvature of i_phi
  subroutine mg$D_set_curvature(tree, i_crv, mg)
    type(a$D_t), intent(inout)  :: tree !< Tree to do multigrid on
    integer, intent(in)        :: i_crv
    type(mg$D_t), intent(in)    :: mg   !< Multigrid options
    integer                    :: i, id, lvl, min_lvl
    real(dp)                   :: dr2

    min_lvl = lbound(tree%lvls, 1)

    do lvl = min_lvl, tree%max_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          dr2 = tree%boxes(id)%dr**2
          call mg%box_op(tree%boxes(id), i_crv, mg)
          call a$D_box_times_cc(tree%boxes(id), dr2, i_crv)
       end do
       !$omp end do nowait
    end do
    !$omp barrier
  end subroutine mg$D_set_curvature

  subroutine fill_gc_phi(boxes, ids, mg)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:)
    type(mg$D_t), intent(in)     :: mg
    integer                     :: i

    !$omp do
    do i = 1, size(ids)
       call a$D_gc_box_sides(boxes, ids(i), mg%i_phi, &
            mg%sides_rb, mg%sides_bc)
    end do
    !$omp end do
  end subroutine fill_gc_phi

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

       ! Store the correction in i_tmp
#if $D == 2
       boxes(id)%cc(:, :, mg%i_tmp) = boxes(id)%cc(:, :, mg%i_phi) - &
            boxes(id)%cc(:, :, mg%i_tmp)
#elif $D == 3
       boxes(id)%cc(:, :, :, mg%i_tmp) = boxes(id)%cc(:, :, :, mg%i_phi) - &
            boxes(id)%cc(:, :, :, mg%i_tmp)
#endif

       do i_c = 1, a$D_num_children
          c_id = boxes(id)%children(i_c)
          if (c_id == a5_no_box) cycle
          call mg%box_corr(boxes(id), boxes(c_id), mg)
       end do
    end do
    !$omp end do
  end subroutine correct_children

  subroutine mg$D_box_corr_lpl(box_p, box_c, mg)
    type(box$D_t), intent(inout) :: box_c
    type(box$D_t), intent(in)    :: box_p
    type(mg$D_t), intent(in)     :: mg
    integer                     :: ix_offset($D), i_phi, i_corr
    integer                     :: nc, i, j, i_c1, i_c2, j_c1, j_c2
#if $D == 3
    integer                     :: k, k_c1, k_c2
#endif

    nc        = box_c%n_cell
    ix_offset = a$D_get_child_offset(box_c)
    i_phi     = mg%i_phi
    i_corr    = mg%i_tmp

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
               + 0.5_dp * box_p%cc(i_c1, j_c1, i_corr) &
               + 0.25_dp * (box_p%cc(i_c2, j_c1, i_corr) &
               +       box_p%cc(i_c1, j_c2, i_corr))
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
  end subroutine mg$D_box_corr_lpl

  subroutine gsrb_boxes(boxes, ids, mg, n_cycle)
    type(box$D_t), intent(inout) :: boxes(:)
    type(mg$D_t), intent(in)     :: mg
    integer, intent(in)         :: ids(:), n_cycle
    integer                     :: n, i

    do n = 1, 2 * n_cycle
       !$omp do
       do i = 1, size(ids)
          call mg%box_gsrb(boxes(ids(i)), n, mg)
       end do
       !$omp end do

       !$omp do
       do i = 1, size(ids)
          call a$D_gc_box_sides(boxes, ids(i), mg%i_phi, &
               mg%sides_rb, mg%sides_bc)
       end do
       !$omp end do
    end do
  end subroutine gsrb_boxes

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator
  subroutine mg$D_box_gsrb_lpl(box, redblack_cntr, mg)
    type(box$D_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    type(mg$D_t), intent(in)     :: mg
    integer                     :: i, i0, j, nc, i_phi, i_rhs
    real(dp)                    :: dx2
#if $D == 3
    integer                     :: k
    real(dp), parameter         :: sixth = 1/6.0_dp
#endif

    dx2   = box%dr**2
    nc    = box%n_cell
    i_phi = mg%i_phi
    i_rhs = mg%i_rhs

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
  end subroutine mg$D_box_gsrb_lpl

  !> Perform Laplacian operator on a box
  subroutine mg$D_box_lpl(box, i_out, mg)
    type(box$D_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    type(mg$D_t), intent(in)     :: mg
    integer                     :: i, j, nc, i_phi
    real(dp)                    :: inv_dr_sq
#if $D == 3
    integer                     :: k
#endif

    nc = box%n_cell
    inv_dr_sq = 1 / box%dr**2
    i_phi = mg%i_phi

#if $D == 2
    do j = 1, nc
       do i = 1, nc
          box%cc(i, j, i_out) = inv_dr_sq * (box%cc(i-1, j, i_phi) + &
               box%cc(i+1, j, i_phi) + box%cc(i, j-1, i_phi) + &
               box%cc(i, j+1, i_phi) - 4 * box%cc(i, j, i_phi))
       end do
    end do
#elif $D == 3
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             box%cc(i, j, k, i_out) = inv_dr_sq * (box%cc(i-1, j, k, i_phi) + &
                  box%cc(i+1, j, k, i_phi) + box%cc(i, j-1, k, i_phi) + &
                  box%cc(i, j+1, k, i_phi) + box%cc(i, j, k-1, i_phi) + &
                  box%cc(i, j, k+1, i_phi) - 6 * box%cc(i, j, k, i_phi))
          end do
       end do
    end do
#endif
  end subroutine mg$D_box_lpl

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

    call mg%box_op(box, mg%i_res, mg)
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
