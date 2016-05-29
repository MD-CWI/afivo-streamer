! This module contains the multigrid routines that come with Afivo
!
! Author: Jannis Teunissen
! License: GPLv3

module m_a$D_multigrid
  use m_a$D_types

  implicit none
  private

  ! The mg module supports different multigrid operators, and uses these tags to
  ! identify boxes / operators
  integer, parameter, public :: mg_normal_box = 1 !< Normal box
  integer, parameter, public :: mg_lsf_box    = 2 !< Box with an internal boundary
  integer, parameter, public :: mg_ceps_box   = 3 !< Box with constant eps /= 1
  integer, parameter, public :: mg_veps_box   = 4 !< Box with varying eps (on face)

  !> Type to store multigrid options in
  type, public :: mg$D_t
     integer :: i_phi        = -1 !< Variable holding solution
     integer :: i_rhs        = -1 !< Variable holding right-hand side
     integer :: i_tmp        = -1 !< Internal variable (holding prev. solution)

     integer :: i_eps        = -1 !< Optional variable (diel. permittivity)
     integer :: i_lsf        = -1 !< Optional variable for level set function
     integer :: i_bval       = -1 !< Optional variable for boundary value

     integer :: n_cycle_down = -1 !< Number of relaxation cycles in downward sweep
     integer :: n_cycle_up   = -1 !< Number of relaxation cycles in upward sweep
     integer :: n_cycle_base = -1 !< Number of relaxation cycles at bottom level

     logical :: initialized  = .false.

     !> Routine to call for filling ghost cells near physical boundaries
     procedure(a$D_subr_bc), pointer, nopass   :: sides_bc => null()

     !> Routine to call for filling ghost cells near refinement boundaries
     procedure(a$D_subr_rb), pointer, nopass   :: sides_rb => null()

     !> Subroutine that performs the (non)linear operator
     procedure(mg$D_box_op), pointer, nopass   :: box_op => null()

     !> Subroutine that performs Gauss-Seidel relaxation on a box
     procedure(mg$D_box_gsrb), pointer, nopass :: box_gsrb => null()

     !> Subroutine that corrects the children of a box
     procedure(mg$D_box_corr), pointer, nopass :: box_corr => null()

     !> Subroutine for restriction
     procedure(mg$D_box_rstr), pointer, nopass :: box_rstr => null()
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

     subroutine mg$D_box_rstr(box_c, box_p, iv, mg)
       import
       type(box$D_t), intent(in)    :: box_c !< Child box to restrict
       type(box$D_t), intent(inout) :: box_p !< Parent box to restrict to
       integer, intent(in)         :: iv    !< Variable to restrict
       type(mg$D_t), intent(in)     :: mg
     end subroutine mg$D_box_rstr
  end interface

  public :: mg$D_init_mg
  public :: mg$D_fas_fmg
  public :: mg$D_fas_vcycle
  public :: mg$D_box_op
  public :: mg$D_box_gsrb
  public :: mg$D_box_corr
  public :: mg$D_box_rstr

  ! Automatic selection of operators
  public :: mg$D_set_box_tag
  public :: mg$D_auto_op
  public :: mg$D_auto_gsrb
  public :: mg$D_auto_corr
  public :: mg$D_auto_rstr

  ! Methods for normal Laplacian
  public :: mg$D_box_lpl
  public :: mg$D_box_gsrb_lpl
  public :: mg$D_box_corr_lpl
  public :: mg$D_box_rstr_lpl
  public :: mg$D_sides_rb

  ! Methods for Laplacian with jump in coefficient between boxes
  public :: mg$D_box_lpld
  public :: mg$D_box_gsrb_lpld
  public :: mg$D_box_corr_lpld

  ! Methods for normal Laplacian with internal boundary conditions (LSF)
  public :: mg$D_box_lpllsf
  public :: mg$D_box_gsrb_lpllsf
  public :: mg$D_box_corr_lpllsf
  public :: mg$D_box_rstr_lpllsf

#if $D == 2
  public :: mg$D_box_clpl
  public :: mg$D_box_gsrb_clpl
  public :: mg$D_box_clpld
  public :: mg$D_box_gsrb_clpld
#endif

contains

  !> Check multigrid options or set them to default
  subroutine mg$D_init_mg(mg)
    type(mg$D_t), intent(inout) :: mg           !< Multigrid options

    if (mg%i_phi < 0)                  stop "mg$D_init_mg: i_phi not set"
    if (mg%i_tmp < 0)                  stop "mg$D_init_mg: i_tmp not set"
    if (mg%i_rhs < 0)                  stop "mg$D_init_mg: i_rhs not set"
    if (mg%i_lsf * mg%i_bval < 0) &
         stop "mg$D_init_mg: you have to set both i_lsf and i_bval"

    if (.not. associated(mg%sides_bc)) stop "mg$D_init_mg: sides_bc not set"

    ! Check whether these are set, otherwise use default
    if (mg%n_cycle_down < 0)           mg%n_cycle_down = 2
    if (mg%n_cycle_up < 0)             mg%n_cycle_up = 2
    if (mg%n_cycle_base < 0)           mg%n_cycle_base = 2

    ! Check whether methods are set, otherwise use default (for laplacian)
    if (.not. associated(mg%sides_rb)) mg%sides_rb => mg$D_sides_rb
    if (.not. associated(mg%box_op))   mg%box_op => mg$D_box_lpl
    if (.not. associated(mg%box_gsrb)) mg%box_gsrb => mg$D_box_gsrb_lpl
    if (.not. associated(mg%box_corr)) mg%box_corr => mg$D_box_corr_lpl
    if (.not. associated(mg%box_rstr)) mg%box_rstr => mg$D_box_rstr_lpl

    mg%initialized = .true.
  end subroutine mg$D_init_mg

  subroutine check_mg(mg)
    type(mg$D_t), intent(in) :: mg           !< Multigrid options
    if (.not. mg%initialized) stop "check_mg: you haven't called mg$D_init_mg"
  end subroutine check_mg

  !> Perform FAS-FMG cycle (full approximation scheme, full multigrid). Note
  !> that this routine needs valid ghost cells (for i_phi) on input, and gives
  !> back valid ghost cells on output
  subroutine mg$D_fas_fmg(tree, mg, set_residual, have_guess)
    use m_a$D_utils, only: a$D_boxes_copy_cc
    type(a$D_t), intent(inout)       :: tree !< Tree to do multigrid on
    type(mg$D_t), intent(in)         :: mg   !< Multigrid options
    logical, intent(in)             :: set_residual !< If true, store residual in i_tmp
    logical, intent(in)             :: have_guess   !< If false, start from phi = 0
    integer                         :: lvl, min_lvl

    call check_mg(mg)           ! Check whether mg options are set

    min_lvl = lbound(tree%lvls, 1)

    if (have_guess) then
       do lvl = tree%highest_lvl,  min_lvl+1, -1
          ! Set rhs on coarse grid and restrict phi
          call set_coarse_phi_rhs(tree, lvl, mg)
       end do
    else
       call init_phi_rhs(tree, mg)
    end if

    do lvl = min_lvl, tree%highest_lvl
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

       ! Perform V-cycle, only set residual on last iteration
       call mg$D_fas_vcycle(tree, mg, lvl, &
            set_residual .and. lvl == tree%highest_lvl)
    end do
  end subroutine mg$D_fas_fmg

  !> Perform FAS V-cycle (full approximation scheme). Note that this routine
  !> needs valid ghost cells (for i_phi) on input, and gives back valid ghost
  !> cells on output
  subroutine mg$D_fas_vcycle(tree, mg, highest_lvl, set_residual)
    type(a$D_t), intent(inout) :: tree !< Tree to do multigrid on
    type(mg$D_t), intent(in)   :: mg   !< Multigrid options
    integer, intent(in)        :: highest_lvl !< Maximum level for V-cycle
    logical, intent(in)        :: set_residual !< If true, store residual in i_tmp
    integer                    :: lvl, min_lvl, i, id

    call check_mg(mg)           ! Check whether mg options are set
    min_lvl = lbound(tree%lvls, 1)

    do lvl = highest_lvl,  min_lvl+1, -1
       ! Downwards relaxation
       call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_down)

       ! Set rhs on coarse grid, restrict phi, and copy i_phi to i_tmp for the
       ! correction later
       call update_coarse(tree, lvl, mg)
    end do

    lvl = min_lvl
    call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_base)

    ! Do the upwards part of the v-cycle in the tree
    do lvl = min_lvl+1, highest_lvl
       ! Correct solution at this lvl using lvl-1 data
       ! phi = phi + prolong(phi_coarse - phi_old_coarse)
       call correct_children(tree%boxes, tree%lvls(lvl-1)%parents, mg)

       ! Have to fill ghost cells after correction
       call fill_gc_phi(tree%boxes, tree%lvls(lvl)%ids, mg)

       ! Upwards relaxation
       call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_up)
    end do

    if (set_residual) then
       !$omp parallel private(lvl, i, id)
       do lvl = min_lvl, highest_lvl
          !$omp do
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             call residual_box(tree%boxes(id), mg)
          end do
          !$omd end do nowait
       end do
       !$omp end parallel
    end if
  end subroutine mg$D_fas_vcycle

  !> Fill ghost cells near refinement boundaries which preserves diffusive fluxes.
  !>
  !> Basically, we extrapolate from the fine cells to a corner point, and then
  !> take the average between this corner point and a coarse neighbor to fill
  !> ghost cells for the fine cells.
  subroutine mg$D_sides_rb(boxes, id, nb, iv)
    use m_a$D_ghostcell, only: a$D_gc_prolong0
    type(box$D_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc, ix, dix, i, di, j, dj
#if $D == 3
    integer                     :: k, dk
#endif

    nc = boxes(id)%n_cell

    if (a$D_neighb_low(nb)) then
       ix = 1
       dix = 1
    else
       ix = nc
       dix = -1
    end if

    call a$D_gc_prolong0(boxes, id, nb, iv)

    select case (a$D_neighb_dim(nb))
#if $D == 2
    case (1)
       i = ix
       di = dix
       do j = 1, nc
          dj = -1 + 2 * iand(j, 1)

          ! Bilinear extrapolation (using 4 points)
          boxes(id)%cc(i-di, j, iv) = 0.5_dp * boxes(id)%cc(i-di, j, iv) + &
               1.125_dp * boxes(id)%cc(i, j, iv) - 0.375_dp * &
               (boxes(id)%cc(i+di, j, iv) + boxes(id)%cc(i, j+dj, iv)) &
               + 0.125_dp * boxes(id)%cc(i+di, j+dj, iv)

          ! Extrapolation using 3 points
          ! boxes(id)%cc(i-di, j, iv) = 0.5_dp * boxes(id)%cc(i-di, j, iv) + &
          !      boxes(id)%cc(i, j, iv) - 0.25_dp * &
          !      (boxes(id)%cc(i+di, j, iv) + boxes(id)%cc(i, j+dj, iv))

          ! Extrapolation using 2 points
          ! boxes(id)%cc(i-di, j, iv) = 0.5_dp * boxes(id)%cc(i-di, j, iv) + &
          !      0.75_dp * boxes(id)%cc(i, j, iv) - 0.25_dp * &
          !      boxes(id)%cc(i+di, j+dj, iv)
       end do
    case (2)
       j = ix
       dj = dix
       do i = 1, nc
          di = -1 + 2 * iand(i, 1)

          ! Bilinear extrapolation (using 4 points)
          boxes(id)%cc(i, j-dj, iv) = 0.5_dp * boxes(id)%cc(i, j-dj, iv) + &
               1.125_dp * boxes(id)%cc(i, j, iv) - 0.375_dp * &
               (boxes(id)%cc(i+di, j, iv) + boxes(id)%cc(i, j+dj, iv)) &
               + 0.125_dp * boxes(id)%cc(i+di, j+dj, iv)

          ! Extrapolation using 2 points
          ! boxes(id)%cc(i, j-dj, iv) = 0.5_dp * boxes(id)%cc(i, j-dj, iv) + &
          !      0.75_dp * boxes(id)%cc(i, j, iv) - 0.25_dp * &
          !      boxes(id)%cc(i+di, j+dj, iv)
       end do
#elif $D == 3
    case (1)
       i = ix
       di = dix
       do k = 1, nc
          dk = -1 + 2 * iand(k, 1)
          do j = 1, nc
             dj = -1 + 2 * iand(j, 1)
             ! Trilinear extrapolation (using 8 points)
             ! boxes(id)%cc(i-di, j, k, iv) = &
             !      0.5_dp * boxes(id)%cc(i-di, j, k, iv) + 0.0625_dp * (&
             !      27 * boxes(id)%cc(i, j, k, iv) &
             !      - 9 * boxes(id)%cc(i+di, j, k, iv) &
             !      - 9 * boxes(id)%cc(i, j+dj, k, iv) &
             !      - 9 * boxes(id)%cc(i, j, k+dk, iv) &
             !      + 3 * boxes(id)%cc(i+di, j+dj, k, iv) &
             !      + 3 * boxes(id)%cc(i+di, j, k+dk, iv) &
             !      + 3 * boxes(id)%cc(i, j+dj, k+dk, iv) &
             !      - 1 * boxes(id)%cc(i+di, j+dj, k+dk, iv))

             ! Extrapolation using 2 points
             boxes(id)%cc(i-di, j, k, iv) = &
                  0.5_dp * boxes(id)%cc(i-di, j, k, iv) + &
                  0.75_dp * boxes(id)%cc(i, j, k, iv) - &
                  0.25_dp * boxes(id)%cc(i+di, j+dj, k+dk, iv)
          end do
       end do
    case (2)
       j = ix
       dj = dix
       do k = 1, nc
          dk = -1 + 2 * iand(k, 1)
          do i = 1, nc
             di = -1 + 2 * iand(i, 1)

             ! boxes(id)%cc(i, j-dj, k, iv) = &
             !      0.5_dp * boxes(id)%cc(i, j-dj, k, iv) + 0.0625_dp * (&
             !      27 * boxes(id)%cc(i, j, k, iv) &
             !      - 9 * boxes(id)%cc(i+di, j, k, iv) &
             !      - 9 * boxes(id)%cc(i, j+dj, k, iv) &
             !      - 9 * boxes(id)%cc(i, j, k+dk, iv) &
             !      + 3 * boxes(id)%cc(i+di, j+dj, k, iv) &
             !      + 3 * boxes(id)%cc(i+di, j, k+dk, iv) &
             !      + 3 * boxes(id)%cc(i, j+dj, k+dk, iv) &
             !      - 1 * boxes(id)%cc(i+di, j+dj, k+dk, iv))

             boxes(id)%cc(i, j-dj, k, iv) = &
                  0.5_dp * boxes(id)%cc(i, j-dj, k, iv) + &
                  0.75_dp * boxes(id)%cc(i, j, k, iv) - &
                  0.25_dp * boxes(id)%cc(i+di, j+dj, k+dk, iv)
          end do
       end do
    case (3)
       k = ix
       dk = dix
       do j = 1, nc
          dj = -1 + 2 * iand(j, 1)
          do i = 1, nc
             di = -1 + 2 * iand(i, 1)

             ! boxes(id)%cc(i, j, k-dk, iv) = &
             !      0.5_dp * boxes(id)%cc(i, j, k-dk, iv) + 0.0625_dp * (&
             !      27 * boxes(id)%cc(i, j, k, iv) &
             !      - 9 * boxes(id)%cc(i+di, j, k, iv) &
             !      - 9 * boxes(id)%cc(i, j+dj, k, iv) &
             !      - 9 * boxes(id)%cc(i, j, k+dk, iv) &
             !      + 3 * boxes(id)%cc(i+di, j+dj, k, iv) &
             !      + 3 * boxes(id)%cc(i+di, j, k+dk, iv) &
             !      + 3 * boxes(id)%cc(i, j+dj, k+dk, iv) &
             !      - 1 * boxes(id)%cc(i+di, j+dj, k+dk, iv))

             boxes(id)%cc(i, j, k-dk, iv) = &
                  0.5_dp * boxes(id)%cc(i, j, k-dk, iv) + &
                  0.75_dp * boxes(id)%cc(i, j, k, iv) - &
                  0.25_dp * boxes(id)%cc(i+di, j+dj, k+dk, iv)
          end do
       end do
#endif
    end select

  end subroutine mg$D_sides_rb

  subroutine fill_gc_phi(boxes, ids, mg)
    use m_a$D_ghostcell, only: a$D_gc_box
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:)
    type(mg$D_t), intent(in)     :: mg
    integer                     :: i

    !$omp parallel do
    do i = 1, size(ids)
       call a$D_gc_box(boxes, ids(i), mg%i_phi, &
            mg%sides_rb, mg%sides_bc)
    end do
    !$omp end parallel do
  end subroutine fill_gc_phi

  ! Sets phi = phi + prolong(phi_coarse - phi_old_coarse)
  subroutine correct_children(boxes, ids, mg)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:)
    type(mg$D_t), intent(in)     :: mg
    integer                     :: i, id, nc, i_c, c_id

    !$omp parallel do private(id, nc, i_c, c_id)
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
          if (c_id == af_no_box) cycle
          call mg%box_corr(boxes(id), boxes(c_id), mg)
       end do
    end do
    !$omp end parallel do
  end subroutine correct_children

  subroutine gsrb_boxes(boxes, ids, mg, n_cycle)
    use m_a$D_ghostcell, only: a$D_gc_box
    type(box$D_t), intent(inout) :: boxes(:)
    type(mg$D_t), intent(in)     :: mg
    integer, intent(in)         :: ids(:), n_cycle
    integer                     :: n, i

    !$omp parallel private(n, i)
    do n = 1, 2 * n_cycle
       !$omp do
       do i = 1, size(ids)
          call mg%box_gsrb(boxes(ids(i)), n, mg)
       end do
       !$omp end do

       !$omp do
       do i = 1, size(ids)
          call a$D_gc_box(boxes, ids(i), mg%i_phi, &
               mg%sides_rb, mg%sides_bc)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine gsrb_boxes

  ! Set rhs on coarse grid, restrict phi, and copy i_phi to i_tmp for the
  ! correction later
  subroutine update_coarse(tree, lvl, mg)
    use m_a$D_utils, only: a$D_n_cell, a$D_box_add_cc, a$D_box_copy_cc
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)        :: lvl
    type(mg$D_t), intent(in)   :: mg
    integer                    :: i, id, p_id, nc
#if $D == 2
    real(dp), allocatable :: tmp(:,:)
#elif $D == 3
    real(dp), allocatable :: tmp(:,:,:)
#endif

    id = tree%lvls(lvl)%ids(1)
    nc = a$D_n_cell(tree, lvl)
#if $D == 2
    allocate(tmp(1:nc, 1:nc))
#elif $D == 3
    allocate(tmp(1:nc, 1:nc, 1:nc))
#endif

    ! Restrict phi and the residual
    !$omp parallel do private(id, p_id, tmp)
    do i = 1, size(tree%lvls(lvl)%ids)
       id = tree%lvls(lvl)%ids(i)
       p_id = tree%boxes(id)%parent

       ! Copy the data currently in i_tmp, and restore it later (i_tmp holds the
       ! previous state of i_phi)
#if $D == 2
       tmp = tree%boxes(id)%cc(1:nc, 1:nc, mg%i_tmp)
#elif $D == 3
       tmp = tree%boxes(id)%cc(1:nc, 1:nc, 1:nc, mg%i_tmp)
#endif
       call residual_box(tree%boxes(id), mg)
       call mg%box_rstr(tree%boxes(id), tree%boxes(p_id), mg%i_tmp, mg)
       call mg%box_rstr(tree%boxes(id), tree%boxes(p_id), mg%i_phi, mg)
#if $D == 2
       tree%boxes(id)%cc(1:nc, 1:nc, mg%i_tmp) = tmp
#elif $D == 3
       tree%boxes(id)%cc(1:nc, 1:nc, 1:nc, mg%i_tmp) = tmp
#endif
    end do
    !$omp end parallel do

    call fill_gc_phi(tree%boxes, tree%lvls(lvl-1)%ids, mg)

    ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined, and
    ! store current coarse phi in tmp.

    !$omp parallel do private(id)
    do i = 1, size(tree%lvls(lvl-1)%parents)
       id = tree%lvls(lvl-1)%parents(i)

       ! Set rhs = L phi
       call mg%box_op(tree%boxes(id), mg%i_rhs, mg)

       ! Add tmp (the fine grid residual) to rhs
       call a$D_box_add_cc(tree%boxes(id), mg%i_tmp, mg%i_rhs)

       ! Story a copy of phi in tmp
       call a$D_box_copy_cc(tree%boxes(id), mg%i_phi, mg%i_tmp)
    end do
    !$omp end parallel do
  end subroutine update_coarse

  !> This routine performs the same as update_coarse, but it ignores the tmp
  !> variable
  subroutine set_coarse_phi_rhs(tree, lvl, mg)
    use m_a$D_utils, only: a$D_box_add_cc
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)        :: lvl
    type(mg$D_t), intent(in)   :: mg
    integer                    :: i, id, p_id

    ! In case ghost cells are not filled, fill them here to be sure
    if (lvl == tree%highest_lvl) &
         call fill_gc_phi(tree%boxes, tree%lvls(lvl)%ids, mg)

    !$omp parallel do private(id, p_id)
    do i = 1, size(tree%lvls(lvl)%ids)
       id = tree%lvls(lvl)%ids(i)
       p_id = tree%boxes(id)%parent

       call residual_box(tree%boxes(id), mg)
       call mg%box_rstr(tree%boxes(id), tree%boxes(p_id), mg%i_tmp, mg)
       call mg%box_rstr(tree%boxes(id), tree%boxes(p_id), mg%i_phi, mg)
    end do
    !$omp end parallel do

    call fill_gc_phi(tree%boxes, tree%lvls(lvl-1)%ids, mg)

    ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined

    !$omp parallel do private(id)
    do i = 1, size(tree%lvls(lvl-1)%parents)
       id = tree%lvls(lvl-1)%parents(i)
       call mg%box_op(tree%boxes(id), mg%i_rhs, mg)
       call a$D_box_add_cc(tree%boxes(id), mg%i_tmp, mg%i_rhs)
    end do
    !$omp end parallel do
  end subroutine set_coarse_phi_rhs

  !> Set the initial guess for phi and restrict the rhs
  subroutine init_phi_rhs(tree, mg)
    use m_a$D_utils, only: a$D_box_clear_cc
    type(a$D_t), intent(inout) :: tree
    type(mg$D_t), intent(in)   :: mg
    integer                    :: i, id, p_id, lvl, min_lvl

    ! Start from phi = 0 and restrict rhs
    min_lvl = lbound(tree%lvls, 1)

    !$omp parallel private(lvl, i, id, p_id)
    do lvl = tree%highest_lvl,  min_lvl+1, -1
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call a$D_box_clear_cc(tree%boxes(id), mg%i_phi)
          if (lvl > min_lvl) then
             p_id = tree%boxes(id)%parent
             call mg%box_rstr(tree%boxes(id), tree%boxes(p_id), mg%i_rhs, mg)
          end if
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine init_phi_rhs

  subroutine residual_box(box, mg)
    type(box$D_t), intent(inout) :: box
    type(mg$D_t), intent(in)     :: mg
    integer                     :: nc

    call mg%box_op(box, mg%i_tmp, mg)
    nc = box%n_cell
#if $D == 2
    box%cc(1:nc, 1:nc, mg%i_tmp) = box%cc(1:nc, 1:nc, mg%i_rhs) &
         - box%cc(1:nc, 1:nc, mg%i_tmp)
#elif $D == 3
    box%cc(1:nc, 1:nc, 1:nc, mg%i_tmp) = box%cc(1:nc, 1:nc, 1:nc, mg%i_rhs) &
         - box%cc(1:nc, 1:nc, 1:nc, mg%i_tmp)
#endif
  end subroutine residual_box

  !> Based on the box type, apply a Gauss-Seidel relaxation scheme
  subroutine mg$D_auto_gsrb(box, redblack_cntr, mg)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: redblack_cntr
    type(mg$D_t), intent(in)     :: mg

    if (box%tag == af_init_tag) call mg$D_set_box_tag(box, mg)

    select case(box%tag)
#if $D == 2
    case (mg_normal_box)
       if (box%coord_t == af_cyl) then
          call mg$D_box_gsrb_clpl(box, redblack_cntr, mg)
       else
          call mg$D_box_gsrb_lpl(box, redblack_cntr, mg)
       end if
    case (mg_lsf_box)
       call mg$D_box_gsrb_lpllsf(box, redblack_cntr, mg)
    case (mg_veps_box, mg_ceps_box)
       if (box%coord_t == af_cyl) then
          call mg$D_box_gsrb_clpld(box, redblack_cntr, mg)
       else
          call mg$D_box_gsrb_lpld(box, redblack_cntr, mg)
       end if
#elif $D == 3
    case (mg_normal_box)
       call mg$D_box_gsrb_lpl(box, redblack_cntr, mg)
    case (mg_lsf_box)
       call mg$D_box_gsrb_lpllsf(box, redblack_cntr, mg)
    case (mg_veps_box, mg_ceps_box)
       call mg$D_box_gsrb_lpld(box, redblack_cntr, mg)
#endif
    end select
  end subroutine mg$D_auto_gsrb

  !> Based on the box type, apply the approriate Laplace operator
  subroutine mg$D_auto_op(box, i_out, mg)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: i_out
    type(mg$D_t), intent(in)     :: mg

    if (box%tag == af_init_tag) call mg$D_set_box_tag(box, mg)

    select case(box%tag)
#if $D == 2
    case (mg_normal_box)
       if (box%coord_t == af_cyl) then
          call mg$D_box_clpl(box, i_out, mg)
       else
          call mg$D_box_lpl(box, i_out, mg)
       end if
    case (mg_lsf_box)
       call mg$D_box_lpllsf(box, i_out, mg)
    case (mg_veps_box, mg_ceps_box)
       if (box%coord_t == af_cyl) then
          call mg$D_box_clpld(box, i_out, mg)
       else
          call mg$D_box_lpld(box, i_out, mg)
       end if
#elif $D == 3
    case (mg_normal_box)
       call mg$D_box_lpl(box, i_out, mg)
    case (mg_lsf_box)
       call mg$D_box_lpllsf(box, i_out, mg)
    case (mg_veps_box, mg_ceps_box)
       call mg$D_box_lpld(box, i_out, mg)
#endif
    end select
  end subroutine mg$D_auto_op

  !> Based on the box type, apply the approriate Laplace operator
  subroutine mg$D_auto_rstr(box_c, box_p, iv, mg)
    type(box$D_t), intent(in)    :: box_c
    type(box$D_t), intent(inout) :: box_p
    integer, intent(in)         :: iv
    type(mg$D_t), intent(in)     :: mg

    ! We can only restrict after gsrb, so tag should always be set
    if (box_c%tag == af_init_tag) stop "mg$D_auto_rstr: box_c tag not set"

    select case(box_c%tag)
    case (mg_normal_box, mg_veps_box, mg_ceps_box)
       call mg$D_box_rstr_lpl(box_c, box_p, iv, mg)
    case (mg_lsf_box)
       call mg$D_box_rstr_lpllsf(box_c, box_p, iv, mg)
    end select
  end subroutine mg$D_auto_rstr

  !> Based on the box type, correct the solution of the children
  subroutine mg$D_auto_corr(box_p, box_c, mg)
    type(box$D_t), intent(inout) :: box_c
    type(box$D_t), intent(in)    :: box_p
    type(mg$D_t), intent(in)     :: mg

    if (box_c%tag == af_init_tag) call mg$D_set_box_tag(box_c, mg)

    select case(box_c%tag)
    case (mg_normal_box)
       call mg$D_box_corr_lpl(box_p, box_c, mg)
    case (mg_lsf_box)
       call mg$D_box_corr_lpllsf(box_p, box_c, mg)
    case (mg_veps_box, mg_ceps_box)
       call mg$D_box_corr_lpld(box_p, box_c, mg)
    end select
  end subroutine mg$D_auto_corr

  subroutine mg$D_set_box_tag(box, mg)
    type(box$D_t), intent(inout) :: box
    type(mg$D_t), intent(in)     :: mg
    real(dp) :: a, b
    logical :: is_lsf, is_deps, is_eps

    is_lsf = .false.
    is_eps = .false.
    is_deps = .false.

    if (mg%i_lsf /= -1) then
#if $D == 2
       is_lsf = minval(box%cc(:, :, mg%i_lsf)) * &
            maxval(box%cc(:, :, mg%i_lsf)) < 0
#elif $D == 3
       is_lsf = minval(box%cc(:, :, :, mg%i_lsf)) * &
            maxval(box%cc(:, :, :, mg%i_lsf)) < 0
#endif
    end if

    if (mg%i_eps /= -1) then
#if $D == 2
       a = minval(box%cc(:, :, mg%i_eps))
       b = maxval(box%cc(:, :, mg%i_eps))
#elif $D == 3
       a = minval(box%cc(:, :, :, mg%i_eps))
       b = maxval(box%cc(:, :, :, mg%i_eps))
#endif
       is_deps = (b > a)
       if (.not. is_deps) is_eps = (a < 1 .or. a > 1)
    end if

    if (count([is_lsf, is_eps, is_deps]) > 1) &
         stop "mg$D_set_box_tag: Cannot set lsf and eps tag for same box"

    box%tag = mg_normal_box
    if (is_lsf) box%tag = mg_lsf_box
    if (is_eps) box%tag = mg_ceps_box
    if (is_deps) box%tag = mg_veps_box
  end subroutine mg$D_set_box_tag

  subroutine mg$D_box_corr_lpl(box_p, box_c, mg)
    use m_a$D_prolong
    type(box$D_t), intent(inout) :: box_c
    type(box$D_t), intent(in)    :: box_p
    type(mg$D_t), intent(in)     :: mg

    call a$D_prolong2(box_p, box_c, mg%i_tmp, mg%i_phi, add=.true.)
  end subroutine mg$D_box_corr_lpl

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

  !> Restriction of child box (box_c) to its parent (box_p)
  subroutine mg$D_box_rstr_lpl(box_c, box_p, iv, mg)
    use m_a$D_restrict, only: a$D_restrict_box
    type(box$D_t), intent(in)      :: box_c         !< Child box to restrict
    type(box$D_t), intent(inout)   :: box_p         !< Parent box to restrict to
    integer, intent(in)           :: iv            !< Variable to restrict
    type(mg$D_t), intent(in)       :: mg

    call a$D_restrict_box(box_c, box_p, iv)
  end subroutine mg$D_box_rstr_lpl

  subroutine mg$D_box_gsrb_lpld(box, redblack_cntr, mg)
    type(box$D_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    type(mg$D_t), intent(in)     :: mg
    integer                     :: i, i0, j, nc, i_phi, i_eps, i_rhs
    real(dp)                    :: dx2, u(2*$D), a0, a(2*$D), c(2*$D)
#if $D == 3
    integer                     :: k
#endif

    dx2   = box%dr**2
    nc    = box%n_cell
    i_phi = mg%i_phi
    i_eps = mg%i_eps
    i_rhs = mg%i_rhs

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
#if $D == 2
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          a0 = box%cc(i, j, i_eps) ! value of eps at i,j
          u(1:2) = box%cc(i-1:i+1:2, j, i_phi) ! values at neighbors
          a(1:2) = box%cc(i-1:i+1:2, j, i_eps)
          u(3:4) = box%cc(i, j-1:j+1:2, i_phi)
          a(3:4) = box%cc(i, j-1:j+1:2, i_eps)
          c(:) = 2 * a0 * a(:) / (a0 + a(:))

          box%cc(i, j, i_phi) = &
               (sum(c(:) * u(:)) - dx2 * box%cc(i, j, i_rhs)) / sum(c(:))
       end do
    end do
#elif $D == 3
    do k = 1, nc
       do j = 1, nc
          i0 = 2 - iand(ieor(redblack_cntr, k+j), 1)
          do i = i0, nc, 2
             a0 = box%cc(i, j, k, i_eps) ! value of eps at i,j
             u(1:2) = box%cc(i-1:i+1:2, j, k, i_phi) ! values at neighbors
             a(1:2) = box%cc(i-1:i+1:2, j, k, i_eps)
             u(3:4) = box%cc(i, j-1:j+1:2, k, i_phi)
             a(3:4) = box%cc(i, j-1:j+1:2, k, i_eps)
             u(5:6) = box%cc(i, j, k-1:k+1:2, i_phi)
             a(5:6) = box%cc(i, j, k-1:k+1:2, i_eps)
             c(:) = 2 * a0 * a(:) / (a0 + a(:))

             box%cc(i, j, k, i_phi) = &
                  (sum(c(:) * u(:)) - dx2 * box%cc(i, j, k, i_rhs)) / sum(c(:))
          end do
       end do
    end do
#endif
  end subroutine mg$D_box_gsrb_lpld

  !> Perform Laplacian operator on a box
  subroutine mg$D_box_lpld(box, i_out, mg)
    type(box$D_t), intent(inout) :: box   !< Box to operate on
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    type(mg$D_t), intent(in)     :: mg
    integer                     :: i, j, nc, i_phi, i_eps
    real(dp)                    :: inv_dr_sq, a0, u0, u(2*$D), a(2*$D)
#if $D == 3
    integer                     :: k
#endif


    nc        = box%n_cell
    inv_dr_sq = 1 / box%dr**2
    i_phi     = mg%i_phi
    i_eps     = mg%i_eps

#if $D == 2
    do j = 1, nc
       do i = 1, nc
          u0 = box%cc(i, j, i_phi)
          a0 = box%cc(i, j, i_eps)
          u(1:2) = box%cc(i-1:i+1:2, j, i_phi)
          u(3:4) = box%cc(i, j-1:j+1:2, i_phi)
          a(1:2) = box%cc(i-1:i+1:2, j, i_eps)
          a(3:4) = box%cc(i, j-1:j+1:2, i_eps)

          box%cc(i, j, i_out) = inv_dr_sq * 2 * &
               sum(a0*a(:)/(a0 + a(:)) * (u(:) - u0))
       end do
    end do
#elif $D == 3
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             u0 = box%cc(i, j, k, i_phi)
             a0 = box%cc(i, j, k, i_eps)
             u(1:2) = box%cc(i-1:i+1:2, j, k, i_phi)
             u(3:4) = box%cc(i, j-1:j+1:2, k, i_phi)
             u(5:6) = box%cc(i, j, k-1:k+1:2, i_phi)
             a(1:2) = box%cc(i-1:i+1:2, j, k, i_eps)
             a(3:4) = box%cc(i, j-1:j+1:2, k, i_eps)
             a(5:6) = box%cc(i, j, k-1:k+1:2, i_eps)

             box%cc(i, j, k, i_out) = inv_dr_sq * 2 * &
                  sum(a0*a(:)/(a0 + a(:)) * (u(:) - u0))
          end do
       end do
    end do
#endif

  end subroutine mg$D_box_lpld

  subroutine mg$D_box_corr_lpld(box_p, box_c, mg)
    type(box$D_t), intent(inout)  :: box_c
    type(box$D_t), intent(in)     :: box_p
    type(mg$D_t), intent(in)      :: mg
    integer                      :: ix_offset($D), i_phi, i_corr, i_eps
    integer                      :: nc, i, j, i_c1, i_c2, j_c1, j_c2
    real(dp)                     :: u0, u($D), a0, a($D)
#if $D == 3
    integer                      :: k, k_c1, k_c2
    real(dp), parameter          :: third = 1/3.0_dp
#endif


    nc = box_c%n_cell
    ix_offset = a$D_get_child_offset(box_c)
    i_phi = mg%i_phi
    i_corr = mg%i_tmp
    i_eps = mg%i_eps

    ! In these loops, we calculate the closest coarse index (_c1), and the
    ! one-but-closest (_c2). The fine cell lies in between.
#if $D == 2
    do j = 1, nc
       j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

          u0 = box_p%cc(i_c1, j_c1, i_corr)
          a0 = box_p%cc(i_c1, j_c1, i_eps)
          u(1) = box_p%cc(i_c2, j_c1, i_corr)
          u(2) = box_p%cc(i_c1, j_c2, i_corr)
          a(1) = box_p%cc(i_c2, j_c1, i_eps)
          a(2) = box_p%cc(i_c1, j_c2, i_eps)

          box_c%cc(i, j, i_phi) = box_c%cc(i, j, i_phi) + 0.5_dp * &
               sum( (a0*u0 + a(:)*u(:)) / (a0 + a(:)) )
       end do
    end do
#elif $D == 3
    do k = 1, nc
       k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
       k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

             u0 = box_p%cc(i_c1, j_c1, k_c1, i_corr)
             u(1) = box_p%cc(i_c2, j_c1, k_c1, i_corr)
             u(2) = box_p%cc(i_c1, j_c2, k_c1, i_corr)
             u(3) = box_p%cc(i_c1, j_c1, k_c2, i_corr)
             a0 = box_p%cc(i_c1, j_c1, k_c1, i_eps)
             a(1) = box_p%cc(i_c2, j_c1, k_c1, i_eps)
             a(2) = box_p%cc(i_c1, j_c2, k_c1, i_eps)
             a(3) = box_p%cc(i_c1, j_c1, k_c2, i_eps)

             box_c%cc(i, j, k, i_phi) = box_c%cc(i, j, k, i_phi) + third * &
                  sum((a0*u0 + a(:) * (1.5_dp * u(:) - 0.5_dp * u0)) / &
                  (a0 + a(:)))
          end do
       end do
    end do
#endif
  end subroutine mg$D_box_corr_lpld

  ! Below: multigrid operators for internal boundary conditions. A level set
  ! function defines the location of the interface(s).

  !> For a point a, compute value and distance (between 0, 1) of a neighbor b.
  subroutine lsf_dist_val(lsf_val_bval_a, lsf_val_bval_b, dist, val)
    !> Level set function at a, value at a, boundary value at a
    real(dp), intent(in)  :: lsf_val_bval_a(3)
    !> Level set function at b, value at b, boundary value at b
    real(dp), intent(in)  :: lsf_val_bval_b(3)
    !> Distance to neighbor point (value between 0 and 1)
    real(dp), intent(out) :: dist
    !> Value at neighbor point
    real(dp), intent(out) :: val
    real(dp)              :: lsf_a, lsf_b, bval_a, bval_b

    lsf_a = lsf_val_bval_a(1)
    lsf_b = lsf_val_bval_b(1)

    if (lsf_a * lsf_b < 0) then
       ! There is a boundary between the points
       dist = lsf_a / (lsf_a - lsf_b)
       bval_a = lsf_val_bval_a(3)
       bval_b = lsf_val_bval_b(3)

       ! Interpolate between boundary values
       val  = bval_a * (1-dist) + bval_b * dist
    else
       ! Simply use the value at b
       dist = 1
       val  = lsf_val_bval_b(2)
    end if
  end subroutine lsf_dist_val

  subroutine mg$D_box_corr_lpllsf(box_p, box_c, mg)
    type(box$D_t), intent(inout) :: box_c
    type(box$D_t), intent(in)    :: box_p
    type(mg$D_t), intent(in)     :: mg
    integer                     :: i_phi, i_corr, i_lsf, ix_offset($D)
    integer                     :: nc, i, j, i_c1, i_c2, j_c1, j_c2
    real(dp)                    :: v_a(3), v_b(3), val($D+1), dist($D+1), c($D+1)
#if $D == 3
    integer                     :: k, k_c1, k_c2
#endif

    nc        = box_c%n_cell
    ix_offset = a$D_get_child_offset(box_c)
    i_phi     = mg%i_phi
    i_corr    = mg%i_tmp
    i_lsf     = mg%i_lsf

    ! In these loops, we calculate the closest coarse index (_c1), and the
    ! one-but-closest (_c2). The fine cell lies in between.
#if $D == 2
    do j = 1, nc
       j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

          v_a(1:2) = box_c%cc(i, j, [i_lsf, i_corr])
          v_a(3) = 0.0_dp       ! Boundary value for correction is 0
          v_b(3) = 0.0_dp       ! Idem
          v_b(1:2) = box_p%cc(i_c1, j_c1, [i_lsf, i_corr])
          call lsf_dist_val(v_a, v_b, dist(1), val(1))
          v_b(1:2) = box_p%cc(i_c2, j_c1, [i_lsf, i_corr])
          call lsf_dist_val(v_a, v_b, dist(2), val(2))
          v_b(1:2) = box_p%cc(i_c1, j_c2, [i_lsf, i_corr])
          call lsf_dist_val(v_a, v_b, dist(3), val(3))

          ! This expresses general interpolation between 3 points (on the lines
          ! between the fine and the 3 coarse values).
          c(1) = 2 * dist(2) * dist(3)
          c(2) = dist(1) * dist(3)
          c(3) = dist(1) * dist(2)
          box_c%cc(i, j, i_phi) = box_c%cc(i, j, i_phi) + sum(c * val)/sum(c)
       end do
    end do
#elif $D == 3
    do k = 1, nc
       k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
       k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
       do j = 1, nc
          j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
          j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
          do i = 1, nc
             i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
             i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

             v_a(1:2) = box_c%cc(i, j, k, [i_lsf, i_corr])
             v_a(3) = 0.0_dp       ! Boundary value for correction is 0
             v_b(3) = 0.0_dp       ! Idem
             v_b(1:2) = box_p%cc(i_c1, j_c1, k_c1, [i_lsf, i_corr])
             call lsf_dist_val(v_a, v_b, dist(1), val(1))
             v_b(1:2) = box_p%cc(i_c2, j_c1, k_c1, [i_lsf, i_corr])
             call lsf_dist_val(v_a, v_b, dist(2), val(2))
             v_b(1:2) = box_p%cc(i_c1, j_c2, k_c1, [i_lsf, i_corr])
             call lsf_dist_val(v_a, v_b, dist(3), val(3))
             v_b(1:2) = box_p%cc(i_c1, j_c1, k_c2, [i_lsf, i_corr])
             call lsf_dist_val(v_a, v_b, dist(4), val(4))

             ! This expresses general interpolation between 4 points (on the lines
             ! between the fine and the 4 coarse values).
             c(1) = dist(2) * dist(3) * dist(4)
             c(2) = dist(1) * dist(3) * dist(4)
             c(3) = dist(1) * dist(2) * dist(4)
             c(4) = dist(1) * dist(2) * dist(3)
             box_c%cc(i, j, k, i_phi) = box_c%cc(i, j, k, i_phi) + &
                  sum(c * val)/sum(c)
          end do
       end do
    end do
#endif
  end subroutine mg$D_box_corr_lpllsf

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator
  subroutine mg$D_box_gsrb_lpllsf(box, redblack_cntr, mg)
    type(box$D_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    type(mg$D_t), intent(in)     :: mg
    integer                     :: i, i0, j, nc, i_phi, i_rhs, i_lsf, i_bval
    real(dp)                    :: dx2, dd(2*$D), val(2*$D), v_a(3), v_b(3)
#if $D == 3
    integer                     :: k
#endif

    dx2    = box%dr**2
    nc     = box%n_cell
    i_phi  = mg%i_phi
    i_rhs  = mg%i_rhs
    i_lsf  = mg%i_lsf
    i_bval = mg%i_bval

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
#if $D == 2
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          v_a = box%cc(i, j, [i_lsf, i_phi, i_bval])
          v_b = box%cc(i-1, j, [i_lsf, i_phi, i_bval])
          call lsf_dist_val(v_a, v_b, dd(1), val(1))
          v_b = box%cc(i+1, j, [i_lsf, i_phi, i_bval])
          call lsf_dist_val(v_a, v_b, dd(2), val(2))
          v_b = box%cc(i, j-1, [i_lsf, i_phi, i_bval])
          call lsf_dist_val(v_a, v_b, dd(3), val(3))
          v_b = box%cc(i, j+1, [i_lsf, i_phi, i_bval])
          call lsf_dist_val(v_a, v_b, dd(4), val(4))

          ! Solve for generalized Laplacian (see routine mg$D_box_lpllsf)
          box%cc(i, j, i_phi) = 1 / &
               (dd(1) * dd(2) + dd(3) * dd(4)) * ( &
               (dd(2) * val(1) + dd(1) * val(2)) * &
               dd(3) * dd(4) / (dd(1) + dd(2)) + &
               (dd(4) * val(3) + dd(3) * val(4)) * &
               dd(1) * dd(2) / (dd(3) + dd(4)) - &
               0.5_dp * product(dd) * dx2 * box%cc(i, j, i_rhs))
       end do
    end do
#elif $D == 3
    do k = 1, nc
       do j = 1, nc
          i0 = 2 - iand(ieor(redblack_cntr, k+j), 1)
          do i = i0, nc, 2
             v_a = box%cc(i, j, k, [i_lsf, i_phi, i_bval])
             v_b = box%cc(i-1, j, k, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(1), val(1))
             v_b = box%cc(i+1, j, k, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(2), val(2))
             v_b = box%cc(i, j-1, k, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(3), val(3))
             v_b = box%cc(i, j+1, k, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(4), val(4))
             v_b = box%cc(i, j, k-1, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(5), val(5))
             v_b = box%cc(i, j, k+1, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(6), val(6))

             ! Solve for generalized Laplacian (see routine mg$D_box_lpllsf)
             box%cc(i, j, k, i_phi) = 1 / (1/(dd(1)*dd(2)) + &
                  1/(dd(3)*dd(4)) + 1/(dd(5)*dd(6))) * ( &
                  (dd(2) * val(1) + dd(1) * val(2)) / &
                  ((dd(1) + dd(2)) * dd(1) * dd(2)) + &
                  (dd(4) * val(3) + dd(3) * val(4)) / &
                  ((dd(3) + dd(4)) * dd(3) * dd(4)) + &
                  (dd(6) * val(5) + dd(5) * val(6)) / &
                  ((dd(5) + dd(6)) * dd(5) * dd(6)) - &
                  0.5_dp * dx2 * box%cc(i, j, k, i_rhs))

          end do
       end do
    end do
#endif
  end subroutine mg$D_box_gsrb_lpllsf

  !> Perform Laplacian operator on a box
  subroutine mg$D_box_lpllsf(box, i_out, mg)
    type(box$D_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    type(mg$D_t), intent(in)     :: mg
    integer                     :: i, j, nc, i_phi, i_lsf, i_bval
    real(dp)                    :: inv_dr_sq, dd(2*$D), val(2*$D)
    real(dp)                    :: f0, v_a(3), v_b(3)
#if $D == 3
    integer                     :: k
#endif

    nc        = box%n_cell
    inv_dr_sq = 1 / box%dr**2
    i_phi     = mg%i_phi
    i_lsf     = mg%i_lsf
    i_bval    = mg%i_bval

#if $D == 2
    do j = 1, nc
       do i = 1, nc
          v_a = box%cc(i, j, [i_lsf, i_phi, i_bval])
          v_b = box%cc(i-1, j, [i_lsf, i_phi, i_bval])
          call lsf_dist_val(v_a, v_b, dd(1), val(1))
          v_b = box%cc(i+1, j, [i_lsf, i_phi, i_bval])
          call lsf_dist_val(v_a, v_b, dd(2), val(2))
          v_b = box%cc(i, j-1, [i_lsf, i_phi, i_bval])
          call lsf_dist_val(v_a, v_b, dd(3), val(3))
          v_b = box%cc(i, j+1, [i_lsf, i_phi, i_bval])
          call lsf_dist_val(v_a, v_b, dd(4), val(4))

          ! Generalized Laplacian for neighbors at distance dd * dx
          f0 = box%cc(i, j, i_phi)
          box%cc(i, j, i_out) = 2 * inv_dr_sq * ( &
               (dd(2) * val(1) + dd(1) * val(2) - (dd(1)+dd(2)) * f0) / &
               ((dd(1) + dd(2)) * dd(1) * dd(2)) + &
               (dd(4) * val(3) + dd(3) * val(4) - (dd(3)+dd(4)) * f0) / &
               ((dd(3) + dd(4)) * dd(3) * dd(4)))
       end do
    end do
#elif $D == 3
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             v_a = box%cc(i, j, k, [i_lsf, i_phi, i_bval])
             v_b = box%cc(i-1, j, k, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(1), val(1))
             v_b = box%cc(i+1, j, k, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(2), val(2))
             v_b = box%cc(i, j-1, k, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(3), val(3))
             v_b = box%cc(i, j+1, k, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(4), val(4))
             v_b = box%cc(i, j, k-1, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(5), val(5))
             v_b = box%cc(i, j, k+1, [i_lsf, i_phi, i_bval])
             call lsf_dist_val(v_a, v_b, dd(6), val(6))

             ! Generalized Laplacian for neighbors at distance dd * dx
             f0 = box%cc(i, j, k, i_phi)
             box%cc(i, j, k, i_out) = 2 * inv_dr_sq * ( &
                  (dd(2) * val(1) + dd(1) * val(2) - (dd(1)+dd(2)) * f0) / &
                  ((dd(1) + dd(2)) * dd(1) * dd(2)) + &
                  (dd(4) * val(3) + dd(3) * val(4) - (dd(3)+dd(4)) * f0) / &
                  ((dd(3) + dd(4)) * dd(3) * dd(4)) + &
                  (dd(6) * val(5) + dd(5) * val(6) - (dd(5)+dd(6)) * f0) / &
                  ((dd(5) + dd(6)) * dd(5) * dd(6)))
          end do
       end do
    end do
#endif
  end subroutine mg$D_box_lpllsf

  !> Restriction of child box (box_c) to its parent (box_p)
  subroutine mg$D_box_rstr_lpllsf(box_c, box_p, iv, mg)
    type(box$D_t), intent(in)      :: box_c         !< Child box to restrict
    type(box$D_t), intent(inout)   :: box_p         !< Parent box to restrict to
    integer, intent(in)           :: iv            !< Variable to restrict
    type(mg$D_t), intent(in)       :: mg
    integer                       :: i, j, i_f, j_f, i_c, j_c
    integer                       :: hnc, ix_offset($D), n_ch
#if $D == 2
    logical                       :: child_mask(2, 2)
#elif $D == 3
    logical                       :: child_mask(2, 2, 2)
    integer                       :: k, k_f, k_c
#endif

    hnc       = ishft(box_c%n_cell, -1) ! n_cell / 2
    ix_offset = a$D_get_child_offset(box_c)

#if $D == 2
    do j = 1, hnc
       j_c = ix_offset(2) + j
       j_f = 2 * j - 1
       do i = 1, hnc
          i_c = ix_offset(1) + i
          i_f = 2 * i - 1

          child_mask = (box_p%cc(i_c, j_c, mg%i_lsf) * &
               box_c%cc(i_f:i_f+1, j_f:j_f+1, mg%i_lsf) > 0)
          n_ch = count(child_mask)

          if (n_ch < a$D_num_children .and. n_ch > 0) then
             box_p%cc(i_c, j_c, iv) = 1 / n_ch * &
                  sum(box_c%cc(i_f:i_f+1, j_f:j_f+1, iv), mask=child_mask)
          else                  ! Take average of children
             box_p%cc(i_c, j_c, iv) = 0.25_dp * &
                  sum(box_c%cc(i_f:i_f+1, j_f:j_f+1, iv))
          end if
       end do
    end do
#elif $D == 3
    do k = 1, hnc
       k_c = ix_offset(3) + k
       k_f = 2 * k - 1
       do j = 1, hnc
          j_c = ix_offset(2) + j
          j_f = 2 * j - 1
          do i = 1, hnc
             i_c = ix_offset(1) + i
             i_f = 2 * i - 1

             child_mask = (box_p%cc(i_c, j_c, k_c, mg%i_lsf) * &
                  box_c%cc(i_f:i_f+1, j_f:j_f+1, k_f:k_f+1, mg%i_lsf) > 0)
             n_ch = count(child_mask)

             if (n_ch < a$D_num_children .and. n_ch > 0) then
                box_p%cc(i_c, j_c, k_c, iv) = 1 / n_ch * &
                     sum(box_c%cc(i_f:i_f+1, j_f:j_f+1, k_f:k_f+1, iv), &
                     mask=child_mask)
             else                  ! Take average of children
                box_p%cc(i_c, j_c, k_c, iv) = 0.125_dp * &
                     sum(box_c%cc(i_f:i_f+1, j_f:j_f+1, k_f:k_f+1, iv))
             end if
          end do
       end do
    end do
#endif
  end subroutine mg$D_box_rstr_lpllsf

#if $D == 2
  !> Perform Gauss-Seidel relaxation on box for a cylindrical Laplacian operator
  subroutine mg2_box_gsrb_clpl(box, redblack_cntr, mg)
    type(box2_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, i0, j, nc, i_phi, i_rhs, ioff
    real(dp)                    :: dx2, rfac(2)

    dx2   = box%dr**2
    nc    = box%n_cell
    i_phi = mg%i_phi
    i_rhs = mg%i_rhs
    ioff  = (box%ix(1)-1) * nc

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          rfac = [i+ioff-1, i+ioff] / (i+ioff-0.5_dp)
          box%cc(i, j, i_phi) = 0.25_dp * ( &
               rfac(1) * box%cc(i-1, j, i_phi) + &
               rfac(2) * box%cc(i+1, j, i_phi) + &
               box%cc(i, j+1, i_phi) + box%cc(i, j-1, i_phi) - &
               dx2 * box%cc(i, j, i_rhs))
       end do
    end do
  end subroutine mg2_box_gsrb_clpl

  !> Perform cylindrical Laplacian operator on a box
  subroutine mg2_box_clpl(box, i_out, mg)
    type(box2_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, j, nc, i_phi, ioff
    real(dp)                    :: inv_dr_sq, rfac(2)

    nc        = box%n_cell
    inv_dr_sq = 1 / box%dr**2
    i_phi     = mg%i_phi
    ioff      = (box%ix(1)-1) * nc

    do j = 1, nc
       do i = 1, nc
          rfac = [i+ioff-1, i+ioff] / (i+ioff-0.5_dp)
          box%cc(i, j, i_out) = ( &
               rfac(1) * box%cc(i-1, j, i_phi) + &
               rfac(2) * box%cc(i+1, j, i_phi) + &
               box%cc(i, j-1, i_phi) + box%cc(i, j+1, i_phi) - &
               4 * box%cc(i, j, i_phi)) * inv_dr_sq
       end do
    end do
  end subroutine mg2_box_clpl

  !> Perform cylindrical Laplacian operator on a box with varying eps
  subroutine mg2_box_clpld(box, i_out, mg)
    type(box2_t), intent(inout) :: box   !< Box to operate on
    integer, intent(in)         :: i_out !< Index of variable to store Laplacian in
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, j, nc, i_phi, i_eps, ioff
    real(dp)                    :: inv_dr_sq, a0, u0, u(4), a(4), rfac(4)

    nc        = box%n_cell
    inv_dr_sq = 1 / box%dr**2
    i_phi     = mg%i_phi
    i_eps     = mg%i_eps
    ioff      = (box%ix(1)-1) * nc

    do j = 1, nc
       do i = 1, nc
          rfac(1:2) = [i+ioff-1, i+ioff] / (i+ioff-0.5_dp)
          rfac(3:4) = 1
          u0 = box%cc(i, j, i_phi)
          a0 = box%cc(i, j, i_eps)
          u(1:2) = box%cc(i-1:i+1:2, j, i_phi)
          u(3:4) = box%cc(i, j-1:j+1:2, i_phi)
          a(1:2) = box%cc(i-1:i+1:2, j, i_eps)
          a(3:4) = box%cc(i, j-1:j+1:2, i_eps)

          box%cc(i, j, i_out) = inv_dr_sq * 2 * &
               sum(rfac*a0*a(:)/(a0 + a(:)) * (u(:) - u0))
       end do
    end do
  end subroutine mg2_box_clpld

  !> Perform Gauss-Seidel relaxation on box for a cylindrical Laplacian operator
  !> with a changing eps
  subroutine mg2_box_gsrb_clpld(box, redblack_cntr, mg)
    type(box2_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: redblack_cntr !< Iteration counter
    type(mg2_t), intent(in)     :: mg
    integer                     :: i, i0, j, nc, i_phi, i_eps, i_rhs, ioff
    real(dp)                    :: dx2, u(4), a0, a(4), c(4), rfac(4)

    dx2   = box%dr**2
    nc    = box%n_cell
    i_phi = mg%i_phi
    i_eps = mg%i_eps
    i_rhs = mg%i_rhs
    ioff  = (box%ix(1)-1) * nc

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          rfac(1:2) = [i+ioff-1, i+ioff] / (i+ioff-0.5_dp)
          rfac(3:4) = 1
          a0 = box%cc(i, j, i_eps) ! value of eps at i,j
          u(1:2) = box%cc(i-1:i+1:2, j, i_phi) ! values at neighbors
          a(1:2) = box%cc(i-1:i+1:2, j, i_eps)
          u(3:4) = box%cc(i, j-1:j+1:2, i_phi)
          a(3:4) = box%cc(i, j-1:j+1:2, i_eps)
          c(:) = 2 * a0 * a(:) / (a0 + a(:))

          box%cc(i, j, i_phi) = (sum(c(:) * rfac * u(:)) &
               - dx2 * box%cc(i, j, i_rhs)) / sum(c(:) * rfac)
       end do
    end do
  end subroutine mg2_box_gsrb_clpld
#endif

end module m_a$D_multigrid
