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

  ! The mg module supports different multigrid operators, and uses these tags to
  ! identify boxes / operators
  integer, parameter, public :: mg_normal_box = 1 !< Normal box
  integer, parameter, public :: mg_lsf_box = 2    !< Box with an internal boundary
  integer, parameter, public :: mg_ceps_box = 3   !< Box with constant eps /= 1
  integer, parameter, public :: mg_veps_box = 4   !< Box with varying eps (on face)

  !> Type to store multigrid options in
  type, public :: mg$D_t
     integer :: i_phi        = -1 !< Variable holding solution
     integer :: i_rhs        = -1 !< Variable holding right-hand side
     integer :: i_tmp        = -1 !< Internal variable (holding prev. solution)

     integer :: i_eps        = -1 !< Optional variable (diel. permittivity)
     integer :: i_lsf        = -1 !< Optional variable for level set function

     integer :: n_cycle_down = -1 !< Number of relaxation cycles in downward sweep
     integer :: n_cycle_up   = -1 !< Number of relaxation cycles in upward sweep
     integer :: n_cycle_base = -1 !< Number of relaxation cycles at bottom level

     real(dp) :: lsf_bnd_val = 0.0_dp !< Boundary value used for the LSF method

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

     subroutine mg$D_box_rstr(box_c, box_p, iv, i_to)
       import
       type(box$D_t), intent(in)     :: box_c !< Child box to restrict
       type(box$D_t), intent(inout)  :: box_p !< Parent box to restrict to
       integer, intent(in)           :: iv    !< Variable to restrict
       integer, intent(in), optional :: i_to  !< Destination (if /= iv)
     end subroutine mg$D_box_rstr
  end interface

  public :: mg$D_init_mg
  public :: mg$D_fas_fmg
  public :: mg$D_fas_vcycle
  public :: mg$D_box_op
  public :: mg$D_box_gsrb
  public :: mg$D_box_corr

  ! Automatic selection of operators
  public :: mg$D_auto_op
  public :: mg$D_auto_gsrb
  public :: mg$D_auto_corr

  ! Methods for normal Laplacian
  public :: mg$D_box_lpl
  public :: mg$D_box_gsrb_lpl
  public :: mg$D_box_corr_lpl

  ! Methods for Laplacian with jump in coefficient between boxes
  public :: mg$D_box_lpld
  public :: mg$D_box_gsrb_lpld
  public :: mg$D_box_corr_lpld

  ! Methods for normal Laplacian with internal boundary conditions (LSF)
  public :: mg$D_box_lpllsf
  public :: mg$D_box_gsrb_lpllsf
  public :: mg$D_box_corr_lpllsf

contains

  !> Check multigrid options or set them to default
  subroutine mg$D_init_mg(mg)
    type(mg$D_t), intent(inout) :: mg           !< Multigrid options

    if (mg%i_phi < 0)                  stop "mg$D_init_mg: i_phi not set"
    if (mg%i_tmp < 0)                  stop "mg$D_init_mg: i_tmp not set"
    if (mg%i_rhs < 0)                  stop "mg$D_init_mg: i_rhs not set"

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
    if (.not. associated(mg%box_rstr)) mg%box_rstr => a$D_restrict_box

    mg%initialized = .true.
  end subroutine mg$D_init_mg

  subroutine check_mg(mg)
    type(mg$D_t), intent(in) :: mg           !< Multigrid options
    if (.not. mg%initialized) stop "check_mg: you haven't called mg$D_init_mg"
  end subroutine check_mg

  !> Perform FAS-FMG cycle (full approximation scheme, full multigrid). Note
  !> that this routine needs valid ghost cells (for i_phi) on input, and gives
  !> back valid ghost cells on output
  subroutine mg$D_fas_fmg(tree, mg, set_residual)
    type(a$D_t), intent(inout)       :: tree !< Tree to do multigrid on
    type(mg$D_t), intent(in)         :: mg   !< Multigrid options
    logical, intent(in)             :: set_residual !< If true, store residual in i_tmp
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
       call mg$D_fas_vcycle(tree, mg, lvl, &
            set_residual .and. lvl == tree%max_lvl) ! Only set residual on last iteration
    end do
  end subroutine mg$D_fas_fmg

  !> Perform FAS V-cycle (full approximation scheme). Note that this routine
  !> needs valid ghost cells (for i_phi) on input, and gives back valid ghost
  !> cells on output
  subroutine mg$D_fas_vcycle(tree, mg, max_lvl, set_residual)
    type(a$D_t), intent(inout) :: tree !< Tree to do multigrid on
    type(mg$D_t), intent(in)   :: mg   !< Multigrid options
    integer, intent(in)        :: max_lvl !< Maximum level for V-cycle
    logical, intent(in)        :: set_residual !< If true, store residual in i_tmp
    integer                    :: lvl, min_lvl, i, id

    call check_mg(mg)           ! Check whether mg options are set
    min_lvl = lbound(tree%lvls, 1)

    do lvl = max_lvl,  min_lvl+1, -1
       ! Downwards relaxation
       call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_down)

       ! Set rhs on coarse grid, restrict phi, and copy i_phi to i_tmp for the
       ! correction later
       call update_coarse(tree, lvl, mg)
    end do

    lvl = min_lvl
    call gsrb_boxes(tree%boxes, tree%lvls(lvl)%ids, mg, mg%n_cycle_base)

    ! Do the upwards part of the v-cycle in the tree
    do lvl = min_lvl+1, max_lvl
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
       do lvl = min_lvl, max_lvl
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

  subroutine fill_gc_phi(boxes, ids, mg)
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: ids(:)
    type(mg$D_t), intent(in)     :: mg
    integer                     :: i

    !$omp parallel do
    do i = 1, size(ids)
       call a$D_gc_box_sides(boxes, ids(i), mg%i_phi, &
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
          if (c_id == a5_no_box) cycle
          call mg%box_corr(boxes(id), boxes(c_id), mg)
       end do
    end do
    !$omp end parallel do
  end subroutine correct_children

  subroutine gsrb_boxes(boxes, ids, mg, n_cycle)
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
          call a$D_gc_box_sides(boxes, ids(i), mg%i_phi, &
               mg%sides_rb, mg%sides_bc)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine gsrb_boxes

  ! Set rhs on coarse grid, restrict phi, and copy i_phi to i_tmp for the
  ! correction later
  subroutine update_coarse(tree, lvl, mg)
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
       call mg%box_rstr(tree%boxes(id), tree%boxes(p_id), mg%i_tmp)
       call mg%box_rstr(tree%boxes(id), tree%boxes(p_id), mg%i_phi)
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
       call mg%box_op(tree%boxes(id), mg%i_rhs, mg)
       call a$D_box_add_cc(tree%boxes(id), mg%i_tmp, mg%i_rhs)
       call a$D_box_copy_cc(tree%boxes(id), mg%i_phi, mg%i_tmp)
    end do
    !$omp end parallel do
  end subroutine update_coarse

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

    if (box%tag == a5_init_tag) call mg$D_set_box_tag(box, mg)

    select case(box%tag)
    case (mg_normal_box)
       call mg$D_box_gsrb_lpl(box, redblack_cntr, mg)
    case (mg_lsf_box)
       call mg$D_box_gsrb_lpllsf(box, redblack_cntr, mg)
    case (mg_veps_box, mg_ceps_box)
       call mg$D_box_gsrb_lpld(box, redblack_cntr, mg)
    end select
  end subroutine mg$D_auto_gsrb

  !> Based on the box type, apply the approriate Laplace operator
  subroutine mg$D_auto_op(box, i_out, mg)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: i_out
    type(mg$D_t), intent(in)     :: mg

    if (box%tag == a5_init_tag) call mg$D_set_box_tag(box, mg)

    select case(box%tag)
    case (mg_normal_box)
       call mg$D_box_lpl(box, i_out, mg)
    case (mg_lsf_box)
       call mg$D_box_lpllsf(box, i_out, mg)
    case (mg_veps_box, mg_ceps_box)
       call mg$D_box_lpld(box, i_out, mg)
    end select
  end subroutine mg$D_auto_op

  !> Based on the box type, correct the solution of the children
  subroutine mg$D_auto_corr(box_p, box_c, mg)
    type(box$D_t), intent(inout) :: box_c
    type(box$D_t), intent(in)    :: box_p
    type(mg$D_t), intent(in)     :: mg

    ! We can only correct after gsrb, so tag should always be set
    if (box_p%tag == a5_init_tag) stop "mg$D_auto_corr: box tag not set"

    select case(box_p%tag)
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

  subroutine lsf_dist_val(lsf_a, v_b, b_value, dist, val)
    real(dp), intent(in)  :: lsf_a, v_b(2), b_value
    real(dp), intent(out) :: dist, val

    ! Determine whether there is a boundary
    if (lsf_a * v_b(2) < 0) then
       dist = lsf_a / (lsf_a - v_b(2))
       val  = b_value
    else
       dist = 1
       val  = v_b(1)
    end if
  end subroutine lsf_dist_val

  subroutine mg$D_box_corr_lpllsf(box_p, box_c, mg)
    type(box$D_t), intent(inout) :: box_c
    type(box$D_t), intent(in)    :: box_p
    type(mg$D_t), intent(in)     :: mg
    integer                     :: i_phi, i_corr, i_lsf, ix_offset($D)
    integer                     :: nc, i, j, i_c1, i_c2, j_c1, j_c2
    real(dp)                    :: bval, val($D+1), dist($D+1), lsf, c($D+1)
#if $D == 3
    integer :: k, k_c1, k_c2
#endif

    nc        = box_c%n_cell
    ix_offset = a$D_get_child_offset(box_c)
    i_phi     = mg%i_phi
    i_corr    = mg%i_tmp
    i_lsf     = mg%i_lsf
    bval      = 0.0_dp          ! For the correction, boundaries are zero

    ! In these loops, we calculate the closest coarse index (_c1), and the
    ! one-but-closest (_c2). The fine cell lies in between.
#if $D == 2
    do j = 1, nc
       j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

          lsf = box_c%cc(i, j, i_lsf)
          call lsf_dist_val(lsf, box_p%cc(i_c1, j_c1, [i_corr, i_lsf]), &
               bval, dist(1), val(1))
          call lsf_dist_val(lsf, box_p%cc(i_c2, j_c1, [i_corr, i_lsf]), &
               bval, dist(2), val(2))
          call lsf_dist_val(lsf, box_p%cc(i_c1, j_c2, [i_corr, i_lsf]), &
               bval, dist(3), val(3))

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

             lsf = box_c%cc(i, j, k, i_lsf)
             call lsf_dist_val(lsf, box_p%cc(i_c1, j_c1, k_c1, [i_corr, i_lsf]), &
                  bval, dist(1), val(1))
             call lsf_dist_val(lsf, box_p%cc(i_c2, j_c1, k_c1, [i_corr, i_lsf]), &
                  bval, dist(2), val(2))
             call lsf_dist_val(lsf, box_p%cc(i_c1, j_c2, k_c1, [i_corr, i_lsf]), &
                  bval, dist(3), val(3))
             call lsf_dist_val(lsf, box_p%cc(i_c1, j_c1, k_c2, [i_corr, i_lsf]), &
                  bval, dist(4), val(4))

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
    integer                     :: i, i0, j, nc, i_phi, i_rhs, i_lsf
    real(dp)                    :: bval, dx2, dd(2*$D), val(2*$D), lsf
#if $D == 3
    integer                     :: k
#endif

    dx2   = box%dr**2
    nc    = box%n_cell
    i_phi = mg%i_phi
    i_rhs = mg%i_rhs
    i_lsf = mg%i_lsf
    bval  = mg%lsf_bnd_val

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
#if $D == 2
    do j = 1, nc
       i0 = 2 - iand(ieor(redblack_cntr, j), 1)
       do i = i0, nc, 2
          lsf = box%cc(i, j, i_lsf)
          call lsf_dist_val(lsf, box%cc(i-1, j, [i_phi, i_lsf]), &
               bval, dd(1), val(1))
          call lsf_dist_val(lsf, box%cc(i+1, j, [i_phi, i_lsf]), &
               bval, dd(2), val(2))
          call lsf_dist_val(lsf, box%cc(i, j-1, [i_phi, i_lsf]), &
               bval, dd(3), val(3))
          call lsf_dist_val(lsf, box%cc(i, j+1, [i_phi, i_lsf]), &
               bval, dd(4), val(4))

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
             lsf = box%cc(i, j, k, i_lsf)
             call lsf_dist_val(lsf, box%cc(i-1, j, k, [i_phi, i_lsf]), &
                  bval, dd(1), val(1))
             call lsf_dist_val(lsf, box%cc(i+1, j, k, [i_phi, i_lsf]), &
                  bval, dd(2), val(2))
             call lsf_dist_val(lsf, box%cc(i, j-1, k, [i_phi, i_lsf]), &
                  bval, dd(3), val(3))
             call lsf_dist_val(lsf, box%cc(i, j+1, k, [i_phi, i_lsf]), &
                  bval, dd(4), val(4))
             call lsf_dist_val(lsf, box%cc(i, j, k-1, [i_phi, i_lsf]), &
                  bval, dd(5), val(5))
             call lsf_dist_val(lsf, box%cc(i, j, k+1, [i_phi, i_lsf]), &
                  bval, dd(6), val(6))

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
    integer                     :: i, j, nc, i_phi, i_lsf
    real(dp)                    :: inv_dr_sq, dd(2*$D), val(2*$D)
    real(dp)                    :: bval, f0, lsf
#if $D == 3
    integer                     :: k
#endif

    nc        = box%n_cell
    inv_dr_sq = 1 / box%dr**2
    i_phi     = mg%i_phi
    i_lsf     = mg%i_lsf
    bval      = mg%lsf_bnd_val

#if $D == 2
    do j = 1, nc
       do i = 1, nc
          lsf = box%cc(i, j, i_lsf)
          call lsf_dist_val(lsf, box%cc(i-1, j, [i_phi, i_lsf]), &
               bval, dd(1), val(1))
          call lsf_dist_val(lsf, box%cc(i+1, j, [i_phi, i_lsf]), &
               bval, dd(2), val(2))
          call lsf_dist_val(lsf, box%cc(i, j-1, [i_phi, i_lsf]), &
               bval, dd(3), val(3))
          call lsf_dist_val(lsf, box%cc(i, j+1, [i_phi, i_lsf]), &
               bval, dd(4), val(4))

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
             lsf = box%cc(i, j, k, i_lsf)
             call lsf_dist_val(lsf, box%cc(i-1, j, k, [i_phi, i_lsf]), &
                  bval, dd(1), val(1))
             call lsf_dist_val(lsf, box%cc(i+1, j, k, [i_phi, i_lsf]), &
                  bval, dd(2), val(2))
             call lsf_dist_val(lsf, box%cc(i, j-1, k, [i_phi, i_lsf]), &
                  bval, dd(3), val(3))
             call lsf_dist_val(lsf, box%cc(i, j+1, k, [i_phi, i_lsf]), &
                  bval, dd(4), val(4))
             call lsf_dist_val(lsf, box%cc(i, j, k-1, [i_phi, i_lsf]), &
                  bval, dd(5), val(5))
             call lsf_dist_val(lsf, box%cc(i, j, k+1, [i_phi, i_lsf]), &
                  bval, dd(6), val(6))

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

end module m_mg_$Dd
