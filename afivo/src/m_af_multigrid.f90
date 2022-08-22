#include "../src/cpp_macros.h"
!> This module contains the geometric multigrid routines that come with Afivo
!>
!> @todo How to use box tag with different types of operators?
module m_af_multigrid
  use m_af_types
  use m_af_stencil
  use m_af_ghostcell
  use m_coarse_solver

  implicit none
  private

  public :: mg_t

  public :: mg_init
  public :: mg_destroy
  public :: mg_fas_fmg
  public :: mg_fas_vcycle

  ! Methods for level set functions
  public :: mg_lsf_dist_linear
  public :: mg_lsf_dist_gss

  ! Set tag for each box, depending on operator
  public :: mg_set_box_tag
  public :: mg_set_operators_lvl
  public :: mg_set_operators_tree

  ! Methods for normal Laplacian
  public :: mg_box_lpl_gradient
  public :: mg_sides_rb

  ! To compute the gradient of the solution
  public :: mg_compute_phi_gradient
  public :: mg_compute_field_norm
  public :: mg_box_field_norm

contains

  !> Check multigrid options or set them to default
  subroutine mg_init(tree, mg)
    use m_af_core, only: af_set_cc_methods
    type(af_t), intent(inout) :: tree !< Tree to do multigrid on
    type(mg_t), intent(inout) :: mg   !< Multigrid options

    if (mg%i_phi < 0) stop "mg_init: i_phi not set"
    if (mg%i_tmp < 0) stop "mg_init: i_tmp not set"
    if (mg%i_rhs < 0) stop "mg_init: i_rhs not set"

    if (.not. associated(mg%sides_bc)) stop "mg_init: sides_bc not set"
    if (.not. tree%ready) error stop "mg_init: tree not initialized"

    ! Check whether methods are set, otherwise use default (for laplacian)
    if (.not. associated(mg%box_op)) mg%box_op => mg_auto_op
    ! if (.not. associated(mg%box_stencil)) mg%box_stencil => mg_auto_stencil
    if (.not. associated(mg%box_gsrb)) mg%box_gsrb => mg_auto_gsrb
    if (.not. associated(mg%box_corr)) mg%box_corr => mg_auto_corr
    if (.not. associated(mg%box_rstr)) mg%box_rstr => mg_auto_rstr
    if (.not. associated(mg%sides_rb)) mg%sides_rb => mg_auto_rb

    ! By default, store new operator and prolongation stencils
    if (mg%operator_key == af_stencil_none) then
       tree%n_stencil_keys_stored = tree%n_stencil_keys_stored + 1
       mg%operator_key = tree%n_stencil_keys_stored
    end if

    if (mg%prolongation_key == af_stencil_none) then
       tree%n_stencil_keys_stored = tree%n_stencil_keys_stored + 1
       mg%prolongation_key = tree%n_stencil_keys_stored
    end if

    if (mg%i_lsf /= -1) then
       error stop "Set tree%mg_i_lsf instead of mg%i_lsf"
    else if (iand(mg%operator_mask, mg_lsf_box) > 0) then
       mg%i_lsf = tree%mg_i_lsf
    end if

    if (mg%i_eps /= -1) then
       error stop "Set tree%mg_i_eps instead of mg%i_eps"
    else if (iand(mg%operator_mask, mg_veps_box+mg_ceps_box) > 0) then
       mg%i_eps = tree%mg_i_eps
    end if

    if (mg%i_lsf /= -1) then
       if (.not. associated(mg%lsf_dist)) &
            mg%lsf_dist => mg_lsf_dist_linear
       if (associated(mg%lsf_dist, mg_lsf_dist_gss) .and. .not. &
            associated(mg%lsf)) then
          error stop "mg_lsf_dist_gss requires mg%lsf to be set"
       end if
    end if

    call mg_set_operators_lvl(tree, mg, 1)
    call coarse_solver_initialize(tree, mg)

    if (mg%i_lsf /= -1) then
       call check_coarse_representation_lsf(tree, mg)
    end if

    ! Set the proper methods for the phi variable
    call af_set_cc_methods(tree, mg%i_phi, mg%sides_bc, mg%sides_rb)

    mg%initialized = .true.

  end subroutine mg_init

  subroutine mg_destroy(mg)
    type(mg_t), intent(inout) :: mg   !< Multigrid options
    if (.not. mg%initialized) error stop "mg%initialized is false"
    call coarse_solver_destroy(mg%csolver)
  end subroutine mg_destroy

  subroutine use_mg(tree, mg)
    type(af_t), intent(inout)      :: tree
    type(mg_t), intent(in), target :: mg !< Multigrid options

    if (.not. mg%initialized) error stop "mg%initialized is false"
    tree%mg_current_operator_mask = mg%operator_mask

    ! Make sure box tags and operators are set
    call mg_set_operators_tree(tree, mg)
  end subroutine use_mg

  subroutine done_with_mg(tree)
    type(af_t), intent(inout)      :: tree
    tree%mg_current_operator_mask = -1
  end subroutine done_with_mg

  !> Perform FAS-FMG cycle (full approximation scheme, full multigrid). Note
  !> that this routine needs valid ghost cells (for i_phi) on input, and gives
  !> back valid ghost cells on output
  subroutine mg_fas_fmg(tree, mg, set_residual, have_guess)
    use m_af_utils, only: af_boxes_copy_cc
    type(af_t), intent(inout)       :: tree !< Tree to do multigrid on
    type(mg_t), intent(inout)         :: mg   !< Multigrid options
    logical, intent(in)             :: set_residual !< If true, store residual in i_tmp
    logical, intent(in)             :: have_guess   !< If false, start from phi = 0
    integer                         :: lvl

    call use_mg(tree, mg)

    if (have_guess) then
       do lvl = tree%highest_lvl,  2, -1
          ! Set rhs on coarse grid and restrict phi
          call set_coarse_phi_rhs(tree, lvl, mg)
       end do
    else
       call init_phi_rhs(tree, mg)
    end if

    ! Handle level 1 grid
    lvl = 1
    call af_boxes_copy_cc(tree%boxes, tree%lvls(lvl)%ids, &
         mg%i_phi, mg%i_tmp)
    call mg_fas_vcycle(tree, mg, set_residual .and. &
         lvl == tree%highest_lvl, lvl, standalone=.false.)

    do lvl = 2, tree%highest_lvl
       ! Store phi_old in tmp
       call af_boxes_copy_cc(tree%boxes, tree%lvls(lvl)%ids, &
            mg%i_phi, mg%i_tmp)

       ! Correct solution at this lvl using lvl-1 data
       ! phi = phi + prolong(phi_coarse - phi_old_coarse)
       call correct_children(tree%boxes, tree%lvls(lvl-1)%parents, mg)

       ! Update ghost cells
       call af_gc_lvl(tree, lvl, [mg%i_phi])

       call mg_fas_vcycle(tree, mg, set_residual .and. &
            lvl == tree%highest_lvl, lvl, standalone=.false.)
    end do

    call done_with_mg(tree)
  end subroutine mg_fas_fmg

  !> Perform FAS V-cycle (full approximation scheme). Note that this routine
  !> needs valid ghost cells (for i_phi) on input, and gives back valid ghost
  !> cells on output
  subroutine mg_fas_vcycle(tree, mg, set_residual, highest_lvl, standalone)
    use m_af_utils, only: af_tree_sum_cc
    type(af_t), intent(inout)     :: tree         !< Tree to do multigrid on
    type(mg_t), intent(in)        :: mg           !< Multigrid options
    logical, intent(in)           :: set_residual !< If true, store residual in i_tmp
    integer, intent(in), optional :: highest_lvl  !< Maximum level for V-cycle
    logical, intent(in), optional :: standalone   !< False if called by other cycle
    integer                       :: lvl, i, id, max_lvl
    logical                       :: by_itself
    real(dp)                      :: sum_phi, mean_phi

    by_itself = .true.; if (present(standalone)) by_itself = standalone

    if (by_itself) then
       call use_mg(tree, mg)
    end if

    max_lvl = tree%highest_lvl
    if (present(highest_lvl)) max_lvl = highest_lvl

    do lvl = max_lvl, 2, -1
       ! Downwards relaxation
       call gsrb_boxes(tree, lvl, mg, mg_cycle_down)

       ! Set rhs on coarse grid, restrict phi, and copy i_phi to i_tmp for the
       ! correction later
       call update_coarse(tree, lvl, mg)
    end do

    call solve_coarse_grid(tree, mg)

    ! Do the upwards part of the v-cycle in the tree
    do lvl = 2, max_lvl
       ! Correct solution at this lvl using lvl-1 data
       ! phi = phi + prolong(phi_coarse - phi_old_coarse)
       call correct_children(tree%boxes, tree%lvls(lvl-1)%parents, mg)

       ! Have to fill ghost cells after correction
       call af_gc_lvl(tree, lvl, [mg%i_phi])

       ! Upwards relaxation
       call gsrb_boxes(tree, lvl, mg, mg_cycle_up)
    end do

    if (set_residual) then
       !$omp parallel private(lvl, i, id)
       do lvl = 1, max_lvl
          !$omp do
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             call residual_box(tree%boxes(id), mg)
          end do
          !$omp end do nowait
       end do
       !$omp end parallel
    end if

    ! Subtract the mean of phi, for example when doing a fully periodic
    ! simulation. Note that the integrated r.h.s. term should also be zero then.
    if (mg%subtract_mean) then
       call af_tree_sum_cc(tree, mg%i_phi, sum_phi)
       mean_phi = sum_phi / af_total_volume(tree)
       !$omp parallel private(lvl, i, id)
       do lvl = 1, max_lvl
          !$omp do
          do i = 1, size(tree%lvls(lvl)%ids)
             id = tree%lvls(lvl)%ids(i)
             tree%boxes(id)%cc(DTIMES(:), mg%i_phi) = &
                  tree%boxes(id)%cc(DTIMES(:), mg%i_phi) - mean_phi
          end do
          !$omp end do nowait
       end do
       !$omp end parallel
    end if

    if (by_itself) then
       call done_with_mg(tree)
    end if

  end subroutine mg_fas_vcycle

  subroutine solve_coarse_grid(tree, mg)
    type(af_t), intent(inout) :: tree !< Tree to do multigrid on
    type(mg_t), intent(in)    :: mg   !< Multigrid options

    call coarse_solver_set_rhs_phi(tree, mg)
    call coarse_solver(mg%csolver)
    call coarse_solver_get_phi(tree, mg)

    ! Set ghost cells for the new coarse grid solution
    call af_gc_lvl(tree, 1, [mg%i_phi])
  end subroutine solve_coarse_grid

  !> Fill ghost cells near refinement boundaries which preserves diffusive fluxes.
  subroutine mg_sides_rb(boxes, id, nb, iv)
    type(box_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)          :: id       !< Id of box
    integer, intent(in)          :: nb       !< Ghost cell direction
    integer, intent(in)          :: iv       !< Ghost cell variable
    integer                      :: nc, ix, dix, IJK, di, co(NDIM)
    integer                      :: hnc, p_id, p_nb_id
#if NDIM == 1
    real(dp)                     :: tmp
    real(dp)                     :: gc
#elif NDIM == 2
    integer                      :: dj
    real(dp)                     :: tmp(0:boxes(id)%n_cell/2+1)
    real(dp)                     :: gc(boxes(id)%n_cell)
#elif NDIM == 3
    integer                      :: dj, dk
    real(dp)                     :: tmp(0:boxes(id)%n_cell/2+1, 0:boxes(id)%n_cell/2+1)
    real(dp)                     :: gc(boxes(id)%n_cell, boxes(id)%n_cell)
#endif
#if NDIM > 1
    real(dp)                     :: grad(NDIM-1)
#endif

    nc = boxes(id)%n_cell
    hnc = nc/2

    p_id = boxes(id)%parent
    p_nb_id = boxes(p_id)%neighbors(nb)
    co = af_get_child_offset(boxes(id))

    associate(box => boxes(p_nb_id))
      ! First fill a temporary array with data next to the fine grid
      select case (nb)
#if NDIM == 1
      case (af_neighb_lowx)
         tmp = box%cc(nc, iv)
      case (af_neighb_highx)
         tmp = box%cc(1, iv)
#elif NDIM == 2
      case (af_neighb_lowx)
         tmp = box%cc(nc, co(2):co(2)+hnc+1, iv)
      case (af_neighb_highx)
         tmp = box%cc(1, co(2):co(2)+hnc+1, iv)
      case (af_neighb_lowy)
         tmp = box%cc(co(1):co(1)+hnc+1, nc, iv)
      case (af_neighb_highy)
         tmp = box%cc(co(1):co(1)+hnc+1, 1, iv)
#elif NDIM == 3
      case (af_neighb_lowx)
         tmp = box%cc(nc, co(2):co(2)+hnc+1, co(3):co(3)+hnc+1, iv)
      case (af_neighb_highx)
         tmp = box%cc(1, co(2):co(2)+hnc+1, co(3):co(3)+hnc+1, iv)
      case (af_neighb_lowy)
         tmp = box%cc(co(1):co(1)+hnc+1, nc, co(3):co(3)+hnc+1, iv)
      case (af_neighb_highy)
         tmp = box%cc(co(1):co(1)+hnc+1, 1, co(3):co(3)+hnc+1, iv)
      case (af_neighb_lowz)
         tmp = box%cc(co(1):co(1)+hnc+1, co(2):co(2)+hnc+1, nc, iv)
      case (af_neighb_highz)
         tmp = box%cc(co(1):co(1)+hnc+1, co(2):co(2)+hnc+1, 1, iv)
#endif
      case default
         error stop "mg_sides_rb: wrong argument for nb"
      end select
    end associate

    ! Now interpolate the coarse grid data to obtain values 'straight' next to
    ! the fine grid points
#if NDIM == 1
    gc = tmp
#elif NDIM == 2
    do i = 1, hnc
       grad(1) = 0.125_dp * (tmp(i+1) - tmp(i-1))
       gc(2*i-1) = tmp(i) - grad(1)
       gc(2*i) = tmp(i) + grad(1)
    end do
#elif NDIM == 3
    do j = 1, hnc
       do i = 1, hnc
          grad(1)          = 0.125_dp * (tmp(i+1, j) - tmp(i-1, j))
          grad(2)          = 0.125_dp * (tmp(i, j+1) - tmp(i, j-1))
          gc(2*i-1, 2*j-1) = tmp(i, j) - grad(1) - grad(2)
          gc(2*i, 2*j-1)   = tmp(i, j) + grad(1) - grad(2)
          gc(2*i-1, 2*j)   = tmp(i, j) - grad(1) + grad(2)
          gc(2*i, 2*j)     = tmp(i, j) + grad(1) + grad(2)
       end do
    end do
#endif

    if (af_neighb_low(nb)) then
       ix = 1
       dix = 1
    else
       ix = nc
       dix = -1
    end if

    select case (af_neighb_dim(nb))
#if NDIM == 1
    case (1)
       i = ix
       di = dix
       boxes(id)%cc(i-di, iv) = 0.5_dp * gc &
            + 0.75_dp * boxes(id)%cc(i, iv) &
            - 0.25_dp * boxes(id)%cc(i+di, iv)
#elif NDIM == 2
    case (1)
       i = ix
       di = dix
       do j = 1, nc
          dj = -1 + 2 * iand(j, 1)
          boxes(id)%cc(i-di, j, iv) = 0.5_dp * gc(j) &
               + 0.75_dp * boxes(id)%cc(i, j, iv) &
               - 0.25_dp * boxes(id)%cc(i+di, j, iv)
       end do
    case (2)
       j = ix
       dj = dix
       do i = 1, nc
          di = -1 + 2 * iand(i, 1)
          boxes(id)%cc(i, j-dj, iv) = 0.5_dp * gc(i) &
               + 0.75_dp * boxes(id)%cc(i, j, iv) &
               - 0.25_dp * boxes(id)%cc(i, j+dj, iv)
       end do
#elif NDIM == 3
    case (1)
       i = ix
       di = dix
       do k = 1, nc
          dk = -1 + 2 * iand(k, 1)
          do j = 1, nc
             dj = -1 + 2 * iand(j, 1)
             boxes(id)%cc(i-di, j, k, iv) = &
                  0.5_dp * gc(j, k) + &
                  0.75_dp * boxes(id)%cc(i, j, k, iv) - &
                  0.25_dp * boxes(id)%cc(i+di, j, k, iv)
          end do
       end do
    case (2)
       j = ix
       dj = dix
       do k = 1, nc
          dk = -1 + 2 * iand(k, 1)
          do i = 1, nc
             di = -1 + 2 * iand(i, 1)
             boxes(id)%cc(i, j-dj, k, iv) = &
                  0.5_dp * gc(i, k) + &
                  0.75_dp * boxes(id)%cc(i, j, k, iv) - &
                  0.25_dp * boxes(id)%cc(i, j+dj, k, iv)
          end do
       end do
    case (3)
       k = ix
       dk = dix
       do j = 1, nc
          dj = -1 + 2 * iand(j, 1)
          do i = 1, nc
             di = -1 + 2 * iand(i, 1)
             boxes(id)%cc(i, j, k-dk, iv) = &
                  0.5_dp * gc(i, j) + &
                  0.75_dp * boxes(id)%cc(i, j, k, iv) - &
                  0.25_dp * boxes(id)%cc(i, j, k+dk, iv)
          end do
       end do
#endif
    end select

  end subroutine mg_sides_rb

  !> Fill ghost cells near refinement boundaries which preserves diffusive
  !> fluxes. This routine does not do interpolation of coarse grid values.
  !> Basically, we extrapolate from the fine cells to a corner point, and then
  !> take the average between this corner point and a coarse neighbor to fill
  !> ghost cells for the fine cells.
  subroutine mg_sides_rb_extrap(boxes, id, nb, iv)
    type(box_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: id        !< Id of box
    integer, intent(in)         :: nb        !< Ghost cell direction
    integer, intent(in)         :: iv        !< Ghost cell variable
    integer                     :: nc, ix, dix, IJK, di
#if NDIM > 1
    integer                     :: dj
#endif
#if NDIM > 2
    integer                     :: dk
#endif

    nc = boxes(id)%n_cell

    if (af_neighb_low(nb)) then
       ix = 1
       dix = 1
    else
       ix = nc
       dix = -1
    end if

    call af_gc_prolong_copy(boxes, id, nb, iv, 0)

    select case (af_neighb_dim(nb))
#if NDIM == 1
    case (1)
       i = ix
       di = dix
       boxes(id)%cc(i-di, iv) = 0.5_dp * boxes(id)%cc(i-di, iv) &
            + 0.75_dp * boxes(id)%cc(i, iv) &
            - 0.25_dp * boxes(id)%cc(i+di, iv)
#elif NDIM == 2
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
#elif NDIM == 3
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

  end subroutine mg_sides_rb_extrap

  ! Sets phi = phi + prolong(phi_coarse - phi_old_coarse)
  subroutine correct_children(boxes, ids, mg)
    type(box_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)         :: ids(:) !< Operate on these boxes
    type(mg_t), intent(in)     :: mg !< Multigrid options
    integer                     :: i, id, nc, i_c, c_id

    !$omp parallel do private(id, nc, i_c, c_id)
    do i = 1, size(ids)
       id = ids(i)
       nc = boxes(id)%n_cell

       ! Store the correction in i_tmp
       boxes(id)%cc(DTIMES(:), mg%i_tmp) = boxes(id)%cc(DTIMES(:), mg%i_phi) - &
            boxes(id)%cc(DTIMES(:), mg%i_tmp)

       do i_c = 1, af_num_children
          c_id = boxes(id)%children(i_c)
          if (c_id == af_no_box) cycle
          call mg%box_corr(boxes(id), boxes(c_id), mg)
       end do
    end do
    !$omp end parallel do
  end subroutine correct_children

  subroutine gsrb_boxes(tree, lvl, mg, type_cycle)
    type(af_t), intent(inout) :: tree       !< Tree containing full grid
    type(mg_t), intent(in)    :: mg         !< Multigrid options
    integer, intent(in)       :: lvl        !< Operate on this refinement level
    integer, intent(in)       :: type_cycle !< Type of cycle to perform
    integer                   :: n, i, n_cycle
    logical                   :: use_corners

    select case (type_cycle)
    case (mg_cycle_down)
       n_cycle = mg%n_cycle_down
    case (mg_cycle_up)
       n_cycle = mg%n_cycle_up
    case default
       error stop "gsrb_boxes: invalid cycle type"
    end select

    associate (ids => tree%lvls(lvl)%ids)
      !$omp parallel private(n, i)
      do n = 1, 2 * n_cycle
         !$omp do
         do i = 1, size(ids)
            call mg%box_gsrb(tree%boxes(ids(i)), n, mg)
         end do
         !$omp end do

         ! If corner ghost cells are not required, only store them during the
         ! final upward cycle
         use_corners = mg%use_corners .or. &
              (type_cycle /= mg_cycle_down .and. n == 2 * n_cycle)

         !$omp do
         do i = 1, size(ids)
            call af_gc_box(tree, ids(i), [mg%i_phi], use_corners)
         end do
         !$omp end do
      end do
      !$omp end parallel
    end associate
  end subroutine gsrb_boxes

  ! Set rhs on coarse grid, restrict phi, and copy i_phi to i_tmp for the
  ! correction later
  subroutine update_coarse(tree, lvl, mg)
    use m_af_utils, only: af_box_add_cc, af_box_copy_cc
    type(af_t), intent(inout) :: tree !< Tree containing full grid
    integer, intent(in)        :: lvl !< Update coarse values at lvl-1
    type(mg_t), intent(in)   :: mg !< Multigrid options
    integer                    :: i, id, p_id, nc
    real(dp), allocatable :: tmp(DTIMES(:))

    id = tree%lvls(lvl)%ids(1)
    nc = tree%n_cell
    allocate(tmp(DTIMES(1:nc)))

    ! Restrict phi and the residual
    !$omp parallel do private(id, p_id, tmp)
    do i = 1, size(tree%lvls(lvl)%ids)
       id = tree%lvls(lvl)%ids(i)
       p_id = tree%boxes(id)%parent

       ! Copy the data currently in i_tmp, and restore it later (i_tmp holds the
       ! previous state of i_phi)
       tmp = tree%boxes(id)%cc(DTIMES(1:nc), mg%i_tmp)
       call residual_box(tree%boxes(id), mg)
       call mg%box_rstr(tree%boxes(id), tree%boxes(p_id), mg%i_tmp, mg)
       call mg%box_rstr(tree%boxes(id), tree%boxes(p_id), mg%i_phi, mg)
       tree%boxes(id)%cc(DTIMES(1:nc), mg%i_tmp) = tmp
    end do
    !$omp end parallel do

    call af_gc_lvl(tree, lvl-1, [mg%i_phi])

    ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined, and
    ! store current coarse phi in tmp.

    !$omp parallel do private(id)
    do i = 1, size(tree%lvls(lvl-1)%parents)
       id = tree%lvls(lvl-1)%parents(i)

       ! Set rhs = L phi
       call mg%box_op(tree%boxes(id), mg%i_rhs, mg)

       ! Add tmp (the fine grid residual) to rhs
       call af_box_add_cc(tree%boxes(id), mg%i_tmp, mg%i_rhs)

       ! Story a copy of phi in tmp
       call af_box_copy_cc(tree%boxes(id), mg%i_phi, mg%i_tmp)
    end do
    !$omp end parallel do
  end subroutine update_coarse

  !> This routine performs the same as update_coarse, but it ignores the tmp
  !> variable
  subroutine set_coarse_phi_rhs(tree, lvl, mg)
    use m_af_utils, only: af_box_add_cc
    type(af_t), intent(inout) :: tree !< Tree containing full grid
    integer, intent(in)        :: lvl !< Update coarse values at lvl-1
    type(mg_t), intent(in)   :: mg !< Multigrid options
    integer                    :: i, id, p_id

    ! Fill ghost cells here to be sure
    if (lvl == tree%highest_lvl) then
       call af_gc_lvl(tree, lvl, [mg%i_phi])
    end if

    !$omp parallel do private(id, p_id)
    do i = 1, size(tree%lvls(lvl)%ids)
       id = tree%lvls(lvl)%ids(i)
       p_id = tree%boxes(id)%parent

       call residual_box(tree%boxes(id), mg)
       call mg%box_rstr(tree%boxes(id), tree%boxes(p_id), mg%i_tmp, mg)
       call mg%box_rstr(tree%boxes(id), tree%boxes(p_id), mg%i_phi, mg)
    end do
    !$omp end parallel do

    call af_gc_lvl(tree, lvl-1, [mg%i_phi])

    ! Set rhs_c = laplacian(phi_c) + restrict(res) where it is refined

    !$omp parallel do private(id)
    do i = 1, size(tree%lvls(lvl-1)%parents)
       id = tree%lvls(lvl-1)%parents(i)
       call mg%box_op(tree%boxes(id), mg%i_rhs, mg)
       call af_box_add_cc(tree%boxes(id), mg%i_tmp, mg%i_rhs)
    end do
    !$omp end parallel do
  end subroutine set_coarse_phi_rhs

  !> Set the initial guess for phi and restrict the rhs
  subroutine init_phi_rhs(tree, mg)
    use m_af_utils, only: af_box_clear_cc
    type(af_t), intent(inout) :: tree !< Full grid
    type(mg_t), intent(in)    :: mg   !< Multigrid options
    integer                   :: i, id, p_id, lvl

    !$omp parallel private(lvl, i, id, p_id)
    do lvl = tree%highest_lvl,  1+1, -1
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call af_box_clear_cc(tree%boxes(id), mg%i_phi)
          if (lvl > 1) then
             p_id = tree%boxes(id)%parent
             call mg%box_rstr(tree%boxes(id), tree%boxes(p_id), mg%i_rhs, mg)
          end if
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine init_phi_rhs

  subroutine residual_box(box, mg)
    type(box_t), intent(inout) :: box !< Operate on this box
    type(mg_t), intent(in)     :: mg !< Multigrid options
    integer                     :: nc

    call mg%box_op(box, mg%i_tmp, mg)
    nc = box%n_cell
    box%cc(DTIMES(1:nc), mg%i_tmp) = box%cc(DTIMES(1:nc), mg%i_rhs) &
         - box%cc(DTIMES(1:nc), mg%i_tmp)
  end subroutine residual_box

  !> Based on the box type, apply a Gauss-Seidel relaxation scheme
  subroutine mg_auto_gsrb(box, redblack_cntr, mg)
    type(box_t), intent(inout) :: box           !< Box to operate on
    integer, intent(in)        :: redblack_cntr !< Iteration count
    type(mg_t), intent(in)     :: mg            !< Multigrid options

    call af_stencil_gsrb_box(box, mg%operator_key, redblack_cntr, &
         mg%i_phi, mg%i_rhs)
  end subroutine mg_auto_gsrb

  !> Store operator stencil for a box
  subroutine mg_store_operator_stencil(box, mg)
    type(box_t), intent(inout) :: box
    type(mg_t), intent(in)     :: mg
    integer                    :: ix

    call af_stencil_prepare_store(box, mg%operator_key, ix)

    select case (mg%operator_type)
    case (mg_normal_operator)
       call mg_box_lpl_stencil(box, mg, ix)
    case (mg_lsf_operator)
       call mg_box_lsf_stencil(box, mg, ix)
    case (mg_eps_operator)
       call mg_box_lpld_stencil(box, mg, ix)
    case (mg_auto_operator)
       ! Use box tag to set operator
       select case (iand(box%tag, mg%operator_mask))
       case (mg_normal_box)
          call mg_box_lpl_stencil(box, mg, ix)
       case (mg_lsf_box)
          call mg_box_lsf_stencil(box, mg, ix)
       case (mg_veps_box, mg_ceps_box)
          call mg_box_lpld_stencil(box, mg, ix)
       case (mg_veps_box+mg_lsf_box, mg_ceps_box+mg_lsf_box)
          call mg_box_lpld_lsf_stencil(box, mg, ix)
       case default
          error stop "mg_store_operator_stencil: unknown box tag"
       end select
    case default
       error stop "mg_store_operator_stencil: unknown mg%operator_type"
    end select

    call af_stencil_check_box(box, mg%operator_key, ix)
  end subroutine mg_store_operator_stencil

  !> Store prolongation stencil for a box
  subroutine mg_store_prolongation_stencil(tree, id, mg)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: id !< Id of box
    type(mg_t), intent(in)    :: mg
    integer                   :: ix, p_id

    call af_stencil_prepare_store(tree%boxes(id), mg%prolongation_key, ix)

    p_id = tree%boxes(id)%parent
    if (p_id <= af_no_box) error stop "Box does not have parent"

    select case (mg%prolongation_type)
    case (mg_prolong_linear)
       call mg_box_prolong_linear_stencil(tree%boxes(id), &
            tree%boxes(p_id), mg, ix)
    case (mg_prolong_sparse)
       call mg_box_prolong_sparse_stencil(tree%boxes(id), &
            tree%boxes(p_id), mg, ix)
    case (mg_prolong_auto)
       ! Use box tag
       associate (box=>tree%boxes(id), box_p=>tree%boxes(p_id))
         select case (iand(box%tag, mg%operator_mask))
         case (mg_normal_box, mg_ceps_box)
            call mg_box_prolong_linear_stencil(box, box_p, mg, ix)
         case (mg_lsf_box, mg_ceps_box+mg_lsf_box)
            if (mg%lsf_use_custom_prolongation) then
               call mg_box_prolong_lsf_stencil(box, box_p, mg, ix)
            else
               call mg_box_prolong_linear_stencil(box, box_p, mg, ix)
            end if
         case (mg_veps_box, mg_veps_box + mg_lsf_box)
            call mg_box_prolong_eps_stencil(box, box_p, mg, ix)
         case default
            error stop "mg_store_prolongation_stencil: unknown box tag"
         end select
       end associate
    case default
       error stop "mg_store_prolongation_stencil: unknown mg%prolongation_type"
    end select

    call af_stencil_check_box(tree%boxes(id), mg%prolongation_key, ix)
  end subroutine mg_store_prolongation_stencil

  !> Based on the box type, apply the approriate operator
  subroutine mg_auto_op(box, i_out, mg)
    type(box_t), intent(inout) :: box   !< Operate on this box
    integer, intent(in)        :: i_out !< Index of output variable
    type(mg_t), intent(in)     :: mg    !< Multigrid options

    call af_stencil_apply_box(box, mg%operator_key, mg%i_phi, i_out)
  end subroutine mg_auto_op

  !> Restriction operator
  subroutine mg_auto_rstr(box_c, box_p, iv, mg)
    type(box_t), intent(in)    :: box_c !< Child box
    type(box_t), intent(inout) :: box_p !< Parent box
    integer, intent(in)        :: iv    !< Index of variable
    type(mg_t), intent(in)     :: mg    !< Multigrid options

    ! Always apply standard restriction
    call mg_box_rstr_lpl(box_c, box_p, iv, mg)
  end subroutine mg_auto_rstr

  !> Set ghost cells near refinement boundaries
  subroutine mg_auto_rb(boxes, id, nb, iv, op_mask)
    type(box_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)        :: id       !< Id of box
    integer, intent(in)        :: nb       !< Ghost cell direction
    integer, intent(in)        :: iv       !< Ghost cell variable
    integer, intent(in)        :: op_mask  !< Operator mask

    select case (iand(boxes(id)%tag, op_mask))
    case (mg_veps_box)
       ! With a variable coefficients, use local extrapolation for ghost cells
       call mg_sides_rb_extrap(boxes, id, nb, iv)
    case default
       call mg_sides_rb(boxes, id, nb, iv)
    end select
  end subroutine mg_auto_rb

  !> Based on the box type, correct the solution of the children
  subroutine mg_auto_corr(box_p, box_c, mg)
    type(box_t), intent(inout) :: box_c !< Child box
    type(box_t), intent(in)    :: box_p !< Parent box
    type(mg_t), intent(in)     :: mg !< Multigrid options

    call af_stencil_prolong_box(box_p, box_c, mg%prolongation_key, &
         mg%i_tmp, mg%i_phi, .true.)
  end subroutine mg_auto_corr

  !> Check where the level-set function could have a root by computing the
  !> numerical gradient
  subroutine get_possible_lsf_root_mask(box, nc, dmax, mg, root_mask)
    type(box_t), intent(in) :: box  !< Box to operate on
    real(dp), intent(in)    :: dmax !< Maximal distance to consider
    integer, intent(in)     :: nc   !< Box size
    type(mg_t), intent(in)  :: mg   !< Multigrid options
    !> Whether there could be a root
    logical, intent(out)    :: root_mask(DTIMES(nc))
    integer                 :: IJK
    real(dp)                :: gradnorm, rr(NDIM)

    if (.not. associated(mg%lsf)) error stop "mg%lsf not set"

    ! Compute gradient
    do KJI_DO(1, nc)
       rr = af_r_cc(box, [IJK])
       gradnorm = norm2(numerical_gradient(mg%lsf, rr))
       root_mask(IJK) = (abs(box%cc(IJK, mg%i_lsf)) < dmax * gradnorm * &
            mg%lsf_gradient_safety_factor)
    end do; CLOSE_DO
  end subroutine get_possible_lsf_root_mask

  !> Check if the level set function could have zeros, and store distances for
  !> the neighboring cells
  subroutine store_lsf_distance_matrix(box, nc, mg, boundary)
    type(box_t), intent(inout) :: box      !< Box to operate on
    integer, intent(in)        :: nc       !< Box size
    type(mg_t), intent(in)     :: mg       !< Multigrid options
    logical, intent(out)       :: boundary !< Whether a boundary is found
    logical                    :: root_mask(DTIMES(nc))
    real(dp)                   :: dd(2*NDIM), a(NDIM)
    integer                    :: ixs(NDIM, nc**NDIM), IJK, ix, n
    integer                    :: n_mask, m
    real(dp)                   :: v(2*NDIM, nc**NDIM)
#if NDIM > 1
    integer                    :: nb, dim, i_step, n_steps
    real(dp)                   :: dist, gradient(NDIM), dvec(NDIM)
    real(dp)                   :: x(NDIM), step_size(NDIM), min_dr

    min_dr = minval(box%dr)
#endif

    n = 0
    boundary = .false.

    call get_possible_lsf_root_mask(box, nc, norm2(box%dr), mg, root_mask)
    n_mask = count(root_mask)

    if (n_mask > 0) then
       ! Store mask of where potential roots are
       call af_stencil_prepare_store(box, mg_lsf_mask_key, ix)
       box%stencils(ix)%stype = stencil_sparse
       box%stencils(ix)%shape = af_stencil_mask
       allocate(box%stencils(ix)%sparse_ix(NDIM, n_mask))

       m = 0
       do KJI_DO(1, nc)
          if (root_mask(IJK)) then
             m = m + 1
             box%stencils(ix)%sparse_ix(:, m) = [IJK]
          end if
       end do; CLOSE_DO
    else
       return
    end if

    ! Compute distances
    do KJI_DO(1, nc)
       if (root_mask(IJK)) then
          a = af_r_cc(box, [IJK])
#if NDIM == 1
          dd(1) = mg%lsf_dist(a, af_r_cc(box, [i-1]), mg)
          dd(2) = mg%lsf_dist(a, af_r_cc(box, [i+1]), mg)
#elif NDIM == 2
          dd(1) = mg%lsf_dist(a, af_r_cc(box, [i-1, j]), mg)
          dd(2) = mg%lsf_dist(a, af_r_cc(box, [i+1, j]), mg)
          dd(3) = mg%lsf_dist(a, af_r_cc(box, [i, j-1]), mg)
          dd(4) = mg%lsf_dist(a, af_r_cc(box, [i, j+1]), mg)
#elif NDIM == 3
          dd(1) = mg%lsf_dist(a, af_r_cc(box, [i-1, j, k]), mg)
          dd(2) = mg%lsf_dist(a, af_r_cc(box, [i+1, j, k]), mg)
          dd(3) = mg%lsf_dist(a, af_r_cc(box, [i, j-1, k]), mg)
          dd(4) = mg%lsf_dist(a, af_r_cc(box, [i, j+1, k]), mg)
          dd(5) = mg%lsf_dist(a, af_r_cc(box, [i, j, k-1]), mg)
          dd(6) = mg%lsf_dist(a, af_r_cc(box, [i, j, k+1]), mg)
#endif

#if NDIM > 1
          ! If no boundaries are found, search along the gradient of the level
          ! set function to check if there is a boundary nearby
          if (min_dr > mg%lsf_length_scale .and. all(dd >= 1)) then
             n_steps = ceiling(min_dr/mg%lsf_length_scale)
             step_size = sign(mg%lsf_length_scale, box%cc(IJK, mg%i_lsf))
             x = a

             do i_step = 1, n_steps
                gradient = numerical_gradient(mg%lsf, x)
                gradient = gradient/max(norm2(gradient), 1e-50_dp) ! unit vector

                ! Search in direction towards zero
                x = x - gradient * step_size
                if (mg%lsf(x) * box%cc(IJK, mg%i_lsf) <= 0) exit
             end do

             dist = mg%lsf_dist(a, x, mg)

             if (dist < 1) then
                ! Rescale
                dist = dist * norm2(x - a)/min_dr
                dvec = x - a

                ! Select closest direction
                dim = maxloc(abs(dvec), dim=1)
                nb = 2 * dim - 1
                if (dvec(dim) > 0) nb = nb + 1

                dd(nb) = dist
             end if
          end if
#endif
       else
          dd(:) = 1.0_dp
       end if

       ! Only store distances for boundaries
       if (any(dd < 1.0_dp)) then
          n = n + 1
          ixs(:, n) = [IJK]
          v(:, n) = dd
       end if
    end do; CLOSE_DO

    if (n > 0) then
       boundary = .true.
       call af_stencil_prepare_store(box, mg_lsf_distance_key, ix)
       box%stencils(ix)%stype = stencil_sparse
       box%stencils(ix)%shape = af_stencil_246
       box%stencils(ix)%sparse_ix = ixs(:, 1:n)
       box%stencils(ix)%sparse_v = v(:, 1:n)
    end if

  end subroutine store_lsf_distance_matrix

  !> Set tag (box type) for a box in the tree
  subroutine mg_set_box_tag(tree, id, mg)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: id
    type(mg_t), intent(in)    :: mg
    logical                   :: has_lsf
    real(dp)                  :: a, b

    associate (box => tree%boxes(id))
      box%tag = mg_normal_box
      if (tree%mg_i_lsf /= -1) then
         call store_lsf_distance_matrix(box, box%n_cell, mg, has_lsf)
         if (has_lsf) then
            ! Add bits indicating a level set function
            box%tag = box%tag + mg_lsf_box
         end if
      end if

      if (tree%mg_i_eps /= -1) then
         a = minval(box%cc(DTIMES(:), tree%mg_i_eps))
         b = maxval(box%cc(DTIMES(:), tree%mg_i_eps))

         if (b > a) then
            ! Variable coefficient
            box%tag = box%tag + mg_veps_box
         else if (maxval(abs([a, b] - 1)) > 1e-8_dp) then
            ! Constant coefficient (but not equal to one)
            box%tag = box%tag + mg_ceps_box
         end if
      end if
    end associate

  end subroutine mg_set_box_tag

  subroutine mg_set_operators_lvl(tree, mg, lvl)
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(in)    :: mg
    integer, intent(in)       :: lvl
    integer                   :: i, id, n

    !$omp parallel do private(i, id, n)
    do i = 1, size(tree%lvls(lvl)%ids)
       id = tree%lvls(lvl)%ids(i)
       associate (box => tree%boxes(id))

         if (box%tag == af_init_tag) then
            ! We could use the tag of the parent box, but sometimes a part of a
            ! sharp feature can only be detected in the children
            call mg_set_box_tag(tree, id, mg)
         end if

         ! Check if stencils are available
         n = af_stencil_index(box, mg%operator_key)
         if (n == af_stencil_none) then
            call mg_store_operator_stencil(box, mg)
            n = af_stencil_index(box, mg%operator_key)
         end if

         ! Right-hand side correction for internal boundaries
         if (allocated(box%stencils(n)%f)) then
            box%stencils(n)%bc_correction = &
                 box%stencils(n)%f * mg%lsf_boundary_value
         end if

         if (lvl > 1) then
            n = af_stencil_index(box, mg%prolongation_key)
            if (n == af_stencil_none) then
               call mg_store_prolongation_stencil(tree, id, mg)
            end if
         end if
       end associate
    end do
    !$omp end parallel do
  end subroutine mg_set_operators_lvl

  subroutine mg_set_operators_tree(tree, mg)
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(in)    :: mg
    integer                   :: lvl

    do lvl = 1, tree%highest_lvl
       call mg_set_operators_lvl(tree, mg, lvl)
    end do
  end subroutine mg_set_operators_tree

  !> Restriction of child box (box_c) to its parent (box_p)
  subroutine mg_box_rstr_lpl(box_c, box_p, iv, mg)
    use m_af_restrict, only: af_restrict_box
    type(box_t), intent(in)      :: box_c         !< Child box to restrict
    type(box_t), intent(inout)   :: box_p         !< Parent box to restrict to
    integer, intent(in)           :: iv            !< Variable to restrict
    type(mg_t), intent(in)       :: mg !< Multigrid options

    if (iv == mg%i_phi) then
       ! Don't use geometry for restriction, since this is inconsistent with the
       ! filling of ghost cells near refinement boundaries
       call af_restrict_box(box_c, box_p, [iv], use_geometry=.false.)
    else
       ! For the right-hand side, use the geometry
       call af_restrict_box(box_c, box_p, [iv], use_geometry=.true.)
    end if
  end subroutine mg_box_rstr_lpl

  !> Store the matrix stencil for each cell of the box. The order of the stencil
  !> is (i, j), (i-1, j), (i+1, j), (i, j-1), (i, j+1) (e.g., -4, 1, 1, 1, 1)
  subroutine mg_box_lpl_stencil(box, mg, ix)
    type(box_t), intent(inout) :: box
    type(mg_t), intent(in)     :: mg
    integer, intent(in)        :: ix !< Stencil index
    real(dp)                   :: inv_dr2(NDIM)
    integer                    :: n_coeff, idim

    box%stencils(ix)%shape    = af_stencil_357
    box%stencils(ix)%stype    = stencil_constant
    box%stencils(ix)%cylindrical_gradient = (box%coord_t == af_cyl)
    n_coeff                   = af_stencil_sizes(af_stencil_357)
    inv_dr2                   = 1 / box%dr**2

    allocate(box%stencils(ix)%c(n_coeff))
    do idim = 1, NDIM
       box%stencils(ix)%c(2*idim:2*idim+1) = inv_dr2(idim)
    end do
    box%stencils(ix)%c(1) = -sum(box%stencils(ix)%c(2:)) - mg%helmholtz_lambda

  end subroutine mg_box_lpl_stencil

  !> Store linear prolongation stencil for standard Laplacian
  subroutine mg_box_prolong_linear_stencil(box, box_p, mg, ix)
    type(box_t), intent(inout) :: box   !< Current box
    type(box_t), intent(in)    :: box_p !< Parent box
    type(mg_t), intent(in)     :: mg
    integer, intent(in)        :: ix    !< Stencil index
    integer                    :: n_coeff

    box%stencils(ix)%shape = af_stencil_p248
    box%stencils(ix)%stype = stencil_constant
    n_coeff                = af_stencil_sizes(af_stencil_p248)

    allocate(box%stencils(ix)%c(n_coeff))
#if NDIM == 1
    box%stencils(ix)%c(:) = [0.75_dp, 0.25_dp]
#elif NDIM == 2
    box%stencils(ix)%c(:) = [9, 3, 3, 1] / 16.0_dp
#elif NDIM == 3
    box%stencils(ix)%c(:) = [27, 9, 9, 3, 9, 3, 3, 1] / 64.0_dp
#endif
  end subroutine mg_box_prolong_linear_stencil

  !> Store sparse linear prolongation stencil for standard Laplacian
  subroutine mg_box_prolong_sparse_stencil(box, box_p, mg, ix)
    type(box_t), intent(inout) :: box   !< Current box
    type(box_t), intent(in)    :: box_p !< Parent box
    type(mg_t), intent(in)     :: mg
    integer, intent(in)        :: ix    !< Stencil index
    integer                    :: n_coeff

    box%stencils(ix)%shape = af_stencil_p234
    box%stencils(ix)%stype = stencil_constant
    n_coeff                = af_stencil_sizes(af_stencil_p234)

    allocate(box%stencils(ix)%c(n_coeff))
#if NDIM == 1
    box%stencils(ix)%c(:) = [0.75_dp, 0.25_dp]
#elif NDIM == 2
    box%stencils(ix)%c(:) = [0.5_dp, 0.25_dp, 0.25_dp]
#elif NDIM == 3
    box%stencils(ix)%c(:) = [0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp]
#endif
  end subroutine mg_box_prolong_sparse_stencil

  !> Store prolongation stencil for standard Laplacian with variable coefficient
  !> that can jump at cell faces
  subroutine mg_box_prolong_eps_stencil(box, box_p, mg, ix)
    type(box_t), intent(inout) :: box   !< Current box
    type(box_t), intent(in)    :: box_p !< Parent box
    type(mg_t), intent(in)     :: mg
    integer, intent(in)        :: ix    !< Stencil index
    real(dp)                   :: a0, a(NDIM)
    integer                    :: n_coeff, i_eps, nc
    logical                    :: success
    integer                    :: IJK, IJK_(c1)
    integer                    :: IJK_(c2), ix_offset(NDIM)
#if NDIM == 3
    real(dp), parameter        :: third = 1/3.0_dp
#endif

    nc                     = box%n_cell
    box%stencils(ix)%shape = af_stencil_p234
    box%stencils(ix)%stype = stencil_variable
    ix_offset              = af_get_child_offset(box)
    n_coeff                = af_stencil_sizes(af_stencil_p234)
    allocate(box%stencils(ix)%v(n_coeff, DTIMES(nc)))

    i_eps = mg%i_eps

    associate (v => box%stencils(ix)%v)
      ! In these loops, we calculate the closest coarse index (_c1), and the
      ! one-but-closest (_c2). The fine cell lies in between.
#if NDIM == 1
      do i = 1, nc
         i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
         i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

         a0 = box_p%cc(i_c1, i_eps)
         a(1) = box_p%cc(i_c2, i_eps)

         ! Get value of phi at coarse cell faces, and average
         v(:, IJK) = [a0, a(1)] / (a0 + a(1))
      end do
#elif NDIM == 2
      do j = 1, nc
         j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
         j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
         do i = 1, nc
            i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
            i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

            a0 = box_p%cc(i_c1, j_c1, i_eps)
            a(1) = box_p%cc(i_c2, j_c1, i_eps)
            a(2) = box_p%cc(i_c1, j_c2, i_eps)

            ! Get value of phi at coarse cell faces, and average
            v(1, IJK) = 0.5_dp * sum(a0 / (a0 + a(:)))
            v(2:, IJK) = 0.5_dp * a(:) / (a0 + a(:))
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
         k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
         do j = 1, nc
            j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
            j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
            do i = 1, nc
               i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
               i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

               a0 = box_p%cc(i_c1, j_c1, k_c1, i_eps)
               a(1) = box_p%cc(i_c2, j_c1, k_c1, i_eps)
               a(2) = box_p%cc(i_c1, j_c2, k_c1, i_eps)
               a(3) = box_p%cc(i_c1, j_c1, k_c2, i_eps)

               ! Get value of phi at coarse cell faces, and average
               v(1, IJK) = third * sum((a0 - 0.5_dp * a(:))/(a0 + a(:)))
               v(2:, IJK) = 0.5_dp * a(:) / (a0 + a(:))
            end do
         end do
      end do
#endif
    end associate

    call af_stencil_try_constant(box, ix, epsilon(1.0_dp), success)

  end subroutine mg_box_prolong_eps_stencil

  !> Store prolongation stencil for standard Laplacian with level set function
  !> for internal boundaries
  subroutine mg_box_prolong_lsf_stencil(box, box_p, mg, ix)
    type(box_t), intent(inout) :: box   !< Current box
    type(box_t), intent(in)    :: box_p !< Parent box
    type(mg_t), intent(in)     :: mg
    integer, intent(in)        :: ix    !< Stencil index
    real(dp)                   :: dd(NDIM+1), a(NDIM)
    integer                    :: n_coeff, i_lsf, nc, ix_mask
    logical                    :: success
    integer                    :: IJK, IJK_(c1), n, n_mask
    integer                    :: IJK_(c2), ix_offset(NDIM)

    nc                     = box%n_cell
    box%stencils(ix)%shape = af_stencil_p234
    box%stencils(ix)%stype = stencil_variable
    ix_offset              = af_get_child_offset(box)
    n_coeff                = af_stencil_sizes(af_stencil_p234)
    i_lsf                  = mg%i_lsf
    allocate(box%stencils(ix)%v(n_coeff, DTIMES(nc)))

    ix_mask = af_stencil_index(box, mg_lsf_mask_key)
    if (ix_mask == af_stencil_none) error stop "No LSF root mask stored"
    n_mask = size(box%stencils(ix_mask)%sparse_ix, 2)

    associate (v => box%stencils(ix)%v)
      ! In these loops, we calculate the closest coarse index (_c1), and the
      ! one-but-closest (_c2). The fine cell lies in between.
#if NDIM == 1
      v(1, :) = 0.75_dp
      v(2, :) = 0.25_dp

      do n = 1, n_mask
         i = box%stencils(ix_mask)%sparse_ix(1, n)
         i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
         i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

         a = af_r_cc(box, [IJK])
         dd(1) = mg%lsf_dist(a, af_r_cc(box_p, [i_c1]), mg)
         dd(2) = mg%lsf_dist(a, af_r_cc(box_p, [i_c2]), mg)

         v(:, IJK) = [3 * dd(2), dd(1)]
         v(:, IJK) = v(:, IJK) / sum(v(:, IJK))
         where (dd < 1) v(:, IJK) = 0
      end do
#elif NDIM == 2
      v(1, :, :) = 0.5_dp
      v(2, :, :) = 0.25_dp
      v(3, :, :) = 0.25_dp

      do n = 1, n_mask
         i = box%stencils(ix_mask)%sparse_ix(1, n)
         j = box%stencils(ix_mask)%sparse_ix(2, n)

         j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
         j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
         i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
         i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

         a = af_r_cc(box, [IJK])
         dd(1) = mg%lsf_dist(a, af_r_cc(box_p, [i_c1, j_c1]), mg)
         dd(2) = mg%lsf_dist(a, af_r_cc(box_p, [i_c2, j_c1]), mg)
         dd(3) = mg%lsf_dist(a, af_r_cc(box_p, [i_c1, j_c2]), mg)

         v(:, IJK) = [2 * dd(2) * dd(3), dd(1) * dd(3), dd(1) * dd(2)]
         v(:, IJK) = v(:, IJK) / sum(v(:, IJK))
         where (dd < 1) v(:, IJK) = 0
      end do
#elif NDIM == 3
      v(:, :, :, :) = 0.25_dp

      do n = 1, n_mask
         i = box%stencils(ix_mask)%sparse_ix(1, n)
         j = box%stencils(ix_mask)%sparse_ix(2, n)
         k = box%stencils(ix_mask)%sparse_ix(3, n)

         k_c1 = ix_offset(3) + ishft(k+1, -1) ! (k+1)/2
         k_c2 = k_c1 + 1 - 2 * iand(k, 1)     ! even: +1, odd: -1
         j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
         j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
         i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
         i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

         a = af_r_cc(box, [IJK])
         dd(1) = mg%lsf_dist(a, af_r_cc(box_p, [i_c1, j_c1, k_c1]), mg)
         dd(2) = mg%lsf_dist(a, af_r_cc(box_p, [i_c2, j_c1, k_c1]), mg)
         dd(3) = mg%lsf_dist(a, af_r_cc(box_p, [i_c1, j_c2, k_c1]), mg)
         dd(4) = mg%lsf_dist(a, af_r_cc(box_p, [i_c1, j_c1, k_c2]), mg)

         v(:, IJK) = [dd(2) * dd(3) * dd(4), &
              dd(1) * dd(3) * dd(4), dd(1) * dd(2) * dd(4), &
              dd(1) * dd(2) * dd(3)]
         v(:, IJK) = v(:, IJK) / sum(v(:, IJK))
         where (dd < 1) v(:, IJK) = 0
      end do
#endif
    end associate

    call af_stencil_try_constant(box, ix, epsilon(1.0_dp), success)

  end subroutine mg_box_prolong_lsf_stencil

  !> Store the matrix stencil for each cell of the box. The order of the stencil
  !> is (i, j), (i-1, j), (i+1, j), (i, j-1), (i, j+1) (e.g., -4, 1, 1, 1, 1)
  subroutine mg_box_lpld_stencil(box, mg, ix)
    type(box_t), intent(inout) :: box
    type(mg_t), intent(in)     :: mg
    integer, intent(in)        :: ix !< Index of the stencil
    integer                    :: IJK, nc, n_coeff
    real(dp)                   :: idr2(2*NDIM), a0, a(2*NDIM)
    logical                    :: success

    nc               = box%n_cell
    idr2(1:2*NDIM:2) = 1/box%dr**2
    idr2(2:2*NDIM:2) = idr2(1:2*NDIM:2)

    box%stencils(ix)%shape = af_stencil_357
    box%stencils(ix)%stype = stencil_variable
    box%stencils(ix)%cylindrical_gradient = (box%coord_t == af_cyl)
    n_coeff                = af_stencil_sizes(af_stencil_357)

    allocate(box%stencils(ix)%v(n_coeff, DTIMES(nc)))

    associate (cc => box%cc, n => mg%i_phi, i_eps => mg%i_eps)
      do KJI_DO(1, nc)
#if NDIM == 1
         a0 = box%cc(i, i_eps)
         a(1:2) = box%cc(i-1:i+1:2, i_eps)
#elif NDIM == 2
         a0 = box%cc(i, j, i_eps)
         a(1:2) = box%cc(i-1:i+1:2, j, i_eps)
         a(3:4) = box%cc(i, j-1:j+1:2, i_eps)
#elif NDIM == 3
         a0 = box%cc(i, j, k, i_eps)
         a(1:2) = box%cc(i-1:i+1:2, j, k, i_eps)
         a(3:4) = box%cc(i, j-1:j+1:2, k, i_eps)
         a(5:6) = box%cc(i, j, k-1:k+1:2, i_eps)
#endif
         box%stencils(ix)%v(2:, IJK) = idr2 * 2 * a0*a(:)/(a0 + a(:))
         box%stencils(ix)%v(1, IJK) = -sum(box%stencils(ix)%v(2:, IJK))
      end do; CLOSE_DO
    end associate

    call af_stencil_try_constant(box, ix, epsilon(1.0_dp), success)

  end subroutine mg_box_lpld_stencil

  !> Store stencil for a box with variable coefficient and level set function
  subroutine mg_box_lpld_lsf_stencil(box, mg, ix)
    type(box_t), intent(inout) :: box
    type(mg_t), intent(in)     :: mg
    integer, intent(in)        :: ix !< Index of the stencil
    integer                    :: IJK, nc, n_coeff
    real(dp)                   :: idr2(2*NDIM), a0, a(2*NDIM), dr2(NDIM)
    integer                    :: n, m, idim, s_ix(NDIM), ix_dist
    real(dp)                   :: dd(2*NDIM)
    real(dp), allocatable      :: all_distances(:, DTIMES(:))
    logical                    :: success

    nc               = box%n_cell
    idr2(1:2*NDIM:2) = 1/box%dr**2
    idr2(2:2*NDIM:2) = idr2(1:2*NDIM:2)
    dr2              = box%dr**2

    box%stencils(ix)%shape = af_stencil_357
    box%stencils(ix)%stype = stencil_variable
    ! A custom correction for cylindrical geometry could be more accurate. This
    ! would require the derivation of a discretization for cells in which both
    ! epsilon changes and a boundary is present.
    box%stencils(ix)%cylindrical_gradient = (box%coord_t == af_cyl)
    n_coeff                = af_stencil_sizes(af_stencil_357)

    allocate(box%stencils(ix)%v(n_coeff, DTIMES(nc)))
    allocate(box%stencils(ix)%f(DTIMES(nc)))
    allocate(box%stencils(ix)%bc_correction(DTIMES(nc)))
    allocate(all_distances(2*NDIM, DTIMES(nc)))

    box%stencils(ix)%f = 0.0_dp
    ! Distance 1 indicates no boundary
    all_distances = 1.0_dp

    ix_dist = af_stencil_index(box, mg_lsf_distance_key)
    if (ix_dist == af_stencil_none) error stop "No distances stored"

    ! Use sparse storage of boundary distances to update all_distances
    do n = 1, size(box%stencils(ix_dist)%sparse_ix, 2)
       s_ix = box%stencils(ix_dist)%sparse_ix(:, n)
       all_distances(:, DINDEX(s_ix)) = &
            box%stencils(ix_dist)%sparse_v(:, n)
    end do

    associate (cc => box%cc, n => mg%i_phi, i_eps => mg%i_eps)
      do KJI_DO(1, nc)
         dd = all_distances(:, IJK)

#if NDIM == 1
         a0 = box%cc(i, i_eps)
         a(1:2) = box%cc(i-1:i+1:2, i_eps)
#elif NDIM == 2
         a0 = box%cc(i, j, i_eps)
         a(1:2) = box%cc(i-1:i+1:2, j, i_eps)
         a(3:4) = box%cc(i, j-1:j+1:2, i_eps)
#elif NDIM == 3
         a0 = box%cc(i, j, k, i_eps)
         a(1:2) = box%cc(i-1:i+1:2, j, k, i_eps)
         a(3:4) = box%cc(i, j-1:j+1:2, k, i_eps)
         a(5:6) = box%cc(i, j, k-1:k+1:2, i_eps)
#endif
         ! Generalized Laplacian for neighbors at distance dd * dx
         do idim = 1, NDIM
            box%stencils(ix)%v(1+2*idim-1, IJK) = 1 / &
                 (0.5_dp * dr2(idim) * (dd(2*idim-1) + dd(2*idim)) * &
                 dd(2*idim-1))
            box%stencils(ix)%v(1+2*idim, IJK) = 1 / &
                 (0.5_dp * dr2(idim) * (dd(2*idim-1) + dd(2*idim)) * &
                 dd(2*idim))
         end do

         ! Account for permittivity
         !> @todo This is not fully accurate when there is both a change in
         !> epsilon and a boundary for the same cell
         box%stencils(ix)%v(2:, IJK) = box%stencils(ix)%v(2:, IJK) * &
              2 * a0*a(:)/(a0 + a(:))

         box%stencils(ix)%v(1, IJK) = -sum(box%stencils(ix)%v(2:, IJK))

         ! Move internal boundaries to right-hand side
         do m = 1, 2 * NDIM
            if (dd(m) < 1.0_dp) then
               box%stencils(ix)%f(IJK) = box%stencils(ix)%f(IJK) - &
                    box%stencils(ix)%v(m+1, IJK)
               box%stencils(ix)%v(m+1, IJK) = 0.0_dp
            end if
         end do
      end do; CLOSE_DO
    end associate

    call af_stencil_try_constant(box, ix, epsilon(1.0_dp), success)

  end subroutine mg_box_lpld_lsf_stencil

  !> Compute distance to boundary starting at point a going to point b, in
  !> the range from [0, 1], with 1 meaning there is no boundary
  function mg_lsf_dist_linear(a, b, mg) result(dist)
    real(dp), intent(in)   :: a(NDIM) !< Start point
    real(dp), intent(in)   :: b(NDIM) !< End point
    type(mg_t), intent(in) :: mg
    real(dp)               :: dist, lsf_a, lsf_b

    lsf_a = mg%lsf(a)
    lsf_b = mg%lsf(b)

    if (lsf_a * lsf_b < 0) then
       ! There is a boundary between the points
       dist = lsf_a / (lsf_a - lsf_b)
       dist = max(dist, mg%lsf_min_rel_distance)
    else
       dist = 1.0_dp
    end if
  end function mg_lsf_dist_linear

  !> Find root of f in the interval [a, b]. If f(a) and f(b) have different
  !> signs, apply bisection directly. Else, first find the (assumed to be)
  !> unique local minimum/maximum to determine a bracket. Return relative
  !> location of root, or 1 if there is no root.
  function mg_lsf_dist_gss(a, b, mg) result(dist)
    real(dp), intent(in)   :: a(NDIM) !< Start point
    real(dp), intent(in)   :: b(NDIM) !< End point
    type(mg_t), intent(in) :: mg
    real(dp)               :: bracket(NDIM, 2), b_new(NDIM)
    real(dp)               :: dist, r_root(NDIM), lsf_a, lsf_b
    integer, parameter     :: max_iter = 100

    lsf_a = mg%lsf(a)
    lsf_b = mg%lsf(b)
    dist  = 1.0_dp

    if (lsf_a * lsf_b <= 0) then
       r_root = bisection(mg%lsf, a, b, mg%lsf_tol, max_iter)
    else
       ! Determine bracket using Golden section search
       bracket = gss(mg%lsf, a, b, minimization=(lsf_a >= 0), &
            tol=mg%lsf_tol, find_bracket=.true.)

       ! Take one of the endpoints of the bracket
       if (mg%lsf(bracket(:, 1)) * lsf_a <= 0) then
          b_new = bracket(:, 1)
       else
          b_new = bracket(:, 2)
       end if

       if (mg%lsf(b_new) * lsf_a > 0) then
          return                ! No root
       else
          r_root = bisection(mg%lsf, a, b_new, mg%lsf_tol, max_iter)
       end if
    end if

    dist = norm2(r_root - a)/norm2(b-a)
    dist = max(dist, mg%lsf_min_rel_distance)
  end function mg_lsf_dist_gss

  !> Simple bisection
  function bisection(f, in_a, in_b, tol, max_iter) result(c)
    procedure(mg_func_lsf) :: f
    real(dp), intent(in)   :: in_a(NDIM), in_b(NDIM), tol
    integer, intent(in)    :: max_iter
    real(dp)               :: a(NDIM), b(NDIM), c(NDIM), fc
    integer                :: n

    a = in_a
    b = in_b

    do n = 1, max_iter
       c = 0.5 * (a + b)
       fc = f(c)
       if (0.5 * norm2(b-a) < tol .or. abs(fc) <= 0) exit

       if (fc * f(a) >= 0) then
          a = c
       else
          b = c
       end if
    end do
  end function bisection

  !> Golden-section search on a line between a and b. Given a function f with a
  !> single local minimum/maximum in the interval [a,b], gss returns a subset
  !> interval [c,d] that contains the minimum/maximum with d-c <= tol. Adapted
  !> from https://en.wikipedia.org/wiki/Golden-section_search
  function gss(f, in_a, in_b, minimization, tol, find_bracket) result(bracket)
    procedure(mg_func_lsf) :: f !< Function to minimize/maximize
    real(dp), intent(in)   :: in_a(NDIM) !< Start coordinate
    real(dp), intent(in)   :: in_b(NDIM) !< End coordinate
    !> Whether to perform minimization or maximization
    logical, intent(in)    :: minimization
    real(dp), intent(in)   :: tol !< Absolute tolerance
    logical, intent(in) :: find_bracket !< Whether to search for a bracket
    real(dp)               :: bracket(NDIM, 2)
    real(dp)               :: a(NDIM), b(NDIM), c(NDIM), d(NDIM)
    real(dp)               :: h(NDIM), yc, yd, ya
    real(dp)               :: invphi, invphi2
    integer                :: n, k

    invphi  = (sqrt(5.0_dp) - 1) / 2 ! 1 / phi
    invphi2 = (3 - sqrt(5.0_dp)) / 2 ! 1 / phi^2

    a = in_a
    b = in_b
    h = b - a

    if (norm2(h) <= tol) then
       bracket(:, 1) = a
       bracket(:, 2) = b
       return
    end if

    ! Required steps to achieve tolerance
    n = int(ceiling(log(tol / norm2(h)) / log(invphi)))

    c = a + invphi2 * h
    d = a + invphi * h
    ya = f(a)                   ! To search for bracket
    yc = f(c)
    yd = f(d)

    do k = 1, n-1
       if ((yc < yd) .eqv. minimization) then
          b = d
          d = c
          yd = yc
          h = invphi * h
          c = a + invphi2 * h
          yc = f(c)
       else
          a = c
          c = d
          yc = yd
          h = invphi * h
          d = a + invphi * h
          yd = f(d)
       end if

       ! Exit early if we are only searching for a bracket
       if (find_bracket .and. ya * yc <= 0 .and. ya * yd <= 0) exit
    end do

    if ((yc < yd) .eqv. minimization) then
       bracket(:, 1) = a
       bracket(:, 2) = d
    else
       bracket(:, 1) = c
       bracket(:, 2) = b
    end if
  end function gss

  !> Store the matrix stencil for each cell of the box. The order of the stencil
  !> is (i, j), (i-1, j), (i+1, j), (i, j-1), (i, j+1) (e.g., -4, 1, 1, 1, 1)
  subroutine mg_box_lsf_stencil(box, mg, ix)
    type(box_t), intent(inout) :: box
    type(mg_t), intent(in)     :: mg
    integer, intent(in)        :: ix !< Index of stencil
    integer                    :: IJK, n, nc, idim, n_coeff
    integer                    :: s_ix(NDIM), ix_dist
    real(dp)                   :: dd(2*NDIM), dr2(NDIM)
    real(dp), allocatable      :: all_distances(:, DTIMES(:))
#if NDIM == 2
    real(dp)                   :: tmp
#endif

    nc = box%n_cell
    dr2 = box%dr**2

    ix_dist = af_stencil_index(box, mg_lsf_distance_key)
    if (ix_dist == af_stencil_none) error stop "No distances stored"

    ! Use stored distances to construct stencil
    box%stencils(ix)%shape = af_stencil_357
    box%stencils(ix)%stype = stencil_variable
    ! Perform a custom correction in cylindrical coordinates
    box%stencils(ix)%cylindrical_gradient = .false.
    n_coeff                = af_stencil_sizes(af_stencil_357)

    allocate(box%stencils(ix)%v(n_coeff, DTIMES(nc)))
    allocate(box%stencils(ix)%f(DTIMES(nc)))
    allocate(box%stencils(ix)%bc_correction(DTIMES(nc)))
    box%stencils(ix)%f = 0.0_dp

    allocate(all_distances(2*NDIM, DTIMES(nc)))

    ! Distance 1 indicates no boundary
    all_distances = 1.0_dp

    ! Use sparse storage of boundary distances to update all_distances
    do n = 1, size(box%stencils(ix_dist)%sparse_ix, 2)
       s_ix = box%stencils(ix_dist)%sparse_ix(:, n)
       all_distances(:, DINDEX(s_ix)) = &
            box%stencils(ix_dist)%sparse_v(:, n)
    end do

    do KJI_DO(1, nc)
       dd = all_distances(:, IJK)

       ! Generalized Laplacian for neighbors at distance dd * dx
       do idim = 1, NDIM
          box%stencils(ix)%v(1+2*idim-1, IJK) = 1 / &
               (0.5_dp * dr2(idim) * (dd(2*idim-1) + dd(2*idim)) * &
               dd(2*idim-1))
          box%stencils(ix)%v(1+2*idim, IJK) = 1 / &
               (0.5_dp * dr2(idim) * (dd(2*idim-1) + dd(2*idim)) * &
               dd(2*idim))
       end do

#if NDIM == 2
       if (box%coord_t == af_cyl) then
          ! Add correction for cylindrical coordinate system. This is a
          ! discretization of 1/r d/dr phi.
          tmp = 1/(box%dr(1) * (dd(1) + dd(2)) * af_cyl_radius_cc(box, i))
          box%stencils(ix)%v(2, IJK) = box%stencils(ix)%v(2, IJK) - tmp
          box%stencils(ix)%v(3, IJK) = box%stencils(ix)%v(3, IJK) + tmp
       end if
#endif
       box%stencils(ix)%v(1, IJK) = -sum(box%stencils(ix)%v(2:, IJK))

       ! Move internal boundaries to right-hand side
       do n = 1, 2 * NDIM
          if (dd(n) < 1.0_dp) then
             box%stencils(ix)%f(IJK) = box%stencils(ix)%f(IJK) - &
                  box%stencils(ix)%v(n+1, IJK)
             box%stencils(ix)%v(n+1, IJK) = 0.0_dp
          end if
       end do
    end do; CLOSE_DO

  end subroutine mg_box_lsf_stencil

  !> Compute the gradient of the potential and store in face-centered variables
  subroutine mg_compute_phi_gradient(tree, mg, i_fc, fac, i_norm)
    type(af_t), intent(inout)     :: tree
    type(mg_t), intent(in)        :: mg
    integer, intent(in)           :: i_fc !< Face-centered indices
    real(dp), intent(in)          :: fac  !< Multiply with this factor
    !> If present, store norm in this cell-centered variable
    integer, intent(in), optional :: i_norm
    integer                       :: lvl, i, id

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          select case (iand(tree%boxes(id)%tag, mg%operator_mask))
          case (mg_normal_box, mg_ceps_box)
             call mg_box_lpl_gradient(tree, id, mg, i_fc, fac)
          case (mg_lsf_box, mg_ceps_box+mg_lsf_box, mg_veps_box+mg_lsf_box)
             ! @todo is this okay for the case mg_veps_box+mg_lsf_box?
             if (af_has_children(tree%boxes(id))) then
                !> @todo Solution on coarse grid can lead to large gradient due
                !> to inconsistencies with level set function
                call mg_box_lpl_gradient(tree, id, mg, i_fc, fac)
             else
                call mg_box_lpllsf_gradient(tree, id, mg, i_fc, fac)
             end if
          case (mg_veps_box)
             ! Should call surface_correct_field_fc afterwards
             call mg_box_lpl_gradient(tree, id, mg, i_fc, fac)
          case (af_init_tag)
             error stop "box tag not set"
          case default
             error stop "unknown box tag"
          end select

          if (present(i_norm)) then
             call mg_box_field_norm(tree, id, i_fc, i_norm)
          end if
       end do
       !$omp end do nowait
    end do
    !$omp end parallel
  end subroutine mg_compute_phi_gradient

  !> Compute the gradient of the potential
  subroutine mg_box_lpl_gradient(tree, id, mg, i_fc, fac)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: id
    type(mg_t), intent(in)    :: mg
    integer, intent(in)       :: i_fc !< Face-centered indices
    real(dp), intent(in)      :: fac  !< Multiply with this factor
    integer                   :: nc, i_phi, i_eps
    real(dp)                  :: inv_dr(NDIM)

    associate(box => tree%boxes(id), cc => tree%boxes(id)%cc)
      nc     = box%n_cell
      i_phi  = mg%i_phi
      inv_dr = fac / box%dr

#if NDIM == 1
      box%fc(1:nc+1, 1, i_fc) = inv_dr(1) * &
           (cc(1:nc+1, i_phi) - cc(0:nc, i_phi))
#elif NDIM == 2
      box%fc(1:nc+1, 1:nc, 1, i_fc) = inv_dr(1) * &
           (cc(1:nc+1, 1:nc, i_phi) - cc(0:nc, 1:nc, i_phi))
      box%fc(1:nc, 1:nc+1, 2, i_fc) = inv_dr(2) * &
           (cc(1:nc, 1:nc+1, i_phi) - cc(1:nc, 0:nc, i_phi))
#elif NDIM == 3
      box%fc(1:nc+1, 1:nc, 1:nc, 1, i_fc) = inv_dr(1) * &
           (cc(1:nc+1, 1:nc, 1:nc, i_phi) - &
           cc(0:nc, 1:nc, 1:nc, i_phi))
      box%fc(1:nc, 1:nc+1, 1:nc, 2, i_fc) = inv_dr(2) * &
           (cc(1:nc, 1:nc+1, 1:nc, i_phi) - &
           cc(1:nc, 0:nc, 1:nc, i_phi))
      box%fc(1:nc, 1:nc, 1:nc+1, 3, i_fc) = inv_dr(3) * &
           (cc(1:nc, 1:nc, 1:nc+1, i_phi) - &
           cc(1:nc, 1:nc, 0:nc, i_phi))
#endif

      if (iand(box%tag, mg%operator_mask) == mg_veps_box) then
         ! Compute fields at the boundaries of the box, where eps can change
         i_eps = mg%i_eps

#if NDIM == 1
         box%fc(1, 1, i_fc) = 2 * inv_dr(1) * &
              (cc(1, i_phi) - cc(0, i_phi)) * &
              cc(0, i_eps) / &
              (cc(1, i_eps) + cc(0, i_eps))
         box%fc(nc+1, 1, i_fc) = 2 * inv_dr(1) * &
              (cc(nc+1, i_phi) - cc(nc, i_phi)) * &
              cc(nc+1, i_eps) / &
              (cc(nc+1, i_eps) + cc(nc, i_eps))
#elif NDIM == 2
         box%fc(1, 1:nc, 1, i_fc) = 2 * inv_dr(1) * &
              (cc(1, 1:nc, i_phi) - cc(0, 1:nc, i_phi)) * &
              cc(0, 1:nc, i_eps) / &
              (cc(1, 1:nc, i_eps) + cc(0, 1:nc, i_eps))
         box%fc(nc+1, 1:nc, 1, i_fc) = 2 * inv_dr(1) * &
              (cc(nc+1, 1:nc, i_phi) - cc(nc, 1:nc, i_phi)) * &
              cc(nc+1, 1:nc, i_eps) / &
              (cc(nc+1, 1:nc, i_eps) + cc(nc, 1:nc, i_eps))
         box%fc(1:nc, 1, 2, i_fc) = 2 * inv_dr(2) * &
              (cc(1:nc, 1, i_phi) - cc(1:nc, 0, i_phi)) * &
              cc(1:nc, 0, i_eps) / &
              (cc(1:nc, 1, i_eps) + cc(1:nc, 0, i_eps))
         box%fc(1:nc, nc+1, 2, i_fc) = 2 * inv_dr(2) * &
              (cc(1:nc, nc+1, i_phi) - cc(1:nc, nc, i_phi)) * &
              cc(1:nc, nc+1, i_eps) / &
              (cc(1:nc, nc+1, i_eps) + cc(1:nc, nc, i_eps))
#elif NDIM == 3
         box%fc(1, 1:nc, 1:nc, 1, i_fc) = 2 * inv_dr(1) * &
              (cc(1, 1:nc, 1:nc, i_phi) - cc(0, 1:nc, 1:nc, i_phi)) * &
              cc(0, 1:nc, 1:nc, i_eps) / &
              (cc(1, 1:nc, 1:nc, i_eps) + cc(0, 1:nc, 1:nc, i_eps))
         box%fc(nc+1, 1:nc, 1:nc, 1, i_fc) = 2 * inv_dr(1) * &
              (cc(nc+1, 1:nc, 1:nc, i_phi) - cc(nc, 1:nc, 1:nc, i_phi)) * &
              cc(nc+1, 1:nc, 1:nc, i_eps) / &
              (cc(nc+1, 1:nc, 1:nc, i_eps) + cc(nc, 1:nc, 1:nc, i_eps))
         box%fc(1:nc, 1, 1:nc, 2, i_fc) = 2 * inv_dr(2) * &
              (cc(1:nc, 1, 1:nc, i_phi) - cc(1:nc, 0, 1:nc, i_phi)) * &
              cc(1:nc, 0, 1:nc, i_eps) / &
              (cc(1:nc, 1, 1:nc, i_eps) + cc(1:nc, 0, 1:nc, i_eps))
         box%fc(1:nc, nc+1, 1:nc, 2, i_fc) = 2 * inv_dr(2) * &
              (cc(1:nc, nc+1, 1:nc, i_phi) - cc(1:nc, nc, 1:nc, i_phi)) * &
              cc(1:nc, nc+1, 1:nc, i_eps) / &
              (cc(1:nc, nc+1, 1:nc, i_eps) + cc(1:nc, nc, 1:nc, i_eps))
         box%fc(1:nc, 1:nc, 1, 3, i_fc) = 2 * inv_dr(3) * &
              (cc(1:nc, 1:nc, 1, i_phi) - cc(1:nc, 1:nc, 0, i_phi)) * &
              cc(1:nc, 1:nc, 0, i_eps) / &
              (cc(1:nc, 1:nc, 1, i_eps) + cc(1:nc, 1:nc, 0, i_eps))
         box%fc(1:nc, 1:nc, nc+1, 3, i_fc) = 2 * inv_dr(3) * &
              (cc(1:nc, 1:nc, nc+1, i_phi) - cc(1:nc, 1:nc, nc, i_phi)) * &
              cc(1:nc, 1:nc, nc+1, i_eps) / &
              (cc(1:nc, 1:nc, nc+1, i_eps) + cc(1:nc, 1:nc, nc, i_eps))
#endif
      end if
    end associate
  end subroutine mg_box_lpl_gradient

  !> Compute norm of face-centered variable
  subroutine mg_compute_field_norm(tree, i_fc, i_norm)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: i_fc   !< Index of face-centered variable
    integer, intent(in)       :: i_norm !< Index of cell-centered variable
    integer                   :: lvl, i, id

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call mg_box_field_norm(tree, id, i_fc, i_norm)
       end do
       !$omp end do nowait
    end do
    !$omp end parallel
  end subroutine mg_compute_field_norm

  subroutine mg_box_field_norm(tree, id, i_fc, i_norm)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: id
    integer, intent(in)       :: i_fc !< Face-centered indices
    !> Store norm in this cell-centered variable
    integer, intent(in)       :: i_norm
    integer                   :: nc

    associate(box => tree%boxes(id))
      nc = box%n_cell
#if NDIM == 1
      box%cc(1:nc, i_norm) = 0.5_dp * sqrt(&
           (box%fc(1:nc, 1, i_fc) + &
           box%fc(2:nc+1, 1, i_fc))**2)
#elif NDIM == 2
      box%cc(1:nc, 1:nc, i_norm) = 0.5_dp * sqrt(&
           (box%fc(1:nc, 1:nc, 1, i_fc) + &
           box%fc(2:nc+1, 1:nc, 1, i_fc))**2 + &
           (box%fc(1:nc, 1:nc, 2, i_fc) + &
           box%fc(1:nc, 2:nc+1, 2, i_fc))**2)
#elif NDIM == 3
      box%cc(1:nc, 1:nc, 1:nc, i_norm) = 0.5_dp * sqrt(&
           (box%fc(1:nc, 1:nc, 1:nc, 1, i_fc) + &
           box%fc(2:nc+1, 1:nc, 1:nc, 1, i_fc))**2 + &
           (box%fc(1:nc, 1:nc, 1:nc, 2, i_fc) + &
           box%fc(1:nc, 2:nc+1, 1:nc, 2, i_fc))**2 + &
           (box%fc(1:nc, 1:nc, 1:nc, 3, i_fc) + &
           box%fc(1:nc, 1:nc, 2:nc+1, 3, i_fc))**2)
#endif
    end associate
  end subroutine mg_box_field_norm

  !> Compute the gradient of the potential with a level set function and store
  !> in face-centered variables. The gradients are computed from the positive
  !> side of the level set function.
  subroutine mg_box_lpllsf_gradient(tree, id, mg, i_fc, fac)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: id
    type(mg_t), intent(in)    :: mg
    integer, intent(in)       :: i_fc !< Face-centered indices
    real(dp), intent(in)      :: fac  !< Multiply with this factor
    integer                   :: IJK
    integer                   :: n, nc, i_phi, ix_dist
    real(dp)                  :: inv_dr(NDIM), dd(2*NDIM)

    ! Compute regular values, correct part of them below
    call mg_box_lpl_gradient(tree, id, mg, i_fc, fac)

    ix_dist = af_stencil_index(tree%boxes(id), mg_lsf_distance_key)
    if (ix_dist == af_stencil_none) error stop "No distances stored"

    associate(box => tree%boxes(id), cc => tree%boxes(id)%cc)
      nc     = box%n_cell
      i_phi  = mg%i_phi
      inv_dr = fac / box%dr

      ! Use sparse storage of boundary distances
      do n = 1, size(box%stencils(ix_dist)%sparse_ix, 2)
         dd = box%stencils(ix_dist)%sparse_v(:, n)

#if NDIM == 1
         i = box%stencils(ix_dist)%sparse_ix(1, n)

         if (dd(1) < 1 .and. cc(IJK, mg%i_lsf) >= 0) then
            box%fc(i, 1, i_fc) = inv_dr(1) * &
                 (cc(i, i_phi) - mg%lsf_boundary_value) / dd(1)
         end if
         if (dd(2) < 1 .and. cc(IJK, mg%i_lsf) >= 0) then
            box%fc(i+1, 1, i_fc) = inv_dr(1) * &
                 (mg%lsf_boundary_value - cc(i, i_phi)) / dd(2)
         end if
#elif NDIM == 2
         i = box%stencils(ix_dist)%sparse_ix(1, n)
         j = box%stencils(ix_dist)%sparse_ix(2, n)

         if (dd(1) < 1 .and. cc(IJK, mg%i_lsf) >= 0) then
            box%fc(i, j, 1, i_fc) = inv_dr(1) * &
                 (cc(i, j, i_phi) - mg%lsf_boundary_value) / dd(1)
         end if
         if (dd(2) < 1 .and. cc(IJK, mg%i_lsf) >= 0) then
            box%fc(i+1, j, 1, i_fc) = inv_dr(1) * &
                 (mg%lsf_boundary_value - cc(i, j, i_phi)) / dd(2)
         end if
         if (dd(3) < 1 .and. cc(IJK, mg%i_lsf) >= 0) then
            box%fc(i, j, 2, i_fc) = inv_dr(2) * &
                 (cc(i, j, i_phi) - mg%lsf_boundary_value) / dd(3)
         end if
         if (dd(4) < 1 .and. cc(IJK, mg%i_lsf) >= 0) then
            box%fc(i, j+1, 2, i_fc) = inv_dr(2) * &
                 (mg%lsf_boundary_value - cc(i, j, i_phi)) / dd(4)
         end if
#elif NDIM == 3
         i = box%stencils(ix_dist)%sparse_ix(1, n)
         j = box%stencils(ix_dist)%sparse_ix(2, n)
         k = box%stencils(ix_dist)%sparse_ix(3, n)

         if (dd(1) < 1 .and. cc(IJK, mg%i_lsf) >= 0) then
            box%fc(i, j, k, 1, i_fc) = inv_dr(1) * &
                 (cc(IJK, i_phi) - mg%lsf_boundary_value) / dd(1)
         end if
         if (dd(2) < 1 .and. cc(IJK, mg%i_lsf) >= 0) then
            box%fc(i+1, j, k, 1, i_fc) = inv_dr(1) * &
                 (mg%lsf_boundary_value - cc(IJK, i_phi)) / dd(2)
         end if
         if (dd(3) < 1 .and. cc(IJK, mg%i_lsf) >= 0) then
            box%fc(i, j, k, 2, i_fc) = inv_dr(2) * &
                 (cc(IJK, i_phi) - mg%lsf_boundary_value) / dd(3)
         end if
         if (dd(4) < 1 .and. cc(IJK, mg%i_lsf) >= 0) then
            box%fc(i, j+1, k, 2, i_fc) = inv_dr(2) * &
                 (mg%lsf_boundary_value - cc(IJK, i_phi)) / dd(4)
         end if
         if (dd(5) < 1 .and. cc(IJK, mg%i_lsf) >= 0) then
            box%fc(i, j, k, 3, i_fc) = inv_dr(3) * &
                 (cc(IJK, i_phi) - mg%lsf_boundary_value) / dd(5)
         end if
         if (dd(6) < 1 .and. cc(IJK, mg%i_lsf) >= 0) then
            box%fc(i, j, k+1, 3, i_fc) = inv_dr(3) * &
                 (mg%lsf_boundary_value - cc(IJK, i_phi)) / dd(6)
         end if
#endif
      end do
    end associate
  end subroutine mg_box_lpllsf_gradient

  !> This method checks whether the level set function is properly defined on
  !> the coarse grid
  subroutine check_coarse_representation_lsf(tree, mg)
    type(af_t), intent(in) :: tree
    type(mg_t), intent(in) :: mg
    integer                :: i, id, ix

    do i = 1, size(tree%lvls(1)%ids)
       id = tree%lvls(1)%ids(i)
       ix = af_stencil_index(tree%boxes(id), mg_lsf_distance_key)
       if (ix > af_stencil_none) exit
    end do

    if (i == size(tree%lvls(1)%ids)+1) then
       print *, "No roots found for level set function on coarse grid"
       print *, "you should probably use a finer coarse grid"
       error stop "level set function not resolved on coarse grid"
    end if

  end subroutine check_coarse_representation_lsf

  !> Get amplitude of numerical gradient of level set function
  function numerical_gradient(f, r) result(gradient)
    procedure(mg_func_lsf) :: f
    real(dp), intent(in)   :: r(NDIM)
    real(dp), parameter    :: sqrteps      = sqrt(epsilon(1.0_dp))
    real(dp), parameter    :: min_stepsize = epsilon(1.0_dp)
    real(dp)               :: r_eval(NDIM), gradient(NDIM)
    real(dp)               :: stepsize(NDIM), flo, fhi
    integer                :: idim

    stepsize = max(min_stepsize, sqrteps * abs(r))
    r_eval = r

    do idim = 1, NDIM
       ! Sample function at (r - step_idim) and (r + step_idim)
       r_eval(idim) = r(idim) - stepsize(idim)
       flo = f(r_eval)

       r_eval(idim) = r(idim) + stepsize(idim)
       fhi = f(r_eval)

       ! Use central difference scheme
       gradient(idim) = (fhi - flo)/(2 * stepsize(idim))

       ! Reset to original coordinate
       r_eval(idim) = r(idim)
    end do
  end function numerical_gradient

end module m_af_multigrid
