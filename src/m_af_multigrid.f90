#include "../src/cpp_macros.h"
!> This module contains the geometric multigrid routines that come with Afivo
module m_af_multigrid
  use m_af_types
  use m_af_stencil
  use m_mg_types
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

  ! Automatic selection of operators
  public :: mg_set_box_tag

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

    ! Check whether these are set, otherwise use default
    if (mg%n_cycle_down < 0) mg%n_cycle_down = 2
    if (mg%n_cycle_up < 0) mg%n_cycle_up = 2

    ! Check whether methods are set, otherwise use default (for laplacian)
    if (.not. associated(mg%box_op)) mg%box_op => mg_auto_op
    ! if (.not. associated(mg%box_stencil)) mg%box_stencil => mg_auto_stencil
    if (.not. associated(mg%box_gsrb)) mg%box_gsrb => mg_auto_gsrb
    if (.not. associated(mg%box_corr)) mg%box_corr => mg_auto_corr
    if (.not. associated(mg%box_rstr)) mg%box_rstr => mg_auto_rstr
    if (.not. associated(mg%sides_rb)) mg%sides_rb => mg_auto_rb

    if (mg%operator_key == af_stencil_none) then
       mg%operator_key = mg_auto_operator
    end if

    if (mg%prolongation_key == af_stencil_none) then
       mg%prolongation_key = mg_auto_prolongation
    end if

    if (mg%i_lsf /= -1 .and. .not. associated(mg%lsf_dist)) then
       mg%lsf_dist => mg_lsf_dist_linear
    end if

    call mg_set_box_tag_lvl(tree, mg, 1)
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
    call check_mg(mg)
    call coarse_solver_destroy(mg%csolver)
  end subroutine mg_destroy

  subroutine check_mg(mg)
    type(mg_t), intent(in) :: mg           !< Multigrid options
    if (.not. mg%initialized) stop "check_mg: you haven't called mg_init"
  end subroutine check_mg

  !> Perform FAS-FMG cycle (full approximation scheme, full multigrid). Note
  !> that this routine needs valid ghost cells (for i_phi) on input, and gives
  !> back valid ghost cells on output
  subroutine mg_fas_fmg(tree, mg, set_residual, have_guess)
    use m_af_ghostcell, only: af_gc_ids
    use m_af_utils, only: af_boxes_copy_cc
    type(af_t), intent(inout)       :: tree !< Tree to do multigrid on
    type(mg_t), intent(inout)         :: mg   !< Multigrid options
    logical, intent(in)             :: set_residual !< If true, store residual in i_tmp
    logical, intent(in)             :: have_guess   !< If false, start from phi = 0
    integer                         :: lvl

    call check_mg(mg)           ! Check whether mg options are set
    call mg_set_box_tag_tree(tree, mg) ! Make sure box tags are set

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
       call af_gc_ids(tree, tree%lvls(lvl)%ids, [mg%i_phi])

       call mg_fas_vcycle(tree, mg, set_residual .and. &
            lvl == tree%highest_lvl, lvl, standalone=.false.)
    end do

  end subroutine mg_fas_fmg

  !> Perform FAS V-cycle (full approximation scheme). Note that this routine
  !> needs valid ghost cells (for i_phi) on input, and gives back valid ghost
  !> cells on output
  subroutine mg_fas_vcycle(tree, mg, set_residual, highest_lvl, standalone)
    use m_af_ghostcell, only: af_gc_ids
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
       call check_mg(mg)        ! Check whether mg options are set
       call mg_set_box_tag_tree(tree, mg)
    end if

    max_lvl = tree%highest_lvl
    if (present(highest_lvl)) max_lvl = highest_lvl

    do lvl = max_lvl, 2, -1
       ! Downwards relaxation
       call gsrb_boxes(tree, tree%lvls(lvl)%ids, mg, mg_cycle_down)

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
       call af_gc_ids(tree, tree%lvls(lvl)%ids, [mg%i_phi])

       ! Upwards relaxation
       call gsrb_boxes(tree, tree%lvls(lvl)%ids, mg, mg_cycle_up)
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

  end subroutine mg_fas_vcycle

  subroutine solve_coarse_grid(tree, mg)
    use m_af_ghostcell, only: af_gc_ids
    type(af_t), intent(inout) :: tree !< Tree to do multigrid on
    type(mg_t), intent(in)    :: mg   !< Multigrid options

    call coarse_solver_set_rhs_phi(tree, mg)
    call coarse_solver(mg%csolver)
    call coarse_solver_get_phi(tree, mg)

    ! Set ghost cells for the new coarse grid solution
    call af_gc_ids(tree, tree%lvls(1)%ids, [mg%i_phi])
  end subroutine solve_coarse_grid

  !> Fill ghost cells near refinement boundaries which preserves diffusive fluxes.
  subroutine mg_sides_rb(boxes, id, nb, iv)
    use m_af_ghostcell, only: af_gc_prolong_copy
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
    use m_af_ghostcell, only: af_gc_prolong_copy
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

    call af_gc_prolong_copy(boxes, id, nb, iv)

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

  subroutine gsrb_boxes(tree, ids, mg, type_cycle)
    use m_af_ghostcell, only: af_gc_box
    type(af_t), intent(inout) :: tree       !< Tree containing full grid
    type(mg_t), intent(in)   :: mg         !< Multigrid options
    integer, intent(in)        :: ids(:)     !< Operate on these boxes
    integer, intent(in)        :: type_cycle !< Type of cycle to perform
    integer                    :: n, i, n_cycle
    logical                    :: use_corners

    select case (type_cycle)
    case (mg_cycle_down)
       n_cycle = mg%n_cycle_down
    case (mg_cycle_up)
       n_cycle = mg%n_cycle_up
    case default
       error stop "gsrb_boxes: invalid cycle type"
    end select

    !$omp parallel private(n, i)
    do n = 1, 2 * n_cycle
       !$omp do
       do i = 1, size(ids)
          call mg%box_gsrb(tree%boxes(ids(i)), n, mg)
       end do
       !$omp end do

       use_corners = mg%use_corners .or. &
            (type_cycle /= mg_cycle_down .and. n == 2 * n_cycle)

       !$omp do
       do i = 1, size(ids)
          call af_gc_box(tree, ids(i), [mg%i_phi], use_corners)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine gsrb_boxes

  ! Set rhs on coarse grid, restrict phi, and copy i_phi to i_tmp for the
  ! correction later
  subroutine update_coarse(tree, lvl, mg)
    use m_af_utils, only: af_box_add_cc, af_box_copy_cc
    use m_af_ghostcell, only: af_gc_ids
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

    call af_gc_ids(tree, tree%lvls(lvl-1)%ids, [mg%i_phi])

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
    use m_af_ghostcell, only: af_gc_ids
    use m_af_utils, only: af_box_add_cc
    type(af_t), intent(inout) :: tree !< Tree containing full grid
    integer, intent(in)        :: lvl !< Update coarse values at lvl-1
    type(mg_t), intent(in)   :: mg !< Multigrid options
    integer                    :: i, id, p_id

    ! Fill ghost cells here to be sure
    if (lvl == tree%highest_lvl) then
       call af_gc_ids(tree, tree%lvls(lvl)%ids, [mg%i_phi])
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

    call af_gc_ids(tree, tree%lvls(lvl-1)%ids, [mg%i_phi])

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

    select case (mg%operator_key)
    case (mg_auto_operator)
       call mg_box_auto_stencil(box, mg, ix)
    case (mg_normal_box)
       call mg_box_lpl_stencil(box, mg, ix)
    case (mg_lsf_box)
       call mg_box_lsf_stencil(box, mg, ix)
    case (mg_veps_box, mg_ceps_box)
       call mg_box_lpld_stencil(box, mg, ix)
    case default
       print *, "mg%operator_key: ", mg%operator_key
       error stop "mg_store_operator_stencil: unknown stencil"
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

    select case (mg%prolongation_key)
    case (mg_auto_prolongation)
       call mg_box_prolong_auto_stencil(tree%boxes(id), &
            tree%boxes(p_id), mg, ix)
    case (mg_prolong_linear)
       call mg_box_prolong_linear_stencil(tree%boxes(id), &
            tree%boxes(p_id), mg, ix)
    case (mg_prolong_sparse)
       call mg_box_prolong_sparse_stencil(tree%boxes(id), &
            tree%boxes(p_id), mg, ix)
    case (mg_prolong_lsf)
       call mg_box_prolong_lsf_stencil(tree%boxes(id), &
            tree%boxes(p_id), mg, ix)
    case (mg_prolong_eps)
       call mg_box_prolong_eps_stencil(tree%boxes(id), &
            tree%boxes(p_id), mg, ix)
    case default
       print *, "mg%prolongation_key: ", mg%prolongation_key
       error stop "mg_store_prolongation_stencil: unknown stencil"
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
  subroutine mg_auto_rb(boxes, id, nb, iv)
    type(box_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)        :: id     !< Id of box
    integer, intent(in)        :: nb     !< Ghost cell direction
    integer, intent(in)        :: iv     !< Ghost cell variable

    select case(boxes(id)%tag)
    case (mg_normal_box, mg_lsf_box, af_init_tag)
       ! Use default method; this is also useful for initial refinement when box
       ! tags are not yet set
       call mg_sides_rb(boxes, id, nb, iv)
    case (mg_veps_box, mg_ceps_box)
       ! With a dielectric, use local extrapolation for ghost cells
       call mg_sides_rb_extrap(boxes, id, nb, iv)
    case default
       error stop "mg_auto_rb: unknown box tag"
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
       gradnorm = numerical_gradient_amplitude(mg%lsf, rr)
       root_mask(IJK) = (abs(box%cc(IJK, mg%i_lsf)) < dmax * gradnorm * &
            mg%lsf_gradient_safety_factor)
    end do; CLOSE_DO
  end subroutine get_possible_lsf_root_mask

  !> Check whether the level-set function could have a root by computing the
  !> numerical gradient
  logical function lsf_root_possible(box, dmax, mg)
    type(box_t), intent(in) :: box  !< Box to operate on
    real(dp), intent(in)    :: dmax !< Maximal distance to consider
    type(mg_t), intent(in)  :: mg   !< Multigrid options
    logical                 :: mask(DTIMES(box%n_cell))

    call get_possible_lsf_root_mask(box, box%n_cell, dmax, mg, mask)
    lsf_root_possible = any(mask)
  end function lsf_root_possible

  subroutine mg_set_box_tag(box, mg)
    type(box_t), intent(inout) :: box !< Box to operate on
    type(mg_t), intent(in)     :: mg  !< Multigrid options
    real(dp)                   :: a, b
    logical                    :: is_lsf, is_deps, is_eps

    is_lsf = .false.
    is_eps = .false.
    is_deps = .false.

    if (mg%i_lsf /= -1) then
       is_lsf = lsf_root_possible(box, norm2(box%dr), mg)
    end if

    if (mg%i_eps /= -1) then
       a = minval(box%cc(DTIMES(:), mg%i_eps))
       b = maxval(box%cc(DTIMES(:), mg%i_eps))
       is_deps = (b > a)
       if (.not. is_deps) is_eps = (a < 1 .or. a > 1)
    end if

    if (count([is_lsf, is_eps, is_deps]) > 1) &
         stop "mg_set_box_tag: Cannot set lsf and eps tag for same box"

    box%tag = mg_normal_box
    if (is_lsf) box%tag = mg_lsf_box
    if (is_eps) box%tag = mg_ceps_box
    if (is_deps) box%tag = mg_veps_box

  end subroutine mg_set_box_tag

  subroutine mg_set_box_tag_lvl(tree, mg, lvl)
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(in)    :: mg
    integer, intent(in)       :: lvl
    integer                   :: i, id, n

    !$omp parallel do private(i, id, n)
    do i = 1, size(tree%lvls(lvl)%ids)
       id = tree%lvls(lvl)%ids(i)
       associate (box => tree%boxes(id))
         if (box%tag == af_init_tag) then
            call mg_set_box_tag(box, mg)
         end if

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
  end subroutine mg_set_box_tag_lvl

  subroutine mg_set_box_tag_tree(tree, mg)
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(in)    :: mg
    integer                   :: lvl

    do lvl = 1, tree%highest_lvl
       call mg_set_box_tag_lvl(tree, mg, lvl)
    end do
  end subroutine mg_set_box_tag_tree

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

  !> Automatically store stencil for a box
  subroutine mg_box_auto_stencil(box, mg, ix)
    type(box_t), intent(inout) :: box
    type(mg_t), intent(in)     :: mg
    integer, intent(in)        :: ix !< Stencil index

    select case(box%tag)
    case (mg_normal_box)
       call mg_box_lpl_stencil(box, mg, ix)
    case (mg_lsf_box)
       call mg_box_lsf_stencil(box, mg, ix)
    case (mg_veps_box, mg_ceps_box)
       call mg_box_lpld_stencil(box, mg, ix)
    case default
       error stop "mg_box_auto_stencil: unknown box tag"
    end select
  end subroutine mg_box_auto_stencil

  !> Store the matrix stencil for each cell of the box. The order of the stencil
  !> is (i, j), (i-1, j), (i+1, j), (i, j-1), (i, j+1) (e.g., -4, 1, 1, 1, 1)
  subroutine mg_box_lpl_stencil(box, mg, ix)
    type(box_t), intent(inout) :: box
    type(mg_t), intent(in)     :: mg
    integer, intent(in)        :: ix !< Stencil index
    real(dp)                   :: inv_dr2(NDIM)
    integer                    :: n_coeff, idim

    box%stencils(ix)%shape    = af_stencil_357
    box%stencils(ix)%constant = .true.
    box%stencils(ix)%cylindrical_gradient = (box%coord_t == af_cyl)
    n_coeff                   = af_stencil_sizes(af_stencil_357)
    inv_dr2                   = 1 / box%dr**2

    allocate(box%stencils(ix)%c(n_coeff))
    do idim = 1, NDIM
       box%stencils(ix)%c(2*idim:2*idim+1) = inv_dr2(idim)
    end do
    box%stencils(ix)%c(1) = -sum(box%stencils(ix)%c(2:)) - mg%helmholtz_lambda

  end subroutine mg_box_lpl_stencil

  !> Automatically store prolongation stencil
  subroutine mg_box_prolong_auto_stencil(box, box_p, mg, ix)
    type(box_t), intent(inout) :: box   !< Current box
    type(box_t), intent(in)    :: box_p !< Parent box
    type(mg_t), intent(in)     :: mg
    integer, intent(in)        :: ix    !< Stencil index

    select case (box%tag)
    case (mg_normal_box)
       call mg_box_prolong_linear_stencil(box, box_p, mg, ix)
    case (mg_lsf_box)
       call mg_box_prolong_lsf_stencil(box, box_p, mg, ix)
    case (mg_ceps_box, mg_veps_box)
       call mg_box_prolong_eps_stencil(box, box_p, mg, ix)
    case default
       error stop "mg_box_prolong_auto_stencil: unknown box tag"
    end select
  end subroutine mg_box_prolong_auto_stencil

  !> Store linear prolongation stencil for standard Laplacian
  subroutine mg_box_prolong_linear_stencil(box, box_p, mg, ix)
    type(box_t), intent(inout) :: box   !< Current box
    type(box_t), intent(in)    :: box_p !< Parent box
    type(mg_t), intent(in)     :: mg
    integer, intent(in)        :: ix    !< Stencil index
    integer                    :: n_coeff

    box%stencils(ix)%shape    = af_stencil_p248
    box%stencils(ix)%constant = .true.
    n_coeff                   = af_stencil_sizes(af_stencil_p248)

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

    box%stencils(ix)%shape    = af_stencil_p234
    box%stencils(ix)%constant = .true.
    n_coeff                   = af_stencil_sizes(af_stencil_p234)

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

    nc                        = box%n_cell
    box%stencils(ix)%shape    = af_stencil_p234
    box%stencils(ix)%constant = .false.
    ix_offset                 = af_get_child_offset(box)
    n_coeff                   = af_stencil_sizes(af_stencil_p234)
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
    real(dp)                   :: dd(NDIM+1)
    logical                    :: root_mask(DTIMES(box%n_cell))
    integer                    :: n_coeff, i_lsf, nc
    logical                    :: has_boundary, success
    integer                    :: IJK, IJK_(c1)
    integer                    :: IJK_(c2), ix_offset(NDIM)

    nc                        = box%n_cell
    box%stencils(ix)%shape    = af_stencil_p234
    box%stencils(ix)%constant = .false.
    ix_offset                 = af_get_child_offset(box)
    n_coeff                   = af_stencil_sizes(af_stencil_p234)
    allocate(box%stencils(ix)%v(n_coeff, DTIMES(nc)))
    has_boundary              = .false.

    call get_possible_lsf_root_mask(box, nc, norm2(box%dr), &
         mg, root_mask)
    i_lsf = mg%i_lsf

    associate (v => box%stencils(ix)%v)
      ! In these loops, we calculate the closest coarse index (_c1), and the
      ! one-but-closest (_c2). The fine cell lies in between.
#if NDIM == 1
      do i = 1, nc
         i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
         i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

         if (root_mask(IJK)) then
            dd(1) = mg%lsf_dist(box, IJK, box_p, i_c1, mg)
            dd(2) = mg%lsf_dist(box, IJK, box_p, i_c2, mg)
         else
            dd(:) = 1.0_dp
         end if

         v(:, IJK) = [3 * dd(2), dd(1)]
         v(:, IJK) = v(:, IJK) / sum(v(:, IJK))
      end do
#elif NDIM == 2
      do j = 1, nc
         j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
         j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
         do i = 1, nc
            i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
            i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

            if (root_mask(IJK)) then
               dd(1) = mg%lsf_dist(box, IJK, box_p, i_c1, j_c1, mg)
               dd(2) = mg%lsf_dist(box, IJK, box_p, i_c2, j_c1, mg)
               dd(3) = mg%lsf_dist(box, IJK, box_p, i_c1, j_c2, mg)
            else
               dd(:) = 1.0_dp
            end if

            v(:, IJK) = [2 * dd(2) * dd(3), dd(1) * dd(2), &
                 dd(1) ** dd(3)]
            v(:, IJK) = v(:, IJK) / sum(v(:, IJK))
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

               if (root_mask(IJK)) then
                  dd(1) = mg%lsf_dist(box, IJK, box_p, i_c1, j_c1, k_c1, mg)
                  dd(2) = mg%lsf_dist(box, IJK, box_p, i_c2, j_c1, k_c1, mg)
                  dd(3) = mg%lsf_dist(box, IJK, box_p, i_c1, j_c2, k_c1, mg)
                  dd(4) = mg%lsf_dist(box, IJK, box_p, i_c1, j_c1, k_c2, mg)
               else
                  dd(:) = 1.0_dp
               end if

               v(:, IJK) = [dd(2) * dd(3) * dd(4), &
                    dd(1) * dd(3) * dd(4), dd(1) * dd(2) * dd(4), &
                    dd(1) * dd(2) * dd(3)]
               v(:, IJK) = v(:, IJK) / sum(v(:, IJK))
            end do
         end do
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

    box%stencils(ix)%shape    = af_stencil_357
    box%stencils(ix)%constant = .false.
    box%stencils(ix)%cylindrical_gradient = (box%coord_t == af_cyl)
    n_coeff                   = af_stencil_sizes(af_stencil_357)

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

  !> Compute distance to boundary starting at point a going to point b, in
  !> the range from [0, 1], with 1 meaning there is no boundary
  function mg_lsf_dist_linear(box_a, IJK_(a), box_b, IJK_(b), mg) result(dist)
    type(box_t), intent(in) :: box_a   !< Box a (start point)
    integer, intent(in)     :: IJK_(a) !< Cell-centered index in box a
    type(box_t), intent(in) :: box_b   ! Box b (end point)
    integer, intent(in)     :: IJK_(b) !< Cell-centered index in box b
    type(mg_t), intent(in)  :: mg
    real(dp)                :: dist

    associate(lsf_a => box_a%cc(IJK_(a), mg%i_lsf), &
         lsf_b => box_b%cc(IJK_(b), mg%i_lsf))
      if (lsf_a * lsf_b < 0) then
         ! There is a boundary between the points
         dist = lsf_a / (lsf_a - lsf_b)
      else
         dist = 1.0_dp
      end if
    end associate
  end function mg_lsf_dist_linear

  !> Find root of f in the interval [a, b]. If f(a) and f(b) have different
  !> signs, apply bisection directly. Else, first find the (assumed to be)
  !> unique local minimum/maximum to determine a bracket. Return relative
  !> location of root, or 1 if there is no root.
  !>
  !> @todo Allow to set tolerance?
  function mg_lsf_dist_gss(box_a, IJK_(a), box_b, IJK_(b), mg) result(dist)
    type(box_t), intent(in) :: box_a   !< Box a (start point)
    integer, intent(in)     :: IJK_(a) !< Cell-centered index in box a
    type(box_t), intent(in) :: box_b   ! Box b (end point)
    integer, intent(in)     :: IJK_(b) !< Cell-centered index in box b
    type(mg_t), intent(in)  :: mg
    real(dp)                :: a(NDIM), b(NDIM), bracket(NDIM, 2)
    real(dp)                :: dist, r_root(NDIM), lsf_a, lsf_b
    real(dp), parameter     :: tol      = 1e-8_dp
    integer, parameter      :: max_iter = 100

    a = af_r_cc(box_a, [IJK_(a)])
    b = af_r_cc(box_b, [IJK_(b)])
    lsf_a = mg%lsf(a)
    lsf_b = mg%lsf(b)

    dist = 1.0_dp

    if (lsf_a * lsf_b < 0) then
       r_root = bisection(mg%lsf, a, b, tol, max_iter)
    else
       ! Determine bracket by finding local minimum/maximum
       bracket = gss(mg%lsf, a, b, &
            minimization=(lsf_a >= 0), tol=tol)

       if (mg%lsf(bracket(:, 1)) * lsf_a > 0) then
          return                ! No root
       else
          r_root = bisection(mg%lsf, a, bracket(:, 1), tol, max_iter)
       end if
    end if

    dist = norm2(r_root - a)/norm2(b-a)
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

  !> Golden-section search. Given a function f with a single local minimum in
  !> the interval [a,b], gss returns a subset interval [c,d] that contains the
  !> minimum with d-c <= tol. Copied from
  !> https://en.wikipedia.org/wiki/Golden-section_search
  function gss(f, in_a, in_b, minimization, tol) result(bracket)
    procedure(mg_func_lsf) :: f
    real(dp), intent(in)  :: in_a(NDIM), in_b(NDIM), tol
    logical, intent(in)   :: minimization
    real(dp)              :: bracket(NDIM, 2)
    real(dp)              :: a(NDIM), b(NDIM), c(NDIM), d(NDIM)
    real(dp)              :: h(NDIM), yc, yd
    real(dp)              :: invphi, invphi2
    integer               :: n, k

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
    end do

    if ((yc < yd) .eqv. minimization) then
       bracket(:, 1) = a
       bracket(:, 2) = d
    else
       bracket(:, 1) = c
       bracket(:, 2) = b
    end if
  end function gss

  !> For a point a, compute value and distance (between 0, 1) of a neighbor b.
  subroutine lsf_dist_val(lsf_a, lsf_b, val_b, boundary_value, dist, val)
    !> Level set function at a
    real(dp), intent(in)  :: lsf_a
    !> Level set function at b
    real(dp), intent(in)  :: lsf_b
    !> Value at b
    real(dp), intent(in)  :: val_b
    !> Boundary value
    real(dp), intent(in)  :: boundary_value
    !> Distance to neighbor point (value between 0 and 1)
    real(dp), intent(out) :: dist
    !> Value at neighbor point
    real(dp), intent(out) :: val

    if (lsf_a * lsf_b < 0) then
       ! There is a boundary between the points
       dist = lsf_a / (lsf_a - lsf_b)
       val  = boundary_value
    else
       ! Simply use the value at b
       dist = 1
       val  = val_b
    end if
  end subroutine lsf_dist_val

  !> Store the matrix stencil for each cell of the box. The order of the stencil
  !> is (i, j), (i-1, j), (i+1, j), (i, j-1), (i, j+1) (e.g., -4, 1, 1, 1, 1)
  subroutine mg_box_lsf_stencil(box, mg, ix)
    type(box_t), intent(inout) :: box
    type(mg_t), intent(in)     :: mg
    integer, intent(in)        :: ix !< Index of stencil
    integer                    :: IJK, n, nc, idim, n_coeff
    real(dp)                   :: dd(2*NDIM), dr2(NDIM)
    logical                    :: success
    logical                    :: root_mask(DTIMES(box%n_cell))
#if NDIM == 2
    real(dp)                   :: tmp
#endif

    nc = box%n_cell
    dr2 = box%dr**2

    box%stencils(ix)%shape    = af_stencil_357
    box%stencils(ix)%constant = .false.
    ! Perform a custom correction in cylindrical coordinates
    box%stencils(ix)%cylindrical_gradient = .false.
    n_coeff                   = af_stencil_sizes(af_stencil_357)

    allocate(box%stencils(ix)%v(n_coeff, DTIMES(nc)))
    allocate(box%stencils(ix)%f(DTIMES(nc)))
    allocate(box%stencils(ix)%bc_correction(DTIMES(nc)))
    box%stencils(ix)%f = 0.0_dp

    call get_possible_lsf_root_mask(box, nc, norm2(box%dr), &
         mg, root_mask)

    do KJI_DO(1, nc)
       if (root_mask(IJK)) then
#if NDIM == 1
          dd(1) = mg%lsf_dist(box, IJK, box, i-1, mg)
          dd(2) = mg%lsf_dist(box, IJK, box, i+1, mg)
#elif NDIM == 2
          dd(1) = mg%lsf_dist(box, IJK, box, i-1, j, mg)
          dd(2) = mg%lsf_dist(box, IJK, box, i+1, j, mg)
          dd(3) = mg%lsf_dist(box, IJK, box, i, j-1, mg)
          dd(4) = mg%lsf_dist(box, IJK, box, i, j+1, mg)
#elif NDIM == 3
          dd(1) = mg%lsf_dist(box, IJK, box, i-1, j, k, mg)
          dd(2) = mg%lsf_dist(box, IJK, box, i+1, j, k, mg)
          dd(3) = mg%lsf_dist(box, IJK, box, i, j-1, k, mg)
          dd(4) = mg%lsf_dist(box, IJK, box, i, j+1, k, mg)
          dd(5) = mg%lsf_dist(box, IJK, box, i, j, k-1, mg)
          dd(6) = mg%lsf_dist(box, IJK, box, i, j, k+1, mg)
#endif
       else
          dd(:) = 1.0_dp
       end if

       ! Generalized Laplacian for neighbors at distance dd * dx
       do idim = 1, NDIM
          box%stencils(ix)%v(1+2*idim-1:1+2*idim, IJK) = &
               [dd(2*idim), dd(2*idim-1)] / &
               (0.5_dp * dr2(idim) * (dd(2*idim-1) + dd(2*idim)) * &
               dd(2*idim-1) * dd(2*idim))
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

    call af_stencil_try_constant(box, ix, epsilon(1.0_dp), success)

    if (success) then
       deallocate(box%stencils(ix)%f)
       deallocate(box%stencils(ix)%bc_correction)
    end if

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
          select case(tree%boxes(id)%tag)
          case (mg_normal_box, mg_ceps_box)
             call mg_box_lpl_gradient(tree, id, mg, i_fc, fac)
          case (mg_lsf_box)
             if (af_has_children(tree%boxes(id))) then
                !> @todo Solution on coarse grid can lead to large gradient due
                !> to inconsistencies with level set function
                call mg_box_lpl_gradient(tree, id, mg, i_fc, fac)
             else
                call mg_box_lpllsf_gradient(tree, id, mg, i_fc, fac)
             end if
          case (mg_veps_box)
             ! Should call dielectric_correct_field_fc afterwards
             call mg_box_lpl_gradient(tree, id, mg, i_fc, fac)
          case (af_init_tag)
             error stop "mg_auto_op: box tag not set"
          case default
             error stop "mg_auto_op: unknown box tag"
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

      if (box%tag == mg_veps_box) then
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
    real(dp)                  :: &
         cc(DTIMES(0:tree%n_cell+1), 2)
    integer                   :: IJK, nc, i_phi, i_lsf
    real(dp)                  :: dd, val, v_a(2), v_b(2)
    integer                   :: grad_sign, nb
    integer                   :: ilo(NDIM), ihi(NDIM)
    integer                   :: olo(NDIM), ohi(NDIM)

    associate(box => tree%boxes(id))
      nc     = box%n_cell
      i_phi  = mg%i_phi
      i_lsf  = mg%i_lsf

      ! Store a copy because we might need to modify phi in ghost cells
      cc = box%cc(DTIMES(:), [i_lsf, i_phi])

      do nb = 1, af_num_neighbors
         ! If lsf value changes sign in ghost cell, use boundary value
         if (af_is_ref_boundary(tree%boxes, id, nb)) then
            call af_get_index_bc_inside(nb, box%n_cell, 1, ilo, ihi)
            call af_get_index_bc_outside(nb, box%n_cell, 1, olo, ohi)
            where (cc(DSLICE(ilo,ihi), 1) * &
                 cc(DSLICE(olo,ohi), 1) <= 0.0_dp)
               cc(DSLICE(olo,ohi), 2) = mg%lsf_boundary_value
            end where
         end if
      end do

#if NDIM == 1
      do i = 1, nc+1
         if (cc(i, 1) > 0) then
            v_a = cc(i, :)
            v_b = cc(i-1, :)
            grad_sign = 1
         else
            v_a = cc(i-1, :)
            v_b = cc(i, :)
            grad_sign = -1
         end if

         call lsf_dist_val(v_a(1), v_b(1), v_b(2), &
              mg%lsf_boundary_value, dd, val)
         box%fc(IJK, 1, i_fc) = grad_sign * fac * (v_a(2) - val) / &
              (box%dr(1) * dd)
      end do
#elif NDIM == 2
      do j = 1, nc
         do i = 1, nc+1
            if (cc(i, j, 1) > 0) then
               v_a = cc(i, j, :)
               v_b = cc(i-1, j, :)
               grad_sign = 1
            else
               v_a = cc(i-1, j, :)
               v_b = cc(i, j, :)
               grad_sign = -1
            end if

            call lsf_dist_val(v_a(1), v_b(1), v_b(2), &
                 mg%lsf_boundary_value, dd, val)
            box%fc(IJK, 1, i_fc) = grad_sign * fac * (v_a(2) - val) / &
                 (box%dr(1) * dd)

            if (cc(j, i, 1) > 0) then
               v_a = cc(j, i, :)
               v_b = cc(j, i-1, :)
               grad_sign = 1
            else
               v_a = cc(j, i-1, :)
               v_b = cc(j, i, :)
               grad_sign = -1
            end if

            call lsf_dist_val(v_a(1), v_b(1), v_b(2), &
                 mg%lsf_boundary_value, dd, val)
            box%fc(j, i, 2, i_fc) = grad_sign * fac * (v_a(2) - val) / &
                 (box%dr(2) * dd)
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            do i = 1, nc+1
               if (cc(i, j, k, 1) > 0) then
                  v_a = cc(i, j, k, :)
                  v_b = cc(i-1, j, k, :)
                  grad_sign = 1
               else
                  v_a = cc(i-1, j, k, :)
                  v_b = cc(i, j, k, :)
                  grad_sign = -1
               end if

               call lsf_dist_val(v_a(1), v_b(1), v_b(2), &
                    mg%lsf_boundary_value, dd, val)
               box%fc(IJK, 1, i_fc) = grad_sign * fac * (v_a(2) - val) / &
                    (box%dr(1) * dd)

               if (cc(j, i, k, 1) > 0) then
                  v_a = cc(j, i, k, :)
                  v_b = cc(j, i-1, k, :)
                  grad_sign = 1
               else
                  v_a = cc(j, i-1, k, :)
                  v_b = cc(j, i, k, :)
                  grad_sign = -1
               end if

               call lsf_dist_val(v_a(1), v_b(1), v_b(2), &
                    mg%lsf_boundary_value, dd, val)
               box%fc(j, i, k, 2, i_fc) = grad_sign * fac * (v_a(2) - val) / &
                    (box%dr(2) * dd)

               if (cc(j, k, i, 1) > 0) then
                  v_a = cc(j, k, i, :)
                  v_b = cc(j, k, i-1, :)
                  grad_sign = 1
               else
                  v_a = cc(j, k, i-1, :)
                  v_b = cc(j, k, i, :)
                  grad_sign = -1
               end if

               call lsf_dist_val(v_a(1), v_b(1), v_b(2), &
                    mg%lsf_boundary_value, dd, val)
               box%fc(j, k, i, 3, i_fc) = grad_sign * fac * (v_a(2) - val) / &
                    (box%dr(3) * dd)
            end do
         end do
      end do
#endif
    end associate
  end subroutine mg_box_lpllsf_gradient

  !> This method checks whether the level set function is properly defined on
  !> the coarse grid
  subroutine check_coarse_representation_lsf(tree, mg)
    type(af_t), intent(in) :: tree
    type(mg_t), intent(in) :: mg
    integer                :: i, id, n_stencils

    do i = 1, size(tree%lvls(1)%ids)
       id = tree%lvls(1)%ids(i)
       n_stencils = tree%boxes(id)%n_stencils
       if (.not. all(tree%boxes(id)%stencils(1:n_stencils)%constant)) exit
    end do

    if (i == size(tree%lvls(1)%ids)+1) then
       print *, "all stencils on coarse grid are constant"
       print *, "you should probably use a finer coarse grid"
       error stop "level set function not resolved on coarse grid"
    end if

  end subroutine check_coarse_representation_lsf

  !> Get amplitude of numerical gradient of level set function
  function numerical_gradient_amplitude(f, r) result(normgrad)
    procedure(mg_func_lsf) :: f
    real(dp), intent(in)   :: r(NDIM)
    real(dp), parameter    :: sqrteps      = sqrt(epsilon(1.0_dp))
    real(dp), parameter    :: min_stepsize = epsilon(1.0_dp)
    real(dp)               :: r_eval(NDIM), gradient(NDIM)
    real(dp)               :: stepsize(NDIM), flo, fhi, normgrad
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

    normgrad = norm2(gradient)
  end function numerical_gradient_amplitude

end module m_af_multigrid
