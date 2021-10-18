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
    if (.not. associated(mg%box_corr)) mg%box_corr => mg_auto_corr_v2
    if (.not. associated(mg%box_rstr)) mg%box_rstr => mg_auto_rstr
    if (.not. associated(mg%sides_rb)) mg%sides_rb => mg_auto_rb

    !> @todo improve how operator is set
    if (mg%operator_key == af_stencil_none) then
       mg%operator_key = mg_normal_box
       if (mg%i_lsf /= -1) mg%operator_key = mg_lsf_box
       if (mg%i_eps /= -1) mg%operator_key = mg_veps_box

       if (mg%i_lsf /= -1 .and. mg%i_eps /= -1) &
            error stop "TODO: implement mixing of dielectrics/electrodes"
    end if

    if (mg%prolongation_key == af_stencil_none) then
       select case (mg%operator_key)
       case (mg_normal_box)
          mg%prolongation_key = mg_prolong_linear
       case (mg_lsf_box)
          mg%prolongation_key = mg_prolong_lsf
       case (mg_veps_box, mg_ceps_box)
          mg%prolongation_key = mg_prolong_eps
       case default
          error stop "Not implemented yet"
       end select
    end if

    if (mg%i_lsf /= -1) then
       call check_coarse_representation_lsf(tree, mg)
    end if

    call mg_set_box_tag_lvl(tree, mg, 1)
    call coarse_solver_initialize(tree, mg)

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
    use m_af_utils, only: af_n_cell, af_box_add_cc, af_box_copy_cc
    use m_af_ghostcell, only: af_gc_ids
    type(af_t), intent(inout) :: tree !< Tree containing full grid
    integer, intent(in)        :: lvl !< Update coarse values at lvl-1
    type(mg_t), intent(in)   :: mg !< Multigrid options
    integer                    :: i, id, p_id, nc
    real(dp), allocatable :: tmp(DTIMES(:))

    id = tree%lvls(lvl)%ids(1)
    nc = af_n_cell(tree, lvl)
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
    case (mg_normal_box)
       call mg_box_lpl_stencil(box, mg, ix)
    case (mg_lsf_box)
       call mg_box_lpllsf_stencil(box, mg, ix)
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

  !> Based on the box type, apply the approriate Laplace operator
  subroutine mg_auto_rstr(box_c, box_p, iv, mg)
    type(box_t), intent(in)    :: box_c !< Child box
    type(box_t), intent(inout) :: box_p !< Parent box
    integer, intent(in)        :: iv    !< Index of variable
    type(mg_t), intent(in)     :: mg    !< Multigrid options

    select case(box_c%tag)
    case (mg_normal_box, mg_veps_box, mg_ceps_box, mg_lsf_box)
       call mg_box_rstr_lpl(box_c, box_p, iv, mg)
    case (af_init_tag)
       error stop "mg_auto_rstr: box_c tag not set"
    case default
       error stop "mg_auto_rstr: unknown box tag"
    end select
  end subroutine mg_auto_rstr

  !> Set ghost cells near refinement boundaries
  subroutine mg_auto_rb(boxes, id, nb, iv)
    type(box_t), intent(inout) :: boxes(:) !< List of all boxes
    integer, intent(in)        :: id     !< Id of box
    integer, intent(in)        :: nb     !< Ghost cell direction
    integer, intent(in)        :: iv     !< Ghost cell variable

    select case(boxes(id)%tag)
    case (mg_normal_box)
       call mg_sides_rb(boxes, id, nb, iv)
    case (mg_veps_box, mg_ceps_box)
       ! With a dielectric, use local extrapolation for ghost cells
       call mg_sides_rb_extrap(boxes, id, nb, iv)
    case (mg_lsf_box)
       call mg_sides_rb(boxes, id, nb, iv)
    case (af_init_tag)
       ! Use default method; this is useful for initial refinement when box tags
       ! are not yet set
       call mg_sides_rb(boxes, id, nb, iv)
    case default
       error stop "mg_auto_rb: unknown box tag"
    end select

  end subroutine mg_auto_rb

  !> Based on the box type, correct the solution of the children
  subroutine mg_auto_corr_v2(box_p, box_c, mg)
    type(box_t), intent(inout) :: box_c !< Child box
    type(box_t), intent(in)    :: box_p !< Parent box
    type(mg_t), intent(in)     :: mg !< Multigrid options

    call af_stencil_prolong_box(box_p, box_c, mg%prolongation_key, &
         mg%i_tmp, mg%i_phi, .true.)
  end subroutine mg_auto_corr_v2

  !> Based on the box type, correct the solution of the children
  subroutine mg_auto_corr(box_p, box_c, mg)
    type(box_t), intent(inout) :: box_c !< Child box
    type(box_t), intent(in)    :: box_p !< Parent box
    type(mg_t), intent(in)     :: mg !< Multigrid options

    select case(box_c%tag)
    case (mg_normal_box)
       call mg_box_corr_lpl(box_p, box_c, mg)
    case (mg_lsf_box)
       call mg_box_corr_lpllsf(box_p, box_c, mg)
    case (mg_veps_box, mg_ceps_box)
       call mg_box_corr_lpld(box_p, box_c, mg)
    case (af_init_tag)
       error stop "mg_auto_corr: box_c tag not set"
    case default
       error stop "mg_auto_corr: unknown box tag"
    end select
  end subroutine mg_auto_corr

  subroutine mg_set_box_tag(box, mg)
    type(box_t), intent(inout) :: box !< Box to operate on
    type(mg_t), intent(in)     :: mg  !< Multigrid options
    real(dp)                   :: a, b
    logical                    :: is_lsf, is_deps, is_eps

    is_lsf = .false.
    is_eps = .false.
    is_deps = .false.

    if (mg%i_lsf /= -1) then
       is_lsf = minval(box%cc(DTIMES(:), mg%i_lsf)) * &
            maxval(box%cc(DTIMES(:), mg%i_lsf)) < 0
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

         ! Store required stencils
         if (mg%operator_key == af_stencil_none) &
              error stop "mg%operator_key not defined"

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

  subroutine mg_box_corr_lpl(box_p, box_c, mg)
    use m_af_prolong
    type(box_t), intent(inout) :: box_c !< Child box
    type(box_t), intent(in)    :: box_p !< Parent box
    type(mg_t), intent(in)     :: mg !< Multigrid options

    call af_prolong_linear(box_p, box_c, mg%i_tmp, mg%i_phi, add=.true.)
  end subroutine mg_box_corr_lpl

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
    real(dp)                   :: dd(NDIM+1), lsf_a
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

    i_lsf = mg%i_lsf

    associate (v => box%stencils(ix)%v)
      ! In these loops, we calculate the closest coarse index (_c1), and the
      ! one-but-closest (_c2). The fine cell lies in between.
#if NDIM == 1
      do i = 1, nc
         i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
         i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

         lsf_a = box%cc(i, i_lsf)
         dd(1) = lsf_dist(lsf_a, box_p%cc(i_c1, i_lsf))
         dd(2) = lsf_dist(lsf_a, box_p%cc(i_c2, i_lsf))
         if (any(dd < 1)) has_boundary = .true.

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

            lsf_a = box%cc(i, j, i_lsf)
            dd(1) = lsf_dist(lsf_a, box_p%cc(i_c1, j_c1, i_lsf))
            dd(2) = lsf_dist(lsf_a, box_p%cc(i_c2, j_c1, i_lsf))
            dd(3) = lsf_dist(lsf_a, box_p%cc(i_c1, j_c2, i_lsf))

            if (any(dd < 1)) has_boundary = .true.

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

               lsf_a = box%cc(i, j, k, i_lsf)
               dd(1) = lsf_dist(lsf_a, box_p%cc(i_c1, j_c1, k_c1, i_lsf))
               dd(2) = lsf_dist(lsf_a, box_p%cc(i_c2, j_c1, k_c1, i_lsf))
               dd(3) = lsf_dist(lsf_a, box_p%cc(i_c1, j_c2, k_c1, i_lsf))
               dd(4) = lsf_dist(lsf_a, box_p%cc(i_c1, j_c1, k_c2, i_lsf))

               if (any(dd < 1)) has_boundary = .true.

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

  !> Correct fine grid values based on the change in the coarse grid, in the
  !> case of a jump in epsilon
  subroutine mg_box_corr_lpld(box_p, box_c, mg)
    type(box_t), intent(inout)  :: box_c !< Child box
    type(box_t), intent(in)     :: box_p !< Parent box
    type(mg_t), intent(in)      :: mg !< Multigrid options
    integer                      :: ix_offset(NDIM), i_phi, i_corr, i_eps
    integer                      :: nc, IJK, IJK_(c1), IJK_(c2)
    real(dp)                     :: u0, u(NDIM), a0, a(NDIM)
#if NDIM == 3
    real(dp), parameter          :: third = 1/3.0_dp
#endif

    nc = box_c%n_cell
    ix_offset = af_get_child_offset(box_c)
    i_phi = mg%i_phi
    i_corr = mg%i_tmp
    i_eps = mg%i_eps

    ! In these loops, we calculate the closest coarse index (_c1), and the
    ! one-but-closest (_c2). The fine cell lies in between.
#if NDIM == 1
    do i = 1, nc
       i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
       i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

       u0 = box_p%cc(i_c1, i_corr)
       a0 = box_p%cc(i_c1, i_eps)
       u(1) = box_p%cc(i_c2, i_corr)
       a(1) = box_p%cc(i_c2, i_eps)

       ! Get value of phi at coarse cell faces, and average
       box_c%cc(i, i_phi) = box_c%cc(i, i_phi) + &
            sum( (a0*u0 + a(:)*u(:)) / (a0 + a(:)) )
    end do
#elif NDIM == 2
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

          ! Get value of phi at coarse cell faces, and average
          box_c%cc(i, j, i_phi) = box_c%cc(i, j, i_phi) + 0.5_dp * &
               sum( (a0*u0 + a(:)*u(:)) / (a0 + a(:)) )
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

             u0 = box_p%cc(i_c1, j_c1, k_c1, i_corr)
             u(1) = box_p%cc(i_c2, j_c1, k_c1, i_corr)
             u(2) = box_p%cc(i_c1, j_c2, k_c1, i_corr)
             u(3) = box_p%cc(i_c1, j_c1, k_c2, i_corr)
             a0 = box_p%cc(i_c1, j_c1, k_c1, i_eps)
             a(1) = box_p%cc(i_c2, j_c1, k_c1, i_eps)
             a(2) = box_p%cc(i_c1, j_c2, k_c1, i_eps)
             a(3) = box_p%cc(i_c1, j_c1, k_c2, i_eps)

             ! Get value of phi at coarse cell faces, and average
             box_c%cc(i, j, k, i_phi) = box_c%cc(i, j, k, i_phi) + third * &
                  sum((a0*u0 + a(:) * (1.5_dp * u(:) - 0.5_dp * u0)) / &
                  (a0 + a(:)))
          end do
       end do
    end do
#endif
  end subroutine mg_box_corr_lpld

  !> For a point a, compute distance (between 0, 1) of a neighbor b.
  elemental function lsf_dist(lsf_a, lsf_b) result(distance)
    !> Level set function at a
    real(dp), intent(in) :: lsf_a
    !> Level set function at b
    real(dp), intent(in) :: lsf_b
    !> Distance to neighbor point b (value between 0 and 1)
    real(dp) :: distance

    if (lsf_a * lsf_b < 0) then
       ! There is a boundary between the points
       distance = lsf_a / (lsf_a - lsf_b)
    else
       distance = 1
    end if
  end function lsf_dist

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

  subroutine mg_box_corr_lpllsf(box_p, box_c, mg)
    type(box_t), intent(inout) :: box_c !< Child box
    type(box_t), intent(in)    :: box_p !< Parent box
    type(mg_t), intent(in)     :: mg    !< Multigrid options
    integer                    :: i_phi, i_corr, i_lsf, ix_offset(NDIM)
    integer                    :: nc
    integer                    :: IJK, IJK_(c1), IJK_(c2)
    real(dp)                   :: lsf_a, v_b(2)
    real(dp)                   :: val(NDIM+1), dist(NDIM+1), c(NDIM+1)

    nc        = box_c%n_cell
    ix_offset = af_get_child_offset(box_c)
    i_phi     = mg%i_phi
    i_corr    = mg%i_tmp
    i_lsf     = mg%i_lsf

    ! In these loops, we calculate the closest coarse index (_c1), and the
    ! one-but-closest (_c2). The fine cell lies in between.
#if NDIM == 1
    do i = 1, nc
       i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
       i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

       lsf_a = box_c%cc(i, i_lsf)
       v_b(1:2) = box_p%cc(i_c1, [i_lsf, i_corr])
       call lsf_dist_val(lsf_a, v_b(1), v_b(2), 0.0_dp, dist(1), val(1))
       v_b(1:2) = box_p%cc(i_c2, [i_lsf, i_corr])
       call lsf_dist_val(lsf_a, v_b(1), v_b(2), 0.0_dp, dist(2), val(2))

       ! Interpolate between the two coarse values, but account for distance to
       ! potential boundary
       c(1) = 3 * dist(2)
       c(2) = dist(1)
       box_c%cc(i, i_phi) = box_c%cc(i, i_phi) + sum(c * val)/sum(c)
    end do
#elif NDIM == 2
    do j = 1, nc
       j_c1 = ix_offset(2) + ishft(j+1, -1) ! (j+1)/2
       j_c2 = j_c1 + 1 - 2 * iand(j, 1)     ! even: +1, odd: -1
       do i = 1, nc
          i_c1 = ix_offset(1) + ishft(i+1, -1) ! (i+1)/2
          i_c2 = i_c1 + 1 - 2 * iand(i, 1)     ! even: +1, odd: -1

          lsf_a = box_c%cc(i, j, i_lsf)
          v_b(1:2) = box_p%cc(i_c1, j_c1, [i_lsf, i_corr])
          call lsf_dist_val(lsf_a, v_b(1), v_b(2), 0.0_dp, dist(1), val(1))
          v_b(1:2) = box_p%cc(i_c2, j_c1, [i_lsf, i_corr])
          call lsf_dist_val(lsf_a, v_b(1), v_b(2), 0.0_dp, dist(2), val(2))
          v_b(1:2) = box_p%cc(i_c1, j_c2, [i_lsf, i_corr])
          call lsf_dist_val(lsf_a, v_b(1), v_b(2), 0.0_dp, dist(3), val(3))

          ! This expresses general interpolation between 3 points (on the lines
          ! between the fine and the 3 coarse values).
          c(1) = 2 * dist(2) * dist(3)
          c(2) = dist(1) * dist(3)
          c(3) = dist(1) * dist(2)
          box_c%cc(i, j, i_phi) = box_c%cc(i, j, i_phi) + sum(c * val)/sum(c)
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

             lsf_a = box_c%cc(i, j, k, i_lsf)
             v_b(1:2) = box_p%cc(i_c1, j_c1, k_c1, [i_lsf, i_corr])
             call lsf_dist_val(lsf_a, v_b(1), v_b(2), 0.0_dp, dist(1), val(1))
             v_b(1:2) = box_p%cc(i_c2, j_c1, k_c1, [i_lsf, i_corr])
             call lsf_dist_val(lsf_a, v_b(1), v_b(2), 0.0_dp, dist(2), val(2))
             v_b(1:2) = box_p%cc(i_c1, j_c2, k_c1, [i_lsf, i_corr])
             call lsf_dist_val(lsf_a, v_b(1), v_b(2), 0.0_dp, dist(3), val(3))
             v_b(1:2) = box_p%cc(i_c1, j_c1, k_c2, [i_lsf, i_corr])
             call lsf_dist_val(lsf_a, v_b(1), v_b(2), 0.0_dp, dist(4), val(4))

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
  end subroutine mg_box_corr_lpllsf

  !> Store the matrix stencil for each cell of the box. The order of the stencil
  !> is (i, j), (i-1, j), (i+1, j), (i, j-1), (i, j+1) (e.g., -4, 1, 1, 1, 1)
  subroutine mg_box_lpllsf_stencil(box, mg, ix)
    type(box_t), intent(inout) :: box
    type(mg_t), intent(in)     :: mg
    integer, intent(in)        :: ix !< Index of stencil
    integer                    :: IJK, n, nc, idim, n_coeff
    real(dp)                   :: lsf_a, dd(2*NDIM), dr2(NDIM)
    logical                    :: has_boundary, success
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

    has_boundary = .false.

    associate (cc => box%cc, i_lsf => mg%i_lsf)
      do KJI_DO(1, nc)
         lsf_a = box%cc(IJK, i_lsf)
#if NDIM == 1
         dd(1) = lsf_dist(lsf_a, box%cc(i-1, i_lsf))
         dd(2) = lsf_dist(lsf_a, box%cc(i+1, i_lsf))
#elif NDIM == 2
         dd(1) = lsf_dist(lsf_a, box%cc(i-1, j, i_lsf))
         dd(2) = lsf_dist(lsf_a, box%cc(i+1, j, i_lsf))
         dd(3) = lsf_dist(lsf_a, box%cc(i, j-1, i_lsf))
         dd(4) = lsf_dist(lsf_a, box%cc(i, j+1, i_lsf))
#elif NDIM == 3
         dd(1) = lsf_dist(lsf_a, box%cc(i-1, j, k, i_lsf))
         dd(2) = lsf_dist(lsf_a, box%cc(i+1, j, k, i_lsf))
         dd(3) = lsf_dist(lsf_a, box%cc(i, j-1, k, i_lsf))
         dd(4) = lsf_dist(lsf_a, box%cc(i, j+1, k, i_lsf))
         dd(5) = lsf_dist(lsf_a, box%cc(i, j, k-1, i_lsf))
         dd(6) = lsf_dist(lsf_a, box%cc(i, j, k+1, i_lsf))
#endif
         if (any(dd < 1)) has_boundary = .true.

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
    end associate

    call af_stencil_try_constant(box, ix, epsilon(1.0_dp), success)

    if (success) then
       deallocate(box%stencils(ix)%f)
       deallocate(box%stencils(ix)%bc_correction)
    end if

  end subroutine mg_box_lpllsf_stencil

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
    integer                :: i, id, nc
    real(dp)               :: max_lsf, min_lsf

    nc      = tree%n_cell
    max_lsf = -1e100_dp
    min_lsf = 1e100_dp

    do i = 1, size(tree%lvls(1)%ids)
       id = tree%lvls(1)%ids(i)
       max_lsf = max(max_lsf, &
            maxval(tree%boxes(id)%cc(DTIMES(1:nc), mg%i_lsf)))
       min_lsf = min(min_lsf, &
            minval(tree%boxes(id)%cc(DTIMES(1:nc), mg%i_lsf)))
    end do

    if (max_lsf * min_lsf >= 0) then
       print *, "Min/max level set function:", min_lsf, max_lsf
       print *, "Perhaps make the coarse grid finer?"
       error stop "Level set function does not change sign on coarse grid"
    end if

  end subroutine check_coarse_representation_lsf

end module m_af_multigrid
