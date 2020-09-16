#include "../src/cpp_macros.h"
!> This module contains the geometric multigrid routines that come with Afivo
module m_af_multigrid
  use m_af_types
  use m_mg_types
  use m_coarse_solver

  implicit none
  private

  public :: mg_t

  public :: mg_init
  public :: mg_destroy
  public :: mg_fas_fmg
  public :: mg_fas_vcycle
  public :: mg_box_op
  public :: mg_box_gsrb
  public :: mg_box_corr
  public :: mg_box_rstr

  ! Automatic selection of operators
  public :: mg_set_box_tag
  public :: mg_auto_op
  public :: mg_auto_stencil
  public :: mg_auto_gsrb
  public :: mg_auto_corr
  public :: mg_auto_rstr

  ! Methods for normal Laplacian
  public :: mg_box_lpl
  public :: mg_box_lpl_stencil
  public :: mg_box_gsrb_lpl
  public :: mg_box_corr_lpl
  public :: mg_box_rstr_lpl
  public :: mg_sides_rb

  ! Methods for Laplacian with jump in coefficient between boxes
  public :: mg_box_lpld
  public :: mg_box_gsrb_lpld
  public :: mg_box_corr_lpld
  public :: mg_box_lpld_stencil

  ! To adjust operator stencils near boundaries
  public :: mg_stencil_handle_boundaries

#if NDIM == 2
  public :: mg_box_clpl
  public :: mg_box_clpl_stencil
  public :: mg_box_gsrb_clpl
  public :: mg_box_clpld
  public :: mg_box_gsrb_clpld
  public :: mg_box_clpld_stencil
#endif

contains

  !> Check multigrid options or set them to default
  subroutine mg_init(tree, mg)
    use m_af_core, only: af_set_cc_methods
    type(af_t), intent(inout) :: tree !< Tree to do multigrid on
    type(mg_t), intent(inout) :: mg   !< Multigrid options

    if (mg%i_phi < 0)                  stop "mg_init: i_phi not set"
    if (mg%i_tmp < 0)                  stop "mg_init: i_tmp not set"
    if (mg%i_rhs < 0)                  stop "mg_init: i_rhs not set"
    if (mg%i_lsf * mg%i_bval < 0) &
         stop "mg_init: you have to set both i_lsf and i_bval"

    if (.not. associated(mg%sides_bc)) stop "mg_init: sides_bc not set"
    if (.not. tree%ready) error stop "mg_init: tree not initialized"

    ! Check whether these are set, otherwise use default
    if (mg%n_cycle_down < 0)           mg%n_cycle_down = 2
    if (mg%n_cycle_up < 0)             mg%n_cycle_up = 2

    ! Check whether methods are set, otherwise use default (for laplacian)
    if (.not. associated(mg%sides_rb)) then
       if (mg%i_eps > 0) then
          ! With a dielectric, use local extrapolation for ghost cells
          mg%sides_rb => mg_sides_rb_extrap
       else
          mg%sides_rb => mg_sides_rb
       end if
    end if

    if (.not. associated(mg%box_op)) mg%box_op => mg_auto_op
    if (.not. associated(mg%box_stencil)) mg%box_stencil => mg_auto_stencil
    if (.not. associated(mg%box_gsrb)) mg%box_gsrb => mg_auto_gsrb
    if (.not. associated(mg%box_corr)) mg%box_corr => mg_auto_corr
    if (.not. associated(mg%box_rstr)) mg%box_rstr => mg_auto_rstr

    call mg_set_box_tag_lvl(tree, mg, 1)
    call coarse_solver_initialize(tree, mg)

    ! Set the proper methods for the phi variable
    call af_set_cc_methods(tree, mg%i_phi, mg%sides_bc, mg%sides_rb)

    mg%initialized = .true.

  end subroutine mg_init

  subroutine mg_destroy(mg)
    type(mg_t), intent(inout) :: mg   !< Multigrid options

    if (.not. mg%initialized) stop "check_mg: you haven't called mg_init"
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

    do lvl = 1, tree%highest_lvl
       ! Store phi_old in tmp
       call af_boxes_copy_cc(tree%boxes, tree%lvls(lvl)%ids, &
            mg%i_phi, mg%i_tmp)

       if (lvl > 1) then
          ! Correct solution at this lvl using lvl-1 data
          ! phi = phi + prolong(phi_coarse - phi_old_coarse)
          call correct_children(tree%boxes, tree%lvls(lvl-1)%parents, mg)

          ! Update ghost cells
          call af_gc_ids(tree, tree%lvls(lvl)%ids, [mg%i_phi])
       end if

       ! Perform V-cycle, only set residual on last iteration
       call mg_fas_vcycle(tree, mg, &
            set_residual .and. lvl == tree%highest_lvl, lvl, standalone=.false.)
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

    call check_mg(mg)           ! Check whether mg options are set

    by_itself = .true.; if (present(standalone)) by_itself = standalone
    if (by_itself) call mg_set_box_tag_tree(tree, mg)

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
          call af_gc_box(tree, ids(i), [mg%i_phi])
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
    type(box_t), intent(inout) :: box !< Box to operate on
    integer, intent(in)         :: redblack_cntr !< Iteration count
    type(mg_t), intent(in)     :: mg !< Multigrid options

    select case(box%tag)
#if NDIM == 2
    case (mg_normal_box)
       if (box%coord_t == af_cyl) then
          call mg_box_gsrb_clpl(box, redblack_cntr, mg)
       else
          call mg_box_gsrb_lpl(box, redblack_cntr, mg)
       end if
    ! case (mg_lsf_box)
    !    call mg_box_gsrb_lpllsf(box, redblack_cntr, mg)
    case (mg_veps_box, mg_ceps_box)
       if (box%coord_t == af_cyl) then
          call mg_box_gsrb_clpld(box, redblack_cntr, mg)
       else
          call mg_box_gsrb_lpld(box, redblack_cntr, mg)
       end if
#else
    case (mg_normal_box)
       call mg_box_gsrb_lpl(box, redblack_cntr, mg)
    ! case (mg_lsf_box)
    !    call mg_box_gsrb_lpllsf(box, redblack_cntr, mg)
    case (mg_veps_box, mg_ceps_box)
       call mg_box_gsrb_lpld(box, redblack_cntr, mg)
#endif
    case default
       error stop "mg_auto_gsrb: unknown box tag"
    end select
  end subroutine mg_auto_gsrb

  !> Based on the box type, apply the approriate operator
  subroutine mg_auto_op(box, i_out, mg)
    type(box_t), intent(inout) :: box   !< Operate on this box
    integer, intent(in)        :: i_out !< Index of output variable
    type(mg_t), intent(in)     :: mg    !< Multigrid options

    select case(box%tag)
#if NDIM == 2
    case (mg_normal_box)
       if (box%coord_t == af_cyl) then
          call mg_box_clpl(box, i_out, mg)
       else
          call mg_box_lpl(box, i_out, mg)
       end if
    ! case (mg_lsf_box)
    !    call mg_box_lpllsf(box, i_out, mg)
    case (mg_veps_box, mg_ceps_box)
       if (box%coord_t == af_cyl) then
          call mg_box_clpld(box, i_out, mg)
       else
          call mg_box_lpld(box, i_out, mg)
       end if
#else
    case (mg_normal_box)
       call mg_box_lpl(box, i_out, mg)
    ! case (mg_lsf_box)
    !    call mg_box_lpllsf(box, i_out, mg)
    case (mg_veps_box, mg_ceps_box)
       call mg_box_lpld(box, i_out, mg)
#endif
    case (af_init_tag)
       error stop "mg_auto_op: box tag not set"
    case default
       error stop "mg_auto_op: unknown box tag"
    end select
  end subroutine mg_auto_op

  !> Based on the box type, use the appropriate stencil method
  subroutine mg_auto_stencil(box, mg, stencil, bc_to_rhs)
    type(box_t), intent(in) :: box
    type(mg_t), intent(in)  :: mg
    real(dp), intent(inout) :: stencil(2*NDIM+1, DTIMES(box%n_cell))
    real(dp), intent(inout) :: bc_to_rhs(box%n_cell**(NDIM-1), af_num_neighbors)

    select case(box%tag)
#if NDIM == 2
    case (mg_normal_box)
       if (box%coord_t == af_cyl) then
          call mg_box_clpl_stencil(box, mg, stencil, bc_to_rhs)
       else
          call mg_box_lpl_stencil(box, mg, stencil, bc_to_rhs)
       end if
    ! case (mg_lsf_box)
    !    call mg_box_lpllsf_stencil(box, mg, stencil, bc_to_rhs)
    case (mg_veps_box, mg_ceps_box)
       if (box%coord_t == af_cyl) then
          call mg_box_clpld_stencil(box, mg, stencil, bc_to_rhs)
       else
          call mg_box_lpld_stencil(box, mg, stencil, bc_to_rhs)
       end if
#else
    case (mg_normal_box)
       call mg_box_lpl_stencil(box, mg, stencil, bc_to_rhs)
    ! case (mg_lsf_box)
    !    call mg_box_lpllsf_stencil(box, mg, stencil, bc_to_rhs)
    case (mg_veps_box, mg_ceps_box)
       call mg_box_lpld_stencil(box, mg, stencil, bc_to_rhs)
#endif
    case (af_init_tag)
       error stop "mg_auto_stencil: box tag not set"
    case default
       error stop "mg_auto_stencil: unknown box tag"
    end select
  end subroutine mg_auto_stencil

  !> Based on the box type, apply the approriate Laplace operator
  subroutine mg_auto_rstr(box_c, box_p, iv, mg)
    type(box_t), intent(in)    :: box_c !< Child box
    type(box_t), intent(inout) :: box_p !< Parent box
    integer, intent(in)        :: iv    !< Index of variable
    type(mg_t), intent(in)     :: mg    !< Multigrid options

    select case(box_c%tag)
    case (mg_normal_box, mg_veps_box, mg_ceps_box)
       call mg_box_rstr_lpl(box_c, box_p, iv, mg)
    ! case (mg_lsf_box)
       !    call mg_box_rstr_lpllsf(box_c, box_p, iv, mg)
    case (af_init_tag)
       error stop "mg_auto_rstr: box_c tag not set"
    case default
       error stop "mg_auto_rstr: unknown box tag"
    end select
  end subroutine mg_auto_rstr

  !> Based on the box type, correct the solution of the children
  subroutine mg_auto_corr(box_p, box_c, mg)
    type(box_t), intent(inout) :: box_c !< Child box
    type(box_t), intent(in)    :: box_p !< Parent box
    type(mg_t), intent(in)     :: mg !< Multigrid options

    select case(box_c%tag)
    case (mg_normal_box)
       call mg_box_corr_lpl(box_p, box_c, mg)
    ! case (mg_lsf_box)
    !    call mg_box_corr_lpllsf(box_p, box_c, mg)
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
    integer                   :: i, id

    !$omp parallel do private(i, id)
    do i = 1, size(tree%lvls(lvl)%ids)
       id = tree%lvls(lvl)%ids(i)
       if (tree%boxes(id)%tag == af_init_tag) then
          call mg_set_box_tag(tree%boxes(id), mg)
       end if
    end do
    !$omp end parallel do
  end subroutine mg_set_box_tag_lvl

  subroutine mg_set_box_tag_tree(tree, mg)
    type(af_t), intent(inout) :: tree
    type(mg_t), intent(in)    :: mg
    integer                   :: lvl, i, id

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          if (tree%boxes(id)%tag == af_init_tag) then
             call mg_set_box_tag(tree%boxes(id), mg)
          end if
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine mg_set_box_tag_tree

  subroutine mg_box_corr_lpl(box_p, box_c, mg)
    use m_af_prolong
    type(box_t), intent(inout) :: box_c !< Child box
    type(box_t), intent(in)    :: box_p !< Parent box
    type(mg_t), intent(in)     :: mg !< Multigrid options

    call af_prolong_linear(box_p, box_c, mg%i_tmp, mg%i_phi, add=.true.)
  end subroutine mg_box_corr_lpl

  !> Perform Gauss-Seidel relaxation on box for a Laplacian operator
  subroutine mg_box_gsrb_lpl(box, redblack_cntr, mg)
    type(box_t), intent(inout) :: box           !< Box to operate on
    integer, intent(in)        :: redblack_cntr !< Iteration counter
    type(mg_t), intent(in)     :: mg            !< Multigrid options
    integer                    :: IJK, i0, nc
    real(dp)                   :: idr2(NDIM), fac

    idr2 = 1 / box%dr**2
    fac  = 1.0_dp / (2 * sum(idr2) + mg%helmholtz_lambda)
    nc   = box%n_cell

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => box%cc, n => mg%i_phi, i_rhs => mg%i_rhs)
#if NDIM == 1
      i0 = 2 - iand(redblack_cntr, 1)
      do i = i0, nc, 2
         cc(i, n) = fac * ( &
              idr2(1) * (cc(i+1, n) + cc(i-1, n)) - &
              cc(i, mg%i_rhs))
      end do
#elif NDIM == 2
      do j = 1, nc
         i0 = 2 - iand(ieor(redblack_cntr, j), 1)
         do i = i0, nc, 2
            cc(i, j, n) = fac * ( &
                 idr2(1) * (cc(i+1, j, n) + cc(i-1, j, n)) + &
                 idr2(2) * (cc(i, j+1, n) + cc(i, j-1, n)) - &
                 cc(i, j, mg%i_rhs))
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            i0 = 2 - iand(ieor(redblack_cntr, k+j), 1)
            do i = i0, nc, 2
               cc(i, j, k, n) = fac * ( &
                    idr2(1) * (cc(i+1, j, k, n) + cc(i-1, j, k, n)) + &
                    idr2(2) * (cc(i, j+1, k, n) + cc(i, j-1, k, n)) + &
                    idr2(3) * (cc(i, j, k+1, n) + cc(i, j, k-1, n)) - &
                    cc(i, j, k, i_rhs))
            end do
         end do
      end do
#endif
    end associate
  end subroutine mg_box_gsrb_lpl

  !> Perform Laplacian operator on a box
  subroutine mg_box_lpl(box, i_out, mg)
    type(box_t), intent(inout) :: box   !< Box to operate on
    integer, intent(in)        :: i_out !< Index of variable to store Laplacian in
    type(mg_t), intent(in)     :: mg    !< Multigrid options
    integer                    :: IJK, nc
    real(dp)                   :: idr2(NDIM)

    nc   = box%n_cell
    idr2 = 1 / box%dr**2

    associate (cc => box%cc, n => mg%i_phi)
      do KJI_DO(1, nc)
#if NDIM == 1
         cc(i, i_out) = &
              idr2(1) * (cc(i-1, n) + cc(i+1, n) - 2 * cc(i, n)) - &
              mg%helmholtz_lambda * cc(i, n)
#elif NDIM == 2
         cc(i, j, i_out) = &
              idr2(1) * (cc(i-1, j, n) + cc(i+1, j, n) - 2 * cc(i, j, n)) + &
              idr2(2) * (cc(i, j-1, n) + cc(i, j+1, n) - 2 * cc(i, j, n)) - &
              mg%helmholtz_lambda * cc(i, j, n)
#elif NDIM == 3
         cc(i, j, k, i_out) = &
              idr2(1) * (cc(i-1, j, k, n) + cc(i+1, j, k, n) &
              - 2 * cc(i, j, k, n)) &
              + idr2(2) * (cc(i, j-1, k, n) + cc(i, j+1, k, n) &
              - 2 * cc(i, j, k, n)) &
              + idr2(3) * (cc(i, j, k-1, n) + cc(i, j, k+1, n) &
              - 2 * cc(i, j, k, n)) - &
              mg%helmholtz_lambda * cc(i, j, k, n)
#endif
      end do; CLOSE_DO
    end associate
  end subroutine mg_box_lpl

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
  subroutine mg_box_lpl_stencil(box, mg, stencil, bc_to_rhs)
    type(box_t), intent(in) :: box
    type(mg_t), intent(in)  :: mg
    real(dp), intent(inout) :: stencil(2*NDIM+1, DTIMES(box%n_cell))
    real(dp), intent(inout) :: bc_to_rhs(box%n_cell**(NDIM-1), af_num_neighbors)
    real(dp)                :: inv_dr2(NDIM)
    integer                 :: idim

    inv_dr2                 = 1 / box%dr**2
    stencil(1, DTIMES(:))   = -2.0_dp * sum(inv_dr2) - mg%helmholtz_lambda
    do idim = 1, NDIM
       stencil(2*idim:2*idim+1, DTIMES(:)) = inv_dr2(idim)
    end do
    call mg_stencil_handle_boundaries(box, mg, stencil, bc_to_rhs)
  end subroutine mg_box_lpl_stencil

  !> Incorporate boundary conditions into stencil
  subroutine mg_stencil_handle_boundaries(box, mg, stencil, bc_to_rhs)
    use m_af_ghostcell, only: af_gc_get_boundary_coords
    type(box_t), intent(in) :: box
    type(mg_t), intent(in)  :: mg
    real(dp), intent(inout) :: stencil(2*NDIM+1, DTIMES(box%n_cell))
    real(dp), intent(inout) :: bc_to_rhs(box%n_cell**(NDIM-1), af_num_neighbors)
    integer                 :: nb, nc, lo(NDIM), hi(NDIM)
    integer                 :: nb_id, nb_dim, bc_type
    real(dp)                :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp)                :: bc_val(box%n_cell**(NDIM-1))

    bc_to_rhs = 0.0_dp
    nc        = box%n_cell

    do nb = 1, af_num_neighbors
       nb_id = box%neighbors(nb)

       if (nb_id < af_no_box) then
          nb_dim = af_neighb_dim(nb)
          call af_gc_get_boundary_coords(box, nb, coords)
          call mg%sides_bc(box, nb, mg%i_phi, coords, bc_val, bc_type)

          ! Determine index range next to boundary
          call af_get_index_bc_inside(nb, nc, 1, lo, hi)

          select case (bc_type)
          case (af_bc_dirichlet)
             ! Dirichlet value at cell face, so compute gradient over h/2
             ! E.g. 1 -2 1 becomes 0 -3 1 for a 1D Laplacian
             ! The boundary condition is incorporated in the right-hand side
             stencil(1, DSLICE(lo, hi)) = &
                  stencil(1, DSLICE(lo, hi)) - &
                  stencil(nb+1, DSLICE(lo, hi))
             bc_to_rhs(:, nb) = pack(-2 * stencil(nb+1, DSLICE(lo, hi)), .true.)
             stencil(nb+1, DSLICE(lo, hi)) = 0.0_dp
          case (af_bc_neumann)
             ! E.g. 1 -2 1 becomes 0 -1 1 for a 1D Laplacian
             stencil(1, DSLICE(lo, hi)) = &
                  stencil(1, DSLICE(lo, hi)) + &
                  stencil(nb+1, DSLICE(lo, hi))
             bc_to_rhs(:, nb) = &
                  -pack(stencil(nb+1, DSLICE(lo, hi)) * &
                  box%dr(nb_dim), .true.) * af_neighb_high_pm(nb)
             stencil(nb+1, DSLICE(lo, hi)) = 0.0_dp
          case default
             error stop "mg_box_lpl_stencil: unsupported boundary condition"
          end select
       end if
    end do

  end subroutine mg_stencil_handle_boundaries

  !> Perform Gauss-Seidel relaxation on a box. Epsilon can have a jump at cell
  !> faces.
  subroutine mg_box_gsrb_lpld(box, redblack_cntr, mg)
    type(box_t), intent(inout) :: box            !< Box to operate on
    integer, intent(in)        :: redblack_cntr !< Iteration counter
    type(mg_t), intent(in)     :: mg             !< Multigrid options
    integer                    :: IJK, i0, nc
    real(dp)                   :: u(2*NDIM), a0
    real(dp)                   :: idr2(2*NDIM), a(2*NDIM), c(2*NDIM)

    idr2(1:2*NDIM:2) = 1/box%dr**2
    idr2(2:2*NDIM:2) = idr2(1:2*NDIM:2)
    nc    = box%n_cell

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => box%cc, n => mg%i_phi, i_rhs => mg%i_rhs, i_eps => mg%i_eps)
#if NDIM == 1
      i0 = 2 - iand(redblack_cntr, 1)
      do i = i0, nc, 2
         a0     = cc(i, i_eps)
         u(1:2) = cc(i-1:i+1:2, n)
         a(1:2) = cc(i-1:i+1:2, i_eps)
         c(:)   = 2 * a0 * a(:) / (a0 + a(:)) * idr2

         cc(i, n) = (sum(c(:) * u(:)) - cc(i, i_rhs)) / sum(c(:))
      end do
#elif NDIM == 2
      do j = 1, nc
         i0 = 2 - iand(ieor(redblack_cntr, j), 1)
         do i = i0, nc, 2
            a0     = cc(i, j, i_eps)
            u(1:2) = cc(i-1:i+1:2, j, n)
            a(1:2) = cc(i-1:i+1:2, j, i_eps)
            u(3:4) = cc(i, j-1:j+1:2, n)
            a(3:4) = cc(i, j-1:j+1:2, i_eps)
            c(:)   = 2 * a0 * a(:) / (a0 + a(:)) * idr2

            cc(i, j, n) = &
                 (sum(c(:) * u(:)) - cc(i, j, i_rhs)) / sum(c(:))
         end do
      end do
#elif NDIM == 3
      do k = 1, nc
         do j = 1, nc
            i0 = 2 - iand(ieor(redblack_cntr, k+j), 1)
            do i = i0, nc, 2
               a0     = cc(i, j, k, i_eps)
               u(1:2) = cc(i-1:i+1:2, j, k, n)
               a(1:2) = cc(i-1:i+1:2, j, k, i_eps)
               u(3:4) = cc(i, j-1:j+1:2, k, n)
               a(3:4) = cc(i, j-1:j+1:2, k, i_eps)
               u(5:6) = cc(i, j, k-1:k+1:2, n)
               a(5:6) = cc(i, j, k-1:k+1:2, i_eps)
               c(:)   = 2 * a0 * a(:) / (a0 + a(:)) * idr2

               cc(i, j, k, n) = (sum(c(:) * u(:)) - &
                    cc(i, j, k, i_rhs)) / sum(c(:))
            end do
         end do
      end do
#endif
    end associate
  end subroutine mg_box_gsrb_lpld

  !> Perform Laplacian operator on a box where epsilon varies on cell faces
  subroutine mg_box_lpld(box, i_out, mg)
    type(box_t), intent(inout) :: box   !< Box to operate on
    integer, intent(in)        :: i_out !< Index of variable to store Laplacian in
    type(mg_t), intent(in)     :: mg    !< Multigrid options
    integer                    :: IJK, nc
    real(dp)                   :: idr2(2*NDIM), a0, u0
    real(dp)                   :: u(2*NDIM), a(2*NDIM)

    nc               = box%n_cell
    idr2(1:2*NDIM:2) = 1/box%dr**2
    idr2(2:2*NDIM:2) = idr2(1:2*NDIM:2)

    associate (cc => box%cc, n => mg%i_phi, i_eps => mg%i_eps)
      do KJI_DO(1, nc)
#if NDIM == 1
         a0     = cc(i, i_eps)
         a(1:2) = cc(i-1:i+1:2, i_eps)
         u0     = cc(i, n)
         u(1:2) = cc(i-1:i+1:2, n)

         cc(i, i_out) = sum(2 * idr2 * &
              a0*a(:)/(a0 + a(:)) * (u(:) - u0))

#elif NDIM == 2
         a0     = cc(i, j, i_eps)
         a(1:2) = cc(i-1:i+1:2, j, i_eps)
         a(3:4) = cc(i, j-1:j+1:2, i_eps)
         u0     = cc(i, j, n)
         u(1:2) = cc(i-1:i+1:2, j, n)
         u(3:4) = cc(i, j-1:j+1:2, n)

         cc(i, j, i_out) = sum(2 * idr2 * &
              a0*a(:)/(a0 + a(:)) * (u(:) - u0))

#elif NDIM == 3
         u0 = cc(i, j, k, n)
         a0 = cc(i, j, k, i_eps)
         u(1:2) = cc(i-1:i+1:2, j, k, n)
         u(3:4) = cc(i, j-1:j+1:2, k, n)
         u(5:6) = cc(i, j, k-1:k+1:2, n)
         a(1:2) = cc(i-1:i+1:2, j, k, i_eps)
         a(3:4) = cc(i, j-1:j+1:2, k, i_eps)
         a(5:6) = cc(i, j, k-1:k+1:2, i_eps)

         cc(i, j, k, i_out) = sum(2 * idr2 * &
              a0*a(:)/(a0 + a(:)) * (u(:) - u0))
#endif
      end do; CLOSE_DO
    end associate
  end subroutine mg_box_lpld

  !> Store the matrix stencil for each cell of the box. The order of the stencil
  !> is (i, j), (i-1, j), (i+1, j), (i, j-1), (i, j+1) (e.g., -4, 1, 1, 1, 1)
  subroutine mg_box_lpld_stencil(box, mg, stencil, bc_to_rhs)
    type(box_t), intent(in) :: box
    type(mg_t), intent(in)  :: mg
    real(dp), intent(inout) :: stencil(2*NDIM+1, DTIMES(box%n_cell))
    real(dp), intent(inout) :: bc_to_rhs(box%n_cell**(NDIM-1), af_num_neighbors)
    integer                 :: IJK, nc
    real(dp)                :: idr2(2*NDIM), a0, a(2*NDIM)

    nc               = box%n_cell
    idr2(1:2*NDIM:2) = 1/box%dr**2
    idr2(2:2*NDIM:2) = idr2(1:2*NDIM:2)

    associate (cc => box%cc, n => mg%i_phi, i_eps => mg%i_eps)
      do KJI_DO(1, nc)
#if NDIM == 1
         a0 = box%cc(i, i_eps)
         a(1:2) = box%cc(i-1:i+1:2, i_eps)

         stencil(2:, IJK) = idr2 * 2 * a0*a(:)/(a0 + a(:))
         stencil(1, IJK) = -sum(stencil(2:, IJK))
#elif NDIM == 2
         a0 = box%cc(i, j, i_eps)
         a(1:2) = box%cc(i-1:i+1:2, j, i_eps)
         a(3:4) = box%cc(i, j-1:j+1:2, i_eps)

         stencil(2:, IJK) = idr2 * 2 * a0*a(:)/(a0 + a(:))
         stencil(1, IJK) = -sum(stencil(2:, IJK))
#elif NDIM == 3
         a0 = box%cc(i, j, k, i_eps)
         a(1:2) = box%cc(i-1:i+1:2, j, k, i_eps)
         a(3:4) = box%cc(i, j-1:j+1:2, k, i_eps)
         a(5:6) = box%cc(i, j, k-1:k+1:2, i_eps)

         stencil(2:, IJK) = idr2 * 2 * a0*a(:)/(a0 + a(:))
         stencil(1, IJK) = -sum(stencil(2:, IJK))
#endif
      end do; CLOSE_DO
    end associate

    call mg_stencil_handle_boundaries(box, mg, stencil, bc_to_rhs)
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

#if NDIM == 2
  !> Perform Gauss-Seidel relaxation on box for a cylindrical Laplacian operator
  subroutine mg_box_gsrb_clpl(box, redblack_cntr, mg)
    type(box_t), intent(inout) :: box           !< Box to operate on
    integer, intent(in)        :: redblack_cntr !< Iteration counter
    type(mg_t), intent(in)     :: mg            !< Multigrid options
    integer                    :: i, i0, j, nc
    real(dp)                   :: idr2(NDIM), fac, rfac(2, box%n_cell)

    nc   = box%n_cell
    idr2 = 1/box%dr**2
    fac  = 1.0_dp / (2 * sum(idr2) + mg%helmholtz_lambda)
    call af_cyl_flux_factors(box, rfac)

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => box%cc, n => mg%i_phi, i_rhs => mg%i_rhs)
      do j = 1, nc
         i0 = 2 - iand(ieor(redblack_cntr, j), 1)
         do i = i0, nc, 2
            cc(i, j, n) = fac * (&
                 idr2(1) * (rfac(1, i) * cc(i-1, j, n) + &
                 rfac(2, i) * cc(i+1, j, n)) + &
                 idr2(2) * (cc(i, j+1, n) + cc(i, j-1, n)) &
                 - cc(i, j, i_rhs))
         end do
      end do
    end associate
  end subroutine mg_box_gsrb_clpl

  !> Perform cylindrical Laplacian operator on a box
  subroutine mg_box_clpl(box, i_out, mg)
    type(box_t), intent(inout) :: box   !< Box to operate on
    integer, intent(in)        :: i_out !< Index of variable to store Laplacian in
    type(mg_t), intent(in)     :: mg    !< Multigrid options
    integer                    :: i, j, nc
    real(dp)                   :: idr2(NDIM), rfac(2, box%n_cell)

    nc    = box%n_cell
    idr2  = 1 / box%dr**2
    call af_cyl_flux_factors(box, rfac)

    associate (cc => box%cc, n => mg%i_phi)
      do j = 1, nc
         do i = 1, nc
            cc(i, j, i_out) = idr2(1) * ( rfac(1, i) * cc(i-1, j, n) + &
                 rfac(2, i) * cc(i+1, j, n) - 2 * cc(i, j, n)) + &
                 idr2(2) * (cc(i, j-1, n) + cc(i, j+1, n) - 2 * cc(i, j, n)) &
                 - mg%helmholtz_lambda * cc(i, j, n)
         end do
      end do
    end associate
  end subroutine mg_box_clpl

  !> Store the matrix stencil for each cell of the box. The order of the stencil
  !> is (i, j), (i-1, j), (i+1, j), (i, j-1), (i, j+1) (e.g., -4, 1, 1, 1, 1)
  subroutine mg_box_clpl_stencil(box, mg, stencil, bc_to_rhs)
    type(box_t), intent(in) :: box
    type(mg_t), intent(in)  :: mg
    real(dp), intent(inout) :: stencil(2*NDIM+1, DTIMES(box%n_cell))
    real(dp), intent(inout) :: bc_to_rhs(box%n_cell**(NDIM-1), af_num_neighbors)
    integer                 :: i, j, nc
    real(dp)                :: idr2(NDIM), rfac(2, box%n_cell)

    nc   = box%n_cell
    idr2 = 1 / box%dr**2
    call af_cyl_flux_factors(box, rfac)

    do j = 1, nc
       do i = 1, nc
          stencil(1, i, j)  = -2 * sum(idr2) - mg%helmholtz_lambda
          stencil(2:3, i, j) = idr2(1) * rfac(:, i)
          stencil(4:5, i, j) = idr2(2)
       end do
    end do

    call mg_stencil_handle_boundaries(box, mg, stencil, bc_to_rhs)
  end subroutine mg_box_clpl_stencil

  !> Perform cylindrical Laplacian operator on a box with varying eps
  subroutine mg_box_clpld(box, i_out, mg)
    type(box_t), intent(inout) :: box   !< Box to operate on
    integer, intent(in)        :: i_out !< Index of variable to store Laplacian in
    type(mg_t), intent(in)     :: mg    !< Multigrid options
    integer                    :: i, j, nc, i_phi, i_eps
    real(dp)                   :: idr2(2*NDIM), a0, u0, u(4), a(4)
    real(dp)                   :: rfac(2, box%n_cell), r_weight(4)

    nc               = box%n_cell
    idr2(1:2*NDIM:2) = 1/box%dr**2
    idr2(2:2*NDIM:2) = idr2(1:2*NDIM:2)
    i_phi            = mg%i_phi
    i_eps            = mg%i_eps
    call af_cyl_flux_factors(box, rfac)

    do j = 1, nc
       do i = 1, nc
          r_weight(1:2) = rfac(:, i)
          r_weight(3:4) = 1.0_dp
          u0 = box%cc(i, j, i_phi)
          a0 = box%cc(i, j, i_eps)
          u(1:2) = box%cc(i-1:i+1:2, j, i_phi)
          u(3:4) = box%cc(i, j-1:j+1:2, i_phi)
          a(1:2) = box%cc(i-1:i+1:2, j, i_eps)
          a(3:4) = box%cc(i, j-1:j+1:2, i_eps)

          box%cc(i, j, i_out) =  2 * &
               sum(idr2 * r_weight*a0*a(:)/(a0 + a(:)) * (u(:) - u0))
       end do
    end do
  end subroutine mg_box_clpld

  !> Store the matrix stencil for each cell of the box. The order of the stencil
  !> is (i, j), (i-1, j), (i+1, j), (i, j-1), (i, j+1) (e.g., -4, 1, 1, 1, 1)
  subroutine mg_box_clpld_stencil(box, mg, stencil, bc_to_rhs)
    type(box_t), intent(in) :: box
    type(mg_t), intent(in)  :: mg
    real(dp), intent(inout) :: stencil(2*NDIM+1, DTIMES(box%n_cell))
    real(dp), intent(inout) :: bc_to_rhs(box%n_cell**(NDIM-1), af_num_neighbors)
    integer                 :: i, j, nc, i_eps
    real(dp)                :: idr2(2*NDIM), a0, a(4)
    real(dp)                :: rfac(2, box%n_cell), r_weight(4)

    nc               = box%n_cell
    i_eps            = mg%i_eps
    idr2(1:2*NDIM:2) = 1/box%dr**2
    idr2(2:2*NDIM:2) = idr2(1:2*NDIM:2)
    call af_cyl_flux_factors(box, rfac)

    do j = 1, nc
       do i = 1, nc
          r_weight(1:2) = rfac(:, i)
          r_weight(3:4) = 1.0_dp
          a0 = box%cc(i, j, i_eps)
          a(1:2) = box%cc(i-1:i+1:2, j, i_eps)
          a(3:4) = box%cc(i, j-1:j+1:2, i_eps)

          stencil(2:, i, j) = 2 * idr2 * r_weight*a0*a(:)/(a0 + a(:))
          stencil(1, i, j)  = -sum(stencil(2:, i, j))
       end do
    end do

    call mg_stencil_handle_boundaries(box, mg, stencil, bc_to_rhs)
  end subroutine mg_box_clpld_stencil

  !> Perform Gauss-Seidel relaxation on box for a cylindrical Laplacian operator
  !> with a changing eps
  subroutine mg_box_gsrb_clpld(box, redblack_cntr, mg)
    type(box_t), intent(inout) :: box           !< Box to operate on
    integer, intent(in)        :: redblack_cntr !< Iteration counter
    type(mg_t), intent(in)     :: mg            !< Multigrid options
    integer                    :: i, i0, j, nc
    real(dp)                   :: u(4), a0, a(4), c(4), idr2(2*NDIM)
    real(dp)                   :: rfac(2, box%n_cell), r_weight(4)

    nc               = box%n_cell
    idr2(1:2*NDIM:2) = 1/box%dr**2
    idr2(2:2*NDIM:2) = idr2(1:2*NDIM:2)
    call af_cyl_flux_factors(box, rfac)

    ! The parity of redblack_cntr determines which cells we use. If
    ! redblack_cntr is even, we use the even cells and vice versa.
    associate (cc => box%cc, n => mg%i_phi, i_eps => mg%i_eps, i_rhs => mg%i_rhs)
      do j = 1, nc
         i0 = 2 - iand(ieor(redblack_cntr, j), 1)
         do i = i0, nc, 2
            r_weight(1:2) = rfac(:, i)
            r_weight(3:4) = 1.0_dp
            a0 = box%cc(i, j, i_eps) ! value of eps at i,j
            u(1:2) = box%cc(i-1:i+1:2, j, n) ! values at neighbors
            a(1:2) = box%cc(i-1:i+1:2, j, i_eps)
            u(3:4) = box%cc(i, j-1:j+1:2, n)
            a(3:4) = box%cc(i, j-1:j+1:2, i_eps)
            c(:) = 2 * a0 * a(:) / (a0 + a(:)) * idr2 * r_weight

            box%cc(i, j, n) = (sum(c(:) *  u(:)) - box%cc(i, j, i_rhs)) / sum(c(:))
         end do
      end do
    end associate
  end subroutine mg_box_gsrb_clpld
#endif

end module m_af_multigrid
