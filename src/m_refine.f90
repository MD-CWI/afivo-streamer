#include "../afivo/src/cpp_macros.h"
module m_refine
  use m_af_all

  implicit none
  private

  ! The refinement buffer width in cells (around flagged cells)
  integer, public, protected :: refine_buffer_width = 4

  ! The number of steps after which the mesh is updated
  integer, public, protected :: refine_per_steps = 2

  ! The grid spacing will always be larger than this value
  real(dp), protected :: refine_min_dx = 1.0e-7_dp

  ! The grid spacing will always be smaller than this value
  real(dp), protected :: refine_max_dx = 1.0e-3_dp

  ! Refine if alpha*dx is larger than this value
  real(dp), protected :: refine_adx = 1.0_dp

  ! For refinement, use alpha(f * E)/f, where f is this factor
  real(dp), protected :: refine_adx_fac = 1.0_dp

  ! Refine if the curvature in phi is larger than this value
  real(dp), protected :: refine_cphi = 1e99_dp

  ! Allow derefinement if the curvature in phi is smaller than this value
  real(dp), protected :: derefine_cphi = 1e99_dp

  ! Only derefine if grid spacing if smaller than this value
  real(dp), protected :: derefine_dx = 1e-4_dp

  ! Refine around initial conditions up to this time
  real(dp), protected :: refine_init_time = 10e-9_dp

  ! Refine until dx is smaller than this factor times the seed width
  real(dp), protected :: refine_init_fac = 0.25_dp

  ! Minimum electron density for adding grid refinement
  real(dp), protected :: refine_min_dens = -1.0e99_dp

  ! Refine a region up to this grid spacing
  real(dp), protected, allocatable :: refine_regions_dr(:)

  ! Refine regions up to this simulation time
  real(dp), protected, allocatable :: refine_regions_tstop(:)

  ! Minimum coordinate of the refinement regions
  real(dp), protected, allocatable :: refine_regions_rmin(:,:)

  ! Maximum coordinate of the refinement regions
  real(dp), protected, allocatable :: refine_regions_rmax(:,:)

  ! Public methods
  public :: refine_initialize
  public :: refine_routine

contains

  subroutine refine_initialize(cfg)
    use m_config
    use m_gas
    type(CFG_t), intent(inout) :: cfg
    integer                    :: n
    real(dp)                   :: vec(NDIM)
    real(dp), allocatable      :: dbuffer(:)

    call CFG_add_get(cfg, "refine_buffer_width", refine_buffer_width, &
         "The refinement buffer width in cells (around flagged cells)")
    call CFG_add_get(cfg, "refine_per_steps", refine_per_steps, &
         "The number of steps after which the mesh is updated")
    call CFG_add_get(cfg, "refine_min_dx", refine_min_dx, &
         "The grid spacing will always be larger than this value")
    call CFG_add_get(cfg, "refine_max_dx", refine_max_dx, &
         "The grid spacing will always be smaller than this value")

    if (refine_min_dx > refine_max_dx) &
         error stop "Cannot have refine_min_dx < refine_max_dx"

    call CFG_add_get(cfg, "refine_adx", refine_adx, &
         "Refine if alpha*dx is larger than this value")
    call CFG_add_get(cfg, "refine_adx_fac", refine_adx_fac, &
         "For refinement, use alpha(f * E)/f, where f is this factor")
    call CFG_add_get(cfg, "refine_cphi", refine_cphi, &
         "Refine if the curvature in phi is larger than this value")
    call CFG_add_get(cfg, "derefine_cphi", derefine_cphi, &
         "Allow derefinement if the curvature in phi is smaller than this value")
    call CFG_add_get(cfg, "derefine_dx", derefine_dx, &
         "Only derefine if grid spacing if smaller than this value")
    call CFG_add_get(cfg, "refine_init_time", refine_init_time, &
         "Refine around initial conditions up to this time")
    call CFG_add_get(cfg, "refine_init_fac", refine_init_fac, &
         "Refine until dx is smaller than this factor times the seed width")
    call CFG_add_get(cfg, "refine_min_dens", refine_min_dens, &
         "Minimum electron density for adding grid refinement")

    call CFG_add(cfg, "refine_regions_dr", [1.0e99_dp], &
         "Refine regions up to this grid spacing", .true.)
    call CFG_add(cfg, "refine_regions_tstop", [-1.0e99_dp], &
         "Refine regions up to this simulation time", .true.)
    vec = 0.0_dp
    call CFG_add(cfg, "refine_regions_rmin", vec, &
         "Minimum coordinate of the refinement regions", .true.)
    call CFG_add(cfg, "refine_regions_rmax", vec, &
         "Maximum coordinate of the refinement regions", .true.)

    call CFG_get_size(cfg, "refine_regions_dr", n)
    allocate(refine_regions_dr(n))
    allocate(refine_regions_tstop(n))
    allocate(refine_regions_rmin(NDIM, n))
    allocate(refine_regions_rmax(NDIM, n))
    allocate(dbuffer(NDIM * n))

    call CFG_get(cfg, "refine_regions_dr", refine_regions_dr)
    call CFG_get(cfg, "refine_regions_tstop", refine_regions_tstop)
    call CFG_get(cfg, "refine_regions_rmin", dbuffer)
    refine_regions_rmin = reshape(dbuffer, [NDIM, n])
    call CFG_get(cfg, "refine_regions_rmax", dbuffer)
    refine_regions_rmax = reshape(dbuffer, [NDIM, n])

  end subroutine refine_initialize

  ! This routine sets the cell refinement flags for box
  subroutine refine_routine(box, cell_flags)
    use m_streamer
    use m_geometry
    use m_init_cond
    use m_gas
    type(box_t), intent(in) :: box
    ! Refinement flags for the cells of the box
    integer, intent(out)     :: &
         cell_flags(DTIMES(box%n_cell))
    integer                  :: IJK, n, nc
    real(dp)                 :: dx, dx2
    real(dp)                 :: alpha, adx, fld, elec_dens
    real(dp)                 :: dist, rmin(NDIM), rmax(NDIM)

    nc      = box%n_cell
    dx      = box%dr
    dx2     = box%dr**2

    do KJI_DO(1,nc)
       fld   = box%cc(IJK, i_electric_fld) * SI_to_Townsend
       alpha = LT_get_col(ST_td_tbl, i_alpha_eff, refine_adx_fac * fld) * &
            gas_number_density / refine_adx_fac
       adx   = box%dr * alpha
       elec_dens = box%cc(IJK, i_electron)

       if (adx  > refine_adx .and. elec_dens > refine_min_dens) then
          cell_flags(IJK) = af_do_ref
       else if (adx < 0.125_dp * refine_adx .and. dx < derefine_dx) then
          cell_flags(IJK) = af_rm_ref
       else
          cell_flags(IJK) = af_keep_ref
       end if

       ! Refine around the initial conditions
       if (ST_time < refine_init_time) then
          do n = 1, init_conds%n_cond
             dist = GM_dist_line(af_r_cc(box, [IJK]), &
                  init_conds%seed_r0(:, n), &
                  init_conds%seed_r1(:, n), NDIM)
             if (dist - init_conds%seed_width(n) < 2 * dx &
                  .and. box%dr > refine_init_fac * &
                  init_conds%seed_width(n)) then
                cell_flags(IJK) = af_do_ref
             end if
          end do
       end if

    end do; CLOSE_DO

    ! Check fixed refinements
    rmin = box%r_min
    rmax = box%r_min + box%dr * box%n_cell

    do n = 1, size(refine_regions_dr)
       if (ST_time <= refine_regions_tstop(n) .and. &
            dx > refine_regions_dr(n) .and. all(&
            rmax >= refine_regions_rmin(:, n) .and. &
            rmin <= refine_regions_rmax(:, n))) then
          ! Mark just the center cell to prevent refining neighbors
          cell_flags(DTIMES(nc/2)) = af_do_ref
       end if
    end do

    ! Make sure we don't have or get a too fine or too coarse grid
    if (dx > refine_max_dx) then
       cell_flags = af_do_ref
    else if (dx < 2 * refine_min_dx) then
       where(cell_flags == af_do_ref) cell_flags = af_keep_ref
    end if

  end subroutine refine_routine

end module m_refine
