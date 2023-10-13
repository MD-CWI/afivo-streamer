#include "../afivo/src/cpp_macros.h"
!> Module that provides routines for Monte Carlo photoionization
module m_photoi_mc
  use m_streamer
  use m_lookup_table
  use m_af_all
  use m_units_constants

  implicit none
  private

  !> Type to quickly look up absorption lengths from a table
  type, public :: phmc_tbl_t
     type(LT_t) :: tbl         !< The lookup table
     real(dp)   :: frac_in_tbl !< Fraction photons in table
  end type phmc_tbl_t

  !> Whether physical photons are used
  logical, public, protected :: phmc_physical_photons = .true.

  !> Whether a constant or adaptive grid spacing is usedn
  logical, protected :: phmc_const_dx = .true.

  !> Minimum grid spacing for photoionization
  real(dp), protected :: phmc_min_dx = 1e-9_dp

  !> At which grid spacing photons are absorbed compared to their mean distance
  real(dp), protected :: phmc_absorp_fac = 0.25_dp

  !> Number of photons to use
  integer, protected :: phmc_num_photons = 5000*1000

  !> Minimal photon weight
  real(dp), protected :: phmc_min_weight = 1.0_dp

  !> Table for photoionization
  type(phmc_tbl_t), public, protected :: phmc_tbl

  ! Public methods
  public :: phmc_initialize
  public :: phmc_print_grid_spacing
  public :: phmc_get_table_air
  public :: phmc_do_absorption
  public :: phmc_absorption_func_air
  public :: phmc_set_src
  ! public :: phe_mc_set_src

contains

  !> Initialize photoionization parameters
  subroutine phmc_initialize(cfg, is_used)
    use m_config
    use m_streamer
    use m_gas
    type(CFG_t), intent(inout) :: cfg !< The configuration for the simulation
    logical, intent(in)        :: is_used
    integer                    :: ix

    call CFG_add_get(cfg, "photoi_mc%physical_photons", &
         phmc_physical_photons, &
         "Whether physical photons are used")
    call CFG_add_get(cfg, "photoi_mc%min_weight", phmc_min_weight, &
         "Minimal photon weight (default: 1.0)")
    call CFG_add_get(cfg, "photoi_mc%const_dx", &
         phmc_const_dx, &
         "Whether a constant grid spacing is used for photoionization")
    call CFG_add_get(cfg, "photoi_mc%min_dx", &
         phmc_min_dx, &
         "Minimum grid spacing for photoionization")
    call CFG_add_get(cfg, "photoi_mc%absorp_fac", phmc_absorp_fac, &
         "At which grid spacing photons are absorbed compared to their mean distance")
    if (phmc_absorp_fac <= 0.0_dp) error stop "photoi_mc%absorp_fac <= 0.0"
    call CFG_add_get(cfg, "photoi_mc%num_photons", phmc_num_photons, &
         "Maximum number of discrete photons to use")
    if (phmc_num_photons < 1) error stop "photoi_mc%num_photons < 1"

    if (is_used) then
       ! Create table for photoionization
       ix = gas_index("O2")
       if (ix == -1) error stop "Photoionization: no oxygen present"
       call phmc_get_table_air(phmc_tbl, gas_fractions(ix) * gas_pressure, &
            2 * maxval(ST_domain_len))

       call phmc_print_grid_spacing(maxval(ST_domain_len/ST_box_size))

       if (ST_use_dielectric) then
          if (.not. phmc_const_dx) &
               error stop "phmc_initialize: with dielectric use const_dx"
          if (phmc_absorp_fac > 1e-9_dp) &
               error stop "Use phmc_absorp_fac <= 1e-9 with dielectric"
       end if
    end if

  end subroutine phmc_initialize

  !> Print the grid spacing used for the absorption of photons
  subroutine phmc_print_grid_spacing(dr_base)
    real(dp), intent(in)           :: dr_base
    real(dp)                       :: lengthscale, dx
    integer                        :: lvl

    ! Get a typical length scale for the absorption of photons
    lengthscale = LT_get_col(phmc_tbl%tbl, 1, phmc_absorp_fac)

    ! Determine at which level we estimate the photoionization source term. This
    ! depends on the typical length scale for absorption.
    lvl = get_lvl_length(dr_base, lengthscale)
    dx  = dr_base * (0.5_dp**(lvl-1))

    if (phmc_const_dx) then
       write(*, "(A,E12.3)") " Monte-Carlo photoionization spacing: ", dx
    else
       write(*, "(A)") " Monte-Carlo photoionization uses adaptive spacing"
    end if
  end subroutine phmc_print_grid_spacing

  !> Compute the photonization table for air. If the absorption function is
  !> f(r), this table contains r as a function of F (the cumulative absorption
  !> function, between 0 and 1). Later on a distance can be sampled by drawing a
  !> uniform 0,1 random number and finding the corresponding distance. The table
  !> constructed up to max_dist; we can ignore photons that fly very far.
  subroutine phmc_get_table_air(phmc_tbl, p_O2, max_dist)
    !< The photonization table
    type(phmc_tbl_t), intent(inout) :: phmc_tbl
    !< Partial pressure of oxygen (bar)
    real(dp), intent(in)              :: p_O2
    !< Maximum distance in lookup table
    real(dp), intent(in)  :: max_dist

    integer               :: n
    integer, parameter    :: tbl_size  = 500
    real(dp), allocatable :: fsum(:), dist(:)
    real(dp)              :: dF, drdF, r, F, Fmax_guess

    ! First estimate which fraction of photons are within max_dist (start with
    ! the upper bound of 1.0)
    Fmax_guess = 1.0_dp

    ! 5 loops should be enough for a good guess
    do n = 1, 5
       ! When Fmax_guess decreases so will dF, since we keep the number of
       ! points the same
       dF = Fmax_guess / (tbl_size-1)
       r  = 0
       F  = 0

       do
          drdF = rk4_drdF(r, dF, p_O2)
          r = r + dF * drdF
          F = F + df

          if (r > max_dist) then
             ! A better estimate for the upper bound
             Fmax_guess = F
             exit
          end if
       end do
    end do

    ! Make arrays larger so that we surely are in range of maxdist
    allocate(fsum(2 * tbl_size))
    allocate(dist(2 * tbl_size))

    ! Now create table
    dF = Fmax_guess / (tbl_size-1)
    dist(1) = 0
    fsum(1) = 0

    ! Compute r(F) for F = dF, 2 dF, 3 dF, ...
    do n = 2, 2 * tbl_size
       drdF = rk4_drdF(dist(n-1), dF, p_O2)
       fsum(n) = fsum(n-1) + dF
       dist(n) = dist(n-1) + dF * drdF
       if (dist(n) > max_dist) exit
    end do

    if (n > tbl_size + 10) &
         stop "phmc_get_table_air: integration accuracy fail"

    if (ST_use_dielectric) then
       ! Photons that fly outside the domain can be absorbed by the dielectric,
       ! so don't apply any tricks
       phmc_tbl%frac_in_tbl = 1.0_dp
    else
       ! Scale table to lie between 0 and 1 (later on we have to correct for this)
       phmc_tbl%frac_in_tbl = fsum(n-1)
       fsum(1:n-1) = fsum(1:n-1) / fsum(n-1)
    end if

    phmc_tbl%tbl = LT_create(0.0_dp, 1.0_dp, tbl_size, 1)
    call LT_set_col(phmc_tbl%tbl, 1, fsum(1:n-1), dist(1:n-1))

    write(*, "(A,E12.3,A)") " Average photon absorption length", &
         1e3_dp * sum(dist(1:n-1))/(n-1), " mm"
  end subroutine phmc_get_table_air

  !> Runge Kutta 4 method to compute dr/dF, where r is the absorption distance
  !> and F is the cumulative absorption function (between 0 and 1)
  real(dp) function rk4_drdF(r, dF, p_O2)
    real(dp), intent(in) :: r    !> Initial point
    real(dp), intent(in) :: dF   !> grid size
    real(dp), intent(in) :: p_O2 !< Partial pressure of oxygen (bar)
    real(dp)             :: drdF
    real(dp)             :: sum_drdF
    real(dp), parameter  :: one_sixth = 1 / 6.0_dp

    ! Step 1 (at initial r)
    drdF = 1 / phmc_absorption_func_air(r, p_O2)
    sum_drdF = drdF

    ! Step 2 (at initial r + dr/2)
    drdF = 1 / phmc_absorption_func_air(r + 0.5_dp * dF * drdF, p_O2)
    sum_drdF = sum_drdF + 2 * drdF

    ! Step 3 (at initial r + dr/2)
    drdF = 1 / phmc_absorption_func_air(r + 0.5_dp * dF * drdF, p_O2)
    sum_drdF = sum_drdF + 2 * drdF

    ! Step 4 (at initial r + dr)
    drdF = 1 / phmc_absorption_func_air(r + dF * drdF, p_O2)
    sum_drdF = sum_drdF + drdF

    ! Combine r derivatives at steps
    rk4_drdF = one_sixth * sum_drdF
  end function rk4_drdF

  !> The absorption function for photons in air according to Zheleznyak's model
  real(dp) function phmc_absorption_func_air(dist, p_O2)
    use m_units_constants
    real(dp), intent(in) :: dist   !< Distance
    real(dp), intent(in) :: p_O2   !< Partial pressure of oxygen (bar)
    real(dp)             :: r
    real(dp), parameter  :: c0 = 3.5_dp / UC_torr_to_bar
    real(dp), parameter  :: c1 = 200 / UC_torr_to_bar
    real(dp), parameter  :: eps = epsilon(1.0_dp)

    r = p_O2 * dist
    if (r * (c0 + c1) < eps) then
       ! Use limit to prevent over/underflow
       phmc_absorption_func_air = (c1 - c0 + 0.5_dp * (c0**2 - c1**2) * r) &
            * p_O2 / log(c1/c0)
    else if (r * c0 > -log(eps)) then
       ! Use limit to prevent over/underflow
       phmc_absorption_func_air = eps
    else
       phmc_absorption_func_air = (exp(-c0 * r) - exp(-c1 * r)) / (dist * log(c1/c0))
    end if
  end function phmc_absorption_func_air

  !> Determine the lowest level at which the grid spacing is smaller than 'length'.
  integer function get_lvl_length(dr_base, length)
    real(dp), intent(in) :: dr_base !< cell spacing at lvl 1
    real(dp), intent(in) :: length  !< Some length
    real(dp), parameter :: invlog2 = 1 / log(2.0_dp)
    real(dp) :: ratio

    ratio = dr_base / length
    if (ratio <= 1) then
       get_lvl_length = 1
    else
       get_lvl_length = 1 + ceiling(log(ratio) * invlog2)
    end if
  end function get_lvl_length

  !> As get_lvl_length but with a random choice between lvl and lvl-1
  integer function get_rlvl_length(dr_base, length, rng)
    use m_random
    real(dp), intent(in) :: dr_base   !< cell spacing at lvl 1
    real(dp), intent(in) :: length    !< Some length
    type(RNG_t), intent(inout) :: rng !< Random number generator
    real(dp), parameter :: invlog2 = 1 / log(2.0_dp)
    real(dp) :: ratio, tmp

    ratio = dr_base / length
    if (ratio <= 1) then
       get_rlvl_length = 1
    else
       tmp = log(ratio) * invlog2
       get_rlvl_length = floor(tmp)
       if (rng%unif_01() < tmp - get_rlvl_length) &
            get_rlvl_length = get_rlvl_length + 1
    end if
  end function get_rlvl_length

  !> Given a list of photon production positions (xyz_in), compute where they
  !> end up (xyz_out).
  subroutine phmc_do_absorption(xyz_in, xyz_out, n_dim, n_photons, tbl, prng)
    use m_lookup_table
    use m_random
    use omp_lib
    integer, intent(in)        :: n_photons !< Number of photons
    !> Input (x,y,z) values
    real(dp), intent(in)       :: xyz_in(3, n_photons)
    !> Output (x,y,z) values
    real(dp), intent(out)      :: xyz_out(3, n_photons)
    integer, intent(in)        :: n_dim     !< 2 or 3 dimensional
    !< Lookup table
    type(LT_t), intent(in)     :: tbl
    type(PRNG_t), intent(inout) :: prng       !< Random number generator
    integer                    :: n, proc_id
    real(dp)                   :: rr, dist

    !$omp parallel private(n, rr, dist, proc_id)
    proc_id = 1+omp_get_thread_num()

    if (n_dim == 2) then
       !$omp do
       do n = 1, n_photons
          rr = prng%rngs(proc_id)%unif_01()
          dist = LT_get_col(tbl, 1, rr)
          ! Pick a random point on a sphere, and ignore the last dimension
          xyz_out(1:3, n) =  xyz_in(1:3, n) + &
               prng%rngs(proc_id)%sphere(dist)
       end do
       !$omp end do
    else if (n_dim == 3) then
       !$omp do
       do n = 1, n_photons
          rr = prng%rngs(proc_id)%unif_01()
          dist = LT_get_col(tbl, 1, rr)
          xyz_out(:, n) =  xyz_in(:, n) + prng%rngs(proc_id)%sphere(dist)
       end do
       !$omp end do
    else
       print *, "phmc_do_absorption: unknown n_dim", n_dim
       stop
    end if
    !$omp end parallel

  end subroutine phmc_do_absorption

  ! !> Given a list of photon production positions (xyz_in), compute where they
  ! !> end up (xyz_out).
  ! subroutine phe_mc_get_endloc(xyz_in, xyz_out, n_dim, n_photons, prng)
  !   use m_lookup_table
  !   use m_random
  !   use omp_lib
  !   integer, intent(in)        :: n_photons !< Number of photons
  !   !> Input (x,y,z) values
  !   real(dp), intent(in)       :: xyz_in(3, n_photons)
  !   !> Output (x,y,z) values
  !   real(dp), intent(out)      :: xyz_out(3, n_photons)
  !   integer, intent(in)        :: n_dim     !< 2 or 3 dimensional
  !   type(PRNG_t), intent(inout) :: prng       !< Random number generator
  !   integer                    :: n, proc_id
  !   real(dp)                   :: rr, dist

  !   !$omp parallel private(n, rr, dist, proc_id)
  !   proc_id = 1+omp_get_thread_num()

  !   if (n_dim == 2) then
  !      !$omp do
  !      do n = 1, n_photons
  !         rr = prng%rngs(proc_id)%unif_01()
  !         dist = 0.00001_dp
  !         ! Pick a random point on a sphere, and ignore the last dimension
  !         xyz_out(1:3, n) =  xyz_in(1:3, n) + &
  !              prng%rngs(proc_id)%sphere(dist)
  !      end do
  !      !$omp end do
  !   else if (n_dim == 3) then
  !      !$omp do
  !      do n = 1, n_photons
  !         rr = prng%rngs(proc_id)%unif_01()
  !         dist = 0.00001_dp
  !         xyz_out(:, n) =  xyz_in(:, n) + prng%rngs(proc_id)%sphere(dist)
  !      end do
  !      !$omp end do
  !   else
  !      print *, "phmc_do_absorption: unknown n_dim", n_dim
  !      stop
  !   end if
  !   !$omp end parallel

  ! end subroutine phe_mc_get_endloc

  !> Set the source term due to photoionization. At most phmc_num_photons
  !> discrete photons are produced.
  subroutine phmc_set_src(tree, rng, i_src, i_photo, use_cyl, dt)
    use m_random
    use m_lookup_table
    use m_dielectric
    use omp_lib

    type(af_t), intent(inout)  :: tree   !< Tree
    type(RNG_t), intent(inout) :: rng    !< Random number generator
    !> Input variable that contains photon production per cell
    integer, intent(in)        :: i_src
    !> Input variable that contains photoionization source rate
    integer, intent(in)        :: i_photo
    !> Use cylindrical coordinates (only in 2D)
    logical, intent(in)        :: use_cyl
    !> Time step, if present use "physical" photons
    real(dp), intent(in), optional :: dt

    integer                     :: lvl, id, nc, min_lvl
    integer                     :: IJK, n, n_used
    integer                     :: proc_id, n_procs
    integer                     :: pho_lvl, max_num_photons
    real(dp)                    :: dr(NDIM), dt_fac, dist, n_produced
    real(dp)                    :: sum_production_rate, pi_lengthscale
    real(dp), allocatable       :: xyz_src(:, :)
    real(dp), allocatable       :: xyz_abs(:, :)
#if NDIM == 2
    real(dp)                    :: tmp, r(3)
    real(dp), parameter         :: pi = acos(-1.0_dp)
#endif
    type(PRNG_t)                :: prng
    type(af_loc_t), allocatable :: ph_loc(:)
    real(dp), parameter         :: small_value = 1e-100_dp

    if (NDIM == 3 .and. use_cyl) error stop "phmc_set_src: use_cyl is .true."

    if (ST_use_dielectric) then
       call dielectric_reset_photons()
    end if

    nc = tree%n_cell

    ! Initialize parallel random numbers
    n_procs = omp_get_max_threads()
    call prng%init_parallel(n_procs, rng)

    ! Compute the sum of the photon production rate
    call af_tree_sum_cc(tree, i_src, sum_production_rate)

    ! dt_fac is used to convert a discrete number of photons to a number of
    ! photons per unit time. Depending on the arguments dt and phmc_num_photons,
    ! 'physical' photons with a weight of phmc_min_weight, or super-photons with
    ! a weight larger (or smaller) than one can be used.
    if (present(dt)) then
       ! Create "physical" photons when less than phmc_num_photons are produced,
       ! otherwise create approximately phmc_num_photons
       n_produced = dt * sum_production_rate / phmc_min_weight

       if (n_produced < phmc_num_photons) then
          dt_fac = dt/phmc_min_weight
       else
          dt_fac = phmc_num_photons / (sum_production_rate + small_value)
       end if
    else
       ! Create approximately phmc_num_photons by setting dt_fac like this
       dt_fac = phmc_num_photons / (sum_production_rate + small_value)
    end if

    ! Over-estimate number of photons because of stochastic production
    max_num_photons = nint(1.2_dp * dt_fac * sum_production_rate + 1000)

    call phmc_generate_photons(tree, dt_fac, i_src, xyz_src, &
         max_num_photons, n_used, prng)

    allocate(xyz_abs(3, n_used))

    if (use_cyl) then           ! 2D only
       ! Get location of absorption. On input, xyz is set to (r, z, 0). On
       ! output, the coordinates thus correspond to (x, z, y)
       call phmc_do_absorption(xyz_src, xyz_abs, 3, n_used, phmc_tbl%tbl, prng)

       !$omp do
       do n = 1, n_used
          ! Convert back to (r,z) coordinates
          xyz_abs(1, n) = sqrt(xyz_abs(1, n)**2 + xyz_abs(3, n)**2)
          xyz_abs(2, n) = xyz_abs(2, n)
       end do
       !$omp end do

       if (ST_use_dielectric) then
          ! Handle photons that collide with dielectrics separately
          call dielectric_photon_absorption(tree, xyz_src, xyz_abs, &
               2, n_used, 1.0_dp/dt_fac)
       end if
    else
       ! Get location of absorbption
       call phmc_do_absorption(xyz_src, xyz_abs, NDIM, n_used, phmc_tbl%tbl, prng)

       if (ST_use_dielectric) then
          ! Handle photons that collide with dielectrics separately
          call dielectric_photon_absorption(tree, xyz_src, xyz_abs, &
               NDIM, n_used, 1.0_dp/dt_fac)
       end if
    end if

    allocate(ph_loc(n_used))

    if (phmc_const_dx) then
       ! Get a typical length scale for the absorption of photons
       pi_lengthscale = LT_get_col(phmc_tbl%tbl, 1, phmc_absorp_fac)

       ! Determine at which level we estimate the photoionization source term. This
       ! depends on the typical length scale for absorption.
       pho_lvl = get_lvl_length(maxval(tree%dr_base), pi_lengthscale)

       !$omp parallel do
       do n = 1, n_used
          ph_loc(n) = af_get_loc(tree, xyz_abs(1:NDIM, n), pho_lvl)
       end do
       !$omp end parallel do
    else
       !$omp parallel private(n, dist, lvl, proc_id)
       proc_id = 1+omp_get_thread_num()
       !$omp do
       do n = 1, n_used
          dist = phmc_absorp_fac * norm2(xyz_abs(1:NDIM, n) - xyz_src(1:NDIM, n))
          if (dist < phmc_min_dx) dist = phmc_min_dx
          lvl = get_rlvl_length(maxval(tree%dr_base), dist, prng%rngs(proc_id))
          ph_loc(n) = af_get_loc(tree, xyz_abs(1:NDIM, n), lvl)
       end do
       !$omp end do
       !$omp end parallel
    end if

    ! Clear variable i_photo, in which we will store the photoionization source term
    call af_tree_clear_cc(tree, i_photo)

    do n = 1, n_used
       id = ph_loc(n)%id
       if (id > af_no_box) then
#if NDIM == 1
          i = ph_loc(n)%ix(1)
          dr = tree%boxes(id)%dr

          tree%boxes(id)%cc(IJK, i_photo) = &
                  tree%boxes(id)%cc(IJK, i_photo) + &
                  phmc_tbl%frac_in_tbl/(dt_fac * product(dr))
#elif NDIM == 2
          i = ph_loc(n)%ix(1)
          j = ph_loc(n)%ix(2)
          dr = tree%boxes(id)%dr

          if (use_cyl) then
             tmp = dt_fac * 2 * pi
             r(1:2) = af_r_cc(tree%boxes(id), [i, j])
             tree%boxes(id)%cc(i, j, i_photo) = &
                  tree%boxes(id)%cc(i, j, i_photo) + &
                  phmc_tbl%frac_in_tbl/(tmp * product(dr) * r(1))
          else
             tree%boxes(id)%cc(IJK, i_photo) = &
                  tree%boxes(id)%cc(IJK, i_photo) + &
                  phmc_tbl%frac_in_tbl/(dt_fac * product(dr))
          end if
#elif NDIM == 3
          i = ph_loc(n)%ix(1)
          j = ph_loc(n)%ix(2)
          k = ph_loc(n)%ix(3)
          dr = tree%boxes(id)%dr

          tree%boxes(id)%cc(IJK, i_photo) = &
               tree%boxes(id)%cc(IJK, i_photo) + &
               phmc_tbl%frac_in_tbl/(dt_fac * product(dr))
#endif
       end if
    end do

    ! Set ghost cells on highest level with photon source
    if (phmc_const_dx) then
       min_lvl = pho_lvl
    else
       min_lvl = 1
    end if

    !$omp parallel private(lvl, i, id)
    ! Prolong to finer grids
    do lvl = min_lvl, tree%highest_lvl-1
       !$omp do
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call af_gc_box(tree, id, [i_photo])
       end do
       !$omp end do

       !$omp do
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call af_prolong_linear_from(tree%boxes, id, i_photo, add=.true.)
       end do
       !$omp end do
    end do
    !$omp end parallel

    call prng%update_seed(rng)
  end subroutine phmc_set_src

!   !> Whole photoemission procedure
!   subroutine phe_mc_set_src(tree, rng, i_src, i_photo, use_cyl, dt)
!     use m_random
!     use m_lookup_table
!     use m_dielectric2
!     use omp_lib

!     type(af_t), intent(inout)  :: tree   !< Tree
!     type(RNG_t), intent(inout) :: rng    !< Random number generator
!     !> Input variable that contains photon production per cell
!     integer, intent(in)        :: i_src
!     !> Input variable that contains photoionization source rate
!     integer, intent(in)        :: i_photo
!     !> Use cylindrical coordinates (only in 2D)
!     logical, intent(in)        :: use_cyl
!     !> Time step, if present use "physical" photons
!     real(dp), intent(in), optional :: dt

!     integer                     :: nc
!     integer                     :: n, n_used
!     integer                     :: proc_id, n_procs
!     integer                     :: pho_lvl, max_num_photons
!     real(dp)                    :: tmp, dt_fac, r(3)
!     real(dp)                    :: sum_production_rate, pi_lengthscale
!     real(dp), allocatable       :: xyz_src(:, :)
!     real(dp), allocatable       :: xyz_abs(:, :)
! #if NDIM == 2
!     real(dp), parameter         :: pi = acos(-1.0_dp)
! #endif
!     type(PRNG_t)                :: prng
!     type(af_loc_t), allocatable :: ph_loc(:)
!     real(dp), parameter         :: small_value = 1e-100_dp

!     if (NDIM == 3 .and. use_cyl) error stop "phmc_set_src: use_cyl is .true."

!     if (ST_use_dielectric) then
!        call dielectric2_reset_photons()
!     end if

!     nc = tree%n_cell

!     ! Initialize parallel random numbers
!     n_procs = omp_get_max_threads()
!     call prng%init_parallel(n_procs, rng)

!     ! Compute the sum of the photon production rate
!     call af_tree_sum_cc(tree, i_src, sum_production_rate)

!     ! dt_fac is used to convert a discrete number of photons to a number of
!     ! photons per unit time. Depending on the arguments dt and phmc_num_photons,
!     ! 'physical' photons with a weight of 1.0, or super-photons with a weight
!     ! larger (or smaller) than one can be used.
!     if (present(dt)) then
!        ! Create "physical" photons when less than phmc_num_photons are produced,
!        ! otherwise create approximately phmc_num_photons
!        dt_fac = min(dt, phmc_num_photons / (sum_production_rate + small_value))
!     else
!        ! Create approximately phmc_num_photons by setting dt_fac like this
!        dt_fac = phmc_num_photons / (sum_production_rate + small_value)
!     end if

!     ! Over-estimate number of photons because of stochastic production
!     max_num_photons = nint(1.2_dp * dt_fac * sum_production_rate + 1000)

!     call phmc_generate_photons(tree, dt_fac, i_src, xyz_src, &
!          max_num_photons, n_used, prng)

!     allocate(xyz_abs(3, n_used))

!     if (use_cyl) then           ! 2D only
!        ! Get location of absorption. On input, xyz is set to (r, z, 0). On
!        ! output, the coordinates thus correspond to (x, z, y)
!        call phe_mc_get_endloc(xyz_src, xyz_abs, 3, n_used, prng)

!        !$omp do
!        do n = 1, n_used
!           ! Convert back to (r,z) coordinates
!           xyz_abs(1, n) = sqrt(xyz_abs(1, n)**2 + xyz_abs(3, n)**2)
!           xyz_abs(2, n) = xyz_abs(2, n)
!           xyz_src(2, n) = xyz_src(3, n)
!        end do
!        !$omp end do

!        if (ST_use_dielectric) then
!           ! Handle photons that collide with dielectrics separately
!           call dielectric2_photon_absorption(tree, xyz_src, xyz_abs, &
!                2, n_used, 1.0_dp/dt_fac)
!        end if
!     else
!        ! Get location of absorbption
!        call phe_mc_get_endloc(xyz_src, xyz_abs, NDIM, n_used, prng)

!        if (ST_use_dielectric) then
!           ! Handle photons that collide with dielectrics separately
!           call dielectric2_photon_absorption(tree, xyz_src, xyz_abs, &
!                NDIM, n_used, 1.0_dp/dt_fac)
!        end if
!     end if

!     call prng%update_seed(rng)
!   end subroutine phe_mc_set_src


  subroutine phmc_generate_photons(tree, dt_fac, i_src, xyz_src, n_max, n_used, prng)
    use omp_lib
    use m_random
    type(af_t), intent(in)               :: tree
    !> Time step to convert photon production rate to number of photons
    real(dp), intent(in)                 :: dt_fac
    !> Index of cell-centered variable containing photon production rate
    integer, intent(in)                  :: i_src
    !> Origin of the produced photons
    real(dp), allocatable, intent(inout) :: xyz_src(:, :)
    !> Maximum number of photons that will be produced
    integer, intent(in)                  :: n_max
    !> Number of photons that have been created
    integer, intent(out)                 :: n_used
    !> Parallel random number generator
    type(PRNG_t), intent(inout)          :: prng

    integer               :: proc_id, n_procs
    integer               :: lvl, ix, id, IJK, n, nc, i_ph
    integer               :: n_create, my_num_photons
    real(dp)              :: tmp, dr(NDIM), r(3)
    integer, allocatable  :: photons_per_thread(:), photon_thread(:)
    integer, allocatable  :: ix_offset(:)
    real(dp), allocatable :: xyz_tmp(:, :)

    nc      = tree%n_cell
    n_procs = omp_get_max_threads()

    allocate(xyz_src(3, n_max))
    allocate(photons_per_thread(n_procs))
    allocate(photon_thread(n_max))

    ! Loop over all leaves and create photons using random numbers
    n_used = 0

    !$omp parallel private(lvl, ix, id, IJK, n, r, dr, i_ph, proc_id, &
    !$omp tmp, n_create, my_num_photons)
    proc_id = 1+omp_get_thread_num()
    my_num_photons = 0

    do lvl = 1, tree%highest_lvl
       dr = af_lvl_dr(tree, lvl)
       !$omp do
       do ix = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(ix)

          do KJI_DO(1,nc)

#if NDIM == 2
             if (tree%boxes(id)%coord_t == af_cyl) then
                tmp = af_cyl_radius_cc(tree%boxes(id), i)
                tmp = dt_fac * 2 * UC_pi * tmp * &
                     tree%boxes(id)%cc(i, j, i_src) * product(dr)
             else
                tmp = dt_fac * tree%boxes(id)%cc(i, j, i_src) * product(dr)
             end if
#else
             tmp = dt_fac * tree%boxes(id)%cc(IJK, i_src) * product(dr)
#endif

             n_create = floor(tmp)

             if (prng%rngs(proc_id)%unif_01() < tmp - n_create) &
                  n_create = n_create + 1

             if (n_create > 0) then
                !$omp critical
                i_ph = n_used
                n_used = n_used + n_create
                !$omp end critical
                my_num_photons = my_num_photons + n_create

                do n = 1, n_create
                   ! Location of production randomly chosen in cell
                   r(1)   = prng%rngs(proc_id)%unif_01()
                   r(2)   = prng%rngs(proc_id)%unif_01()
#if NDIM == 2
                   xyz_src(1:NDIM, i_ph+n) = af_rr_cc(tree%boxes(id), &
                        [IJK] - 0.5_dp + r(1:NDIM))
                   xyz_src(3, i_ph+n) = 0.0_dp
#elif NDIM == 3
                   r(3)   = prng%rngs(proc_id)%unif_01()
                   xyz_src(:, i_ph+n) = af_rr_cc(tree%boxes(id), &
                        [IJK] - 0.5_dp + r(1:NDIM))
#endif
                   photon_thread(i_ph+n) = proc_id
                end do
             end if
          end do; CLOSE_DO
       end do
       !$omp end do nowait
    end do

    photons_per_thread(proc_id) = my_num_photons
    !$omp end parallel

    ! Sort the xyz_src array so that runs are deterministic (the order in which
    ! the threads write is not deterministic)
    allocate(ix_offset(n_procs))
    allocate(xyz_tmp(3, n_used))

    ix_offset(1) = 0
    do n = 2, n_procs
       ix_offset(n) = ix_offset(n-1) + photons_per_thread(n-1)
    end do

    photons_per_thread = 0
    do n = 1, n_used
       i = photon_thread(n)
       photons_per_thread(i) = photons_per_thread(i) + 1
       ix = ix_offset(i) + photons_per_thread(i)
       xyz_tmp(:, ix) = xyz_src(:, n)
    end do
    xyz_src(:, 1:n_used) = xyz_tmp(:, 1:n_used)

  end subroutine phmc_generate_photons

end module m_photoi_mc
