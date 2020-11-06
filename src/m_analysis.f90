#include "../afivo/src/cpp_macros.h"
!> Module with routines to help analyze simulations
module m_analysis
  use m_af_all
  use m_types

  implicit none
  private

  ! Public methods
  public :: analysis_get_maxima
  public :: analysis_get_cross
  public :: analysis_zmin_zmax_threshold
  public :: analysis_max_var_region
  public :: analysis_max_var_product

contains

  !> Find at most n_max maxima of a variable. Maxima are determined by looking
  !> at the direct neighbors.
  subroutine analysis_get_maxima(tree, iv, threshold, n_max, &
       coord_val, n_found)
    type(af_t), intent(in)  :: tree
    integer, intent(in)     :: iv                       !< Index of variable
    real(dp), intent(in)    :: threshold                !< Consider maxima above this value
    integer, intent(in)     :: n_max                    !< Maximum number of maxima
    real(dp), intent(inout) :: coord_val(NDIM+1, n_max) !< Up to n_max maxima and their coordinates
    integer, intent(out)    :: n_found                  !< Number of maxima found

    integer  :: IJK, nc, n, lvl, id
    real(dp) :: val, neighbs(2*NDIM)

    n_found = 0
    nc      = tree%n_cell

    !$omp parallel private(lvl, n, id, IJK, val, neighbs)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do n = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(n)

          associate (cc => tree%boxes(id)%cc)
            do KJI_DO(1,nc)
               val = cc(IJK, iv)
#if NDIM == 1
               neighbs = [cc(i-1, iv), cc(i+1, iv)]
#elif NDIM == 2
               neighbs = [cc(i-1, j, iv), cc(i+1, j, iv), &
                    cc(i, j-1, iv), cc(i, j+1, iv)]
#elif NDIM == 3
               neighbs = [cc(i-1, j, k, iv), cc(i+1, j, k, iv), &
                    cc(i, j-1, k, iv), cc(i, j+1, k, iv), &
                    cc(i, j, k-1, iv), cc(i, j, k+1, iv)]
#endif
               ! The value has to be strictly larger than at least one of the
               ! neigbors, and not smaller than the neighbors
               if (val > threshold .and. all(val >= neighbs) &
                    .and. any(val > neighbs)) then
                  ! New maximum found
                  !$omp critical
                  n_found = n_found + 1
                  if (n_found <= n_max) then
                     coord_val(NDIM+1, n_found) = val
                     coord_val(1:NDIM, n_found) = af_r_cc(tree%boxes(id), [IJK])
                  end if
                  !$omp end critical
               end if
            end do; CLOSE_DO
          end associate

       end do
       !$omp end do
    end do
    !$omp end parallel

  end subroutine analysis_get_maxima

  !> Find minimum and maximum z coordinate where a variable exceeds a threshold
  subroutine analysis_zmin_zmax_threshold(tree, iv, threshold, limits, z_minmax)
    type(af_t), intent(in) :: tree
    integer, intent(in)    :: iv !< Index of variable
    !> Threshold for variable
    real(dp), intent(in)   :: threshold
    !> Limits for min/max
    real(dp), intent(in)   :: limits(2)
    !> Minimum/maximum z coordinate above threshold
    real(dp), intent(out)  :: z_minmax(2)

    call af_reduction_vec(tree, box_minmax_z, reduce_minmax, &
         limits, z_minmax, 2)

  contains

    ! Find cell with min/max z coordinate that has a density exceeding a
    ! threshold
    function box_minmax_z(box, n_vals) result(vec)
      type(box_t), intent(in) :: box
      integer, intent(in)     :: n_vals
      real(dp)                :: vec(n_vals)

      integer  :: i, j, n, nc, ix(NDIM)
      logical  :: above
      real(dp) :: r(NDIM)

      nc = box%n_cell
      i = -1
      j = -1

      do n = 1, nc
#if NDIM == 1
         above = box%cc(n, iv) > threshold
#elif NDIM == 2
         above = maxval(box%cc(1:nc, n, iv)) > threshold
#elif NDIM == 3
         above = maxval(box%cc(1:nc, 1:nc, n, iv)) > threshold
#endif
         if (above) then
            if (i == -1) i = n
            j = n
         end if
      end do

      vec = [1e100_dp, -1e100_dp]
      ix = 1
      if (i /= -1) then
         ix(NDIM) = i
         r = af_r_cc(box, ix)
         vec(1) = r(NDIM)
      end if

      if (j /= -1) then
         ix(NDIM) = i
         r = af_r_cc(box, ix)
         vec(2) = r(NDIM)
      end if
    end function box_minmax_z

    !> Reduction method (e.g., min, max, sum)
    function reduce_minmax(vec_1, vec_2, n_vals) result(vec)
      integer, intent(in)  :: n_vals
      real(dp), intent(in) :: vec_1(n_vals), vec_2(n_vals)
      real(dp)             :: vec(n_vals)
      vec(1) = min(vec_1(1), vec_2(1))
      vec(2) = max(vec_1(2), vec_2(2))
    end function reduce_minmax

  end subroutine analysis_zmin_zmax_threshold

  !> Find maximal value for boxes that are (at least partially) in the rectangle
  !> from r0 to r1
  subroutine analysis_max_var_region(tree, iv, r0, r1, max_value, loc)
    type(af_t), intent(in)      :: tree
    integer, intent(in)         :: iv        !< Index of variable
    real(dp), intent(in)        :: r0(NDIM)  !< Minimum coordinates
    real(dp), intent(in)        :: r1(NDIM)  !< Maximum coordinates
    real(dp), intent(out)       :: max_value !< Found maximum
    type(af_loc_t), intent(out) :: loc

    call af_reduction_loc(tree, iv, box_max_region, reduce_max, &
         -1e100_dp, max_value, loc)

  contains

    subroutine box_max_region(box, iv, val, ix)
      type(box_t), intent(in) :: box
      integer, intent(in)     :: iv
      real(dp), intent(out)   :: val
      integer, intent(out)    :: ix(NDIM)
      integer                 :: nc
      real(dp)                :: r_max(NDIM)

      nc = box%n_cell
      r_max = box%r_min + box%n_cell * box%dr

      if (any(box%r_min > r1) .or. any(r_max < r0)) then
         ix = -1
         val = -1e100_dp
      else
         ix = maxloc(box%cc(DTIMES(1:nc), iv))
         val = box%cc(DINDEX(ix), iv)
      end if
    end subroutine box_max_region

  end subroutine analysis_max_var_region

  subroutine analysis_max_var_product(tree, ivs, max_value, loc)
    type(af_t), intent(in)      :: tree
    integer, intent(in)         :: ivs(:)    !< Indices of variables
    real(dp), intent(out)       :: max_value !< Found maximum
    type(af_loc_t), intent(out) :: loc

    call af_reduction_loc(tree, -1, box_max_product, reduce_max, &
         -1e100_dp, max_value, loc)

  contains

    subroutine box_max_product(box, iv, val, ix)
      type(box_t), intent(in) :: box
      integer, intent(in)     :: iv
      real(dp), intent(out)   :: val
      integer, intent(out)    :: ix(NDIM)
      integer                 :: nc

      nc  = box%n_cell
      ix  = maxloc(product(box%cc(DTIMES(1:nc), ivs), dim=NDIM+1))
      val = product(box%cc(DINDEX(ix), ivs))
    end subroutine box_max_product

  end subroutine analysis_max_var_product

  real(dp) function reduce_max(a, b)
    real(dp), intent(in) :: a, b
    reduce_max = max(a, b)
  end function reduce_max

  !> Get integrated quantities of an axisymmetric streamer at a z-coordinate
  subroutine analysis_get_cross(tree, rmax, z, elec_dens, charge_dens, current_dens)
    use m_lookup_table
    use m_transport_data
    use m_gas
    use m_streamer
    use m_units_constants
    use m_chemistry
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: rmax         !< Integrate up to this radius
    real(dp), intent(in)   :: z            !< z-coordinate
    real(dp), intent(out)  :: elec_dens    !< electron density
    real(dp), intent(out)  :: charge_dens  !< charge density
    real(dp), intent(out)  :: current_dens !< current density

    real(dp) :: ne_fld_rhs(3), mu, Td, r, dr, N_inv
    real(dp) :: d_elec_dens, d_charge_dens, d_current_dens
    logical  :: success
    integer  :: id_guess, i, m, n

    id_guess     = -1
    elec_dens    = 0.0_dp
    charge_dens  = 0.0_dp
    current_dens = 0.0_dp
    N_inv = 1.0_dp/gas_number_density
    dr = af_min_dr(tree)
    m = int(rmax/dr) + 1

    do i = 1, m
       r = i * rmax / (m + 1)
#if NDIM == 2
       ne_fld_rhs = af_interp1(tree, [r, z], [i_electron, i_electric_fld, i_rhs], &
            success, id_guess)
       if (.not. success) error stop "unsuccessful interp1"
#endif
       Td = ne_fld_rhs(2) * SI_to_Townsend * N_inv
       mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv       
       d_elec_dens = ne_fld_rhs(1) * 2.0_dp * UC_pi * r * dr
       d_charge_dens = ne_fld_rhs(3) * UC_eps0 * 2.0_dp * UC_pi * r * dr / UC_elec_charge
       d_current_dens = ne_fld_rhs(2) * mu * ne_fld_rhs(1) * 2.0_dp * UC_pi * r * dr * UC_elem_charge

       ! Update total
       elec_dens = elec_dens + d_elec_dens
       charge_dens = charge_dens + d_charge_dens
       current_dens = current_dens + d_current_dens
    end do

  end subroutine analysis_get_cross

end module m_analysis
