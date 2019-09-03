#include "../afivo/src/cpp_macros.h"
!> Module with routines to help analyze simulations
module m_analysis
  use m_af_all
  use m_types

  implicit none
  private

  ! Public methods
  public :: analysis_get_maxima
  public :: analysis_get_sigma
  public :: analysis_get_cross

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
#if NDIM == 2
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

  subroutine sigma_calculator(boxes, id, nc)
    use m_lookup_table
    use m_units_constants
    use m_gas
    use m_transport_data
    use m_streamer
    type(box_t), intent(inout) :: boxes(:)
    integer, intent(in)        :: id, nc

    real(dp) :: ne_fld(2), mu, Td, N_inv
    integer  :: n, m

    N_inv = 1/gas_number_density

#if NDIM == 2
    do n = 1, nc
      do m = 1, nc
        ne_fld = boxes(id)%cc(n, m, [i_electron, i_electric_fld])
        Td = ne_fld(2) * SI_to_Townsend * N_inv
        mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv
        boxes(id)%cc(n , m, i_conductivity) = mu * ne_fld(1) * UC_elec_charge
      end do
    end do
#endif

  end subroutine sigma_calculator

  subroutine analysis_get_sigma(tree)
    type(af_t), intent(inout) :: tree

    integer     :: lvl, id, i, nc

    nc = tree%n_cell

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
      !$omp do
      do i = 1, size(tree%lvls(lvl)%leaves)
        id = tree%lvls(lvl)%leaves(i)
        call sigma_calculator(tree%boxes, id, nc)
      end do
      !$omp end do
    end do
    !$omp end parallel
  end subroutine analysis_get_sigma

  !> Get the conductivity and densities of an axisymmetric streamer at a z-coordinate
  subroutine analysis_get_cross(tree, rmax, z, sigma, elec_dens, charge_dens, current_dens)
    use m_lookup_table
    use m_transport_data
    use m_gas
    use m_streamer
    use m_units_constants
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: rmax         !< Integrate up to this radius
    real(dp), intent(in)   :: z            !< z-coordinate
    real(dp), intent(out)  :: sigma        !< The conductivity
    real(dp), intent(out)  :: elec_dens    !< electron density
    real(dp), intent(out)  :: charge_dens  !< charge density
    real(dp), intent(out)  :: current_dens !< current density

    real(dp) :: ne_fld_rhs(3), mu, Td, r, dr, N_inv
    real(dp) :: d_sigma, d_elec_dens, d_charge_dens, d_current_dens
    logical  :: success
    integer  :: id_guess, i, m

    id_guess     = -1
    sigma        = 0.0_dp
    elec_dens    = 0.0_dp
    charge_dens  = 0.0_dp
    current_dens = 0.0_dp
    N_inv = 1/gas_number_density
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
       d_sigma = mu * ne_fld_rhs(1) * 2 * UC_pi * r * dr
       d_elec_dens = ne_fld_rhs(1) * 2 * UC_pi * r * dr
       d_charge_dens = ne_fld_rhs(3) * UC_eps0 * 2 * UC_pi * r * dr / UC_elec_charge
       d_current_dens = ne_fld_rhs(2) * mu * ne_fld_rhs(1) * 2 * UC_pi * r * dr

       ! Update total
       sigma = sigma + d_sigma
       elec_dens = elec_dens + d_elec_dens
       charge_dens = charge_dens + d_charge_dens
       current_dens = current_dens + d_current_dens
    end do

  end subroutine analysis_get_cross

end module m_analysis
