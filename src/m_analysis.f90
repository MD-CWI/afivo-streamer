#include "../afivo/src/cpp_macros.h"
!> Module with routines to help analyze simulations
module m_analysis
  use m_af_all
  use m_types

  implicit none
  private

  ! Public methods
  public :: analysis_get_maxima

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

  !> Get the conductivity of an axisymmetric streamer at a z-coordinate
  subroutine analysis_get_sigma_cyl(tree, rmax, z, sigma)
    ! Include some modules
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: rmax  !< Integrate up to this radius
    real(dp), intent(in)   :: z     !< z-coordinate
    real(dp), intent(out)  :: sigma !< The conductivity

    real(dp) :: ne_fld(2)
    logical  :: success
    integer  :: id_guess

    id_guess = -1
    sigma    = 0.0_dp

    ! Interpolate up to rmax
    do ...
       r = ...
       ne_fld = af_interp1(tree, [r, z], [i_electron, i_electric_fld], &
            success, id_guess)
       if (.not. success) error stop "oh oh"

       ! Get mobility (see e.g. m_fluid_lfa)
       Td = ... ! Townsend
       ! N_inv is inverse gas number density
       mu = LT_get_col(td_tbl, td_mobility, Td) * N_inv 

       ! Multiply with something like 2 pi r dr

       ! Update total
       sigma = sigma + ...
    end do

  end subroutine analysis_get_sigma_cyl

end module m_analysis
