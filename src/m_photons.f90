module m_photons

  implicit none
  private

  integer, parameter                     :: dp = kind(0.0d0)

  ! Public methods
  public :: PH_get_tbl_air
  public :: PH_do_absorp

contains

  function PH_get_tbl_air(p_O2) result(tbl)
    use m_lookup_table
    real(dp), intent(in) :: p_O2     !< Partial pressure of oxygen (bar)
    type(LT_table_t)     :: tbl      !< The lookup table

    integer, parameter   :: tbl_size  = 1000
    integer              :: n
    real(dp), parameter  :: huge_dist = 1e9_dp, sixth = 1/6.0_dp
    real(dp), allocatable :: fsum(:), dist(:)
    real(dp)             :: dr, r_prev, incr

    allocate(fsum(tbl_size))
    allocate(dist(tbl_size))
    tbl     = LT_create(0.0_dp, 1.0_dp, tbl_size, 1)
    dr      = 1e-10_dp
    dist(1) = 0
    fsum(1) = 0
    r_prev = epsilon(1.0_dp)

    ! Use Simpson's rule for integration
    do n = 2, tbl_size
       incr = dr * sixth * (absfunc_air(r_prev, p_O2) + &
            4 * absfunc_air(r_prev+0.5_dp*dr, p_O2) + &
            absfunc_air(r_prev+dr, p_O2))

       do while (incr > 2.0_dp / tbl_size)
          dr = dr * 0.5_dp
          incr = dr * sixth * (absfunc_air(r_prev, p_O2) + &
               4 * absfunc_air(r_prev+0.5_dp*dr, p_O2) + &
               absfunc_air(r_prev+dr, p_O2))
       end do

       do while (incr < 1.0_dp / (tbl_size-1) .and. dr < huge_dist)
          dr = dr * 2
          incr = dr * sixth * (absfunc_air(r_prev, p_O2) + &
               4 * absfunc_air(r_prev+0.5_dp*dr, p_O2) + &
               absfunc_air(r_prev+dr, p_O2))
       end do

       dist(n) = r_prev + dr
       fsum(n) = fsum(n-1) + incr

       ! The integral will not be completely accurate, and small errors can lead
       ! to big problems (dr -> infinity). So, stop when dr gets too large.
       if ((n > 2 .and. dr > 10 * r_prev) .or. fsum(n) > 1)  then
          fsum(n:) = 1
          dist(n:) = dist(n-1)
          exit
       end if

       r_prev = r_prev + dr
    end do

    call LT_set_col(tbl, 1, fsum, dist)
  end function PH_get_tbl_air

  real(dp) function absfunc_air(dist, p_O2)
    use m_units_constants
    real(dp), intent(in) :: dist, p_O2
    real(dp), parameter  :: chi_min  = 3.5_dp / UC_torr_to_bar
    real(dp), parameter  :: chi_max  = 200 / UC_torr_to_bar

    absfunc_air = (exp(-chi_min * p_O2 * dist) - &
         exp(-chi_max * p_O2 * dist)) / (dist * log(chi_max/chi_min))
  end function absfunc_air

  subroutine PH_do_absorp(xyz_in, xyz_out, n_photons, tbl, rng)
    use m_lookup_table
    use m_random
    integer, intent(in)          :: n_photons
    real(dp), intent(in)         :: xyz_in(3, n_photons)
    real(dp), intent(out)        :: xyz_out(3, n_photons)
    type(LT_table_t), intent(in) :: tbl
    type(RNG_t), intent(inout)   :: rng
    integer                      :: n
    real(dp)                     :: rr, dist

    do n = 1, n_photons
       rr = rng%uni_01()
       dist = LT_get_col(tbl, 1, rr)
       xyz_out(:, n) =  xyz_in(:, n) + rng%sphere(dist)
    end do
  end subroutine PH_do_absorp

  ! real(dp), parameter :: p_quench = 30.0D0 * UC_torr_to_bar
  ! Compute quench factor, because some excited species will be quenched by
  ! collisions, preventing the emission of a UV photon
  ! quench_fac = p_quench / (p_bar + p_quench)

end module m_photons
