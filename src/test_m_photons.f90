program test_m_photons
  use m_photons
  use m_random
  use m_lookup_table

  implicit none
  integer, parameter :: dp = kind(0.0d0)

  integer :: n_photons = 10*1000*1000
  real(dp), allocatable :: xyz_in(:, :), xyz_out(:, :)
  real(dp) :: t1, t2

  type(RNG_t) :: rng
  type(LT_table_t) :: tbl

  allocate(xyz_in(3, n_photons))
  allocate(xyz_out(3, n_photons))

  tbl = PH_get_tbl_air(0.2_dp)
  xyz_in = 0

  call cpu_time(t1)
  call PH_do_absorp(xyz_in, xyz_out, n_photons, tbl, rng)
  call cpu_time(t2)

  print *, "time per photon", (t2-t1)/n_photons * 1e9_dp, "ns"
  print *, "avg dist", sum(norm2(xyz_out-xyz_in, 1))/n_photons
end program test_m_photons
