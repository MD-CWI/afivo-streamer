program test_m_photons
  use omp_lib
  use m_photons
  use m_random
  use m_lookup_table

  implicit none
  integer, parameter :: dp = kind(0.0d0)

  integer :: n_photons = 10*1000*1000
  real(dp), allocatable :: xyz_in(:, :), xyz_out(:, :)
  real(dp) :: t1, t2

  type(RNG_t) :: rng
  type(photoi_tbl_t) :: photoi_tbl

  allocate(xyz_in(3, n_photons))
  allocate(xyz_out(3, n_photons))

  call photoi_get_table_air(photoi_tbl, 0.2_dp, 10e-3_dp)
  xyz_in = 0

  t1 = omp_get_wtime()
  call photoi_do_absorption(xyz_in, xyz_out, 3, n_photons, photoi_tbl%tbl, rng)
  t2 = omp_get_wtime()

  print *, "time per photon", (t2-t1)/n_photons * 1e9_dp, "ns"
  print *, "avg dist", sum(norm2(xyz_out-xyz_in, 1))/n_photons
end program test_m_photons
