program test_parallel
  use m_random
  use omp_lib

  implicit none
  integer, parameter    :: dp        = kind(0.0d0)
  integer, parameter    :: i8        = selected_int_kind(18)
  integer, parameter    :: n_samples = 100*1000*1000
  integer               :: n, tid
  integer               :: time_start, time_end, count_rate
  real(dp)              :: mean, variance
  real(dp), allocatable :: rand_results(:)
  type(RNG_t)           :: rng
  type(PRNG_t)          :: prng

  print *, "Testing parallel random number generation"
  print *, "Number of threads", omp_get_max_threads()

  allocate(rand_results(n_samples))

  call rng%set_seed([89732_i8, 1892342989_i8])
  call prng%init_parallel(omp_get_max_threads(), rng)
  call system_clock(time_start, count_rate)

  !$omp parallel private(n, tid, rng)
  tid = omp_get_thread_num() + 1
  !$omp do
  do n = 1, n_samples
     rand_results(n) = prng%rngs(tid)%unif_01()
  end do
  !$omp end do
  !$omp end parallel

  ! Update the rng seed afterwards, using the 'evolved' prng state
  call prng%update_seed(rng)

  call system_clock(time_end)

  mean     = sum(rand_results) / n_samples
  variance = sum((rand_results - mean)**2) / n_samples

  print *, ""
  print *, "For uniform random numbers (unif_01), the result is:"
  print *, "nanoseconds per number (upper bound)", &
       (1.0e9_dp/count_rate) * (time_end - time_start) / n_samples
  print *, "mean/<mean>", mean/0.5_dp
  print *, "std dev/<std dev>", sqrt(variance)*sqrt(12.0_dp)

end program test_parallel
