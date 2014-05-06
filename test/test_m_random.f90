program test_m_random
   use m_random

   implicit none
   integer, parameter  :: dp = kind(0.0d0)
   integer, parameter  :: n_samples = 10*1000*1000
   integer             :: nn, rng_seed
   real(dp)            :: rand_value, sum_sq, new_mean, old_mean, variance
   real(dp), parameter :: pi        = acos(-1.0_dp)
   real(dp)            :: time_start, time_end

   print *, "Testing implementation of m_random.f90"
   print *, "This is just checking whether everything works, and by no means"
   print *, "a test of the 'randomness' of the pseudo random number generator."
   print *, "For these tests, ", n_samples, " values are used"

   call system_clock(count=rng_seed)
   call RNG_set_seed(rng_seed)

   sum_sq = 0.0_dp
   old_mean = 0.0_dp
   call cpu_time(time_start)
   do nn = 1, n_samples
      call random_number(rand_value)
      new_mean = old_mean + (rand_value - old_mean) / nn
      sum_sq = sum_sq + (rand_value - old_mean) * (rand_value - new_mean)
      old_mean = new_mean
   end do
   call cpu_time(time_end)
   variance = sum_sq / n_samples

   print *, ""
   print *, "For uniform random numbers, the result is:"
   print *, "nanoseconds per number (upper bound)", 1.0e9_dp * (time_end - time_start) / n_samples
   print *, "mean/<mean>", new_mean/0.5_dp
   print *, "std dev/<std dev>", sqrt(variance)*sqrt(12.0_dp)

   sum_sq = 0.0_dp
   old_mean = 0.0_dp
   call cpu_time(time_start)
   do nn = 1, n_samples
      rand_value = 1 + RNG_normal()
      new_mean = old_mean + (rand_value - old_mean) / nn
      sum_sq = sum_sq + (rand_value - old_mean) * (rand_value - new_mean)
      old_mean = new_mean
   end do
   call cpu_time(time_end)
   variance = sum_sq / n_samples

   print *, ""
   print *, "For normal/Gaussian random numbers, the result is:"
   print *, "nanoseconds per number (upper bound)", 1.0e9_dp * (time_end - time_start) / n_samples
   print *, "mean/<mean>", new_mean/1.0_dp ! Above we add one to RNG_normal()
   print *, "std dev/<std dev>", sqrt(variance)

   sum_sq = 0.0_dp
   old_mean = 0.0_dp
   call cpu_time(time_start)
   do nn = 1, n_samples
      rand_value = RNG_poisson(pi)
      new_mean = old_mean + (rand_value - old_mean) / nn
      sum_sq = sum_sq + (rand_value - old_mean) * (rand_value - new_mean)
      old_mean = new_mean
   end do
   call cpu_time(time_end)
   variance = sum_sq / n_samples

   print *, ""
   print *, "For Poisson random numbers, the result is:"
   print *, "nanoseconds per number (upper bound)", 1.0e9_dp * (time_end - time_start) / n_samples
   print *, "mean/<mean>", new_mean/pi ! Above we add one to RNG_normal()
   print *, "std dev/<std dev>", sqrt(variance/pi)

end program test_m_random
