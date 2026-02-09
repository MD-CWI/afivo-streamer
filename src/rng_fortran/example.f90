program example
  use m_random

  ! Define double precision real type
  integer, parameter :: dp = kind(0.0d0)

  ! Create a random number generator
  type(RNG_t) :: rng

  ! Set the initial seed for the generator
  call rng%set_random_seed()

  ! Print some random numbers
  print *, "Uniform random number: ", rng%unif_01()
  print *, "8-byte random integer: ", rng%int_64()
  print *, "4-byte random integer: ", rng%int_32()
  print *, "Two normal numbers:    ", rng%two_normals()
  print *, "Poisson(10.0) deviate: ", rng%poisson(10.0_dp)
  print *, "Point on unit circle:  ", rng%circle(1.0_dp)
  print *, "Point on unit sphere:  ", rng%sphere(1.0_dp)
end program example
