program example
  use m_random

  ! Define 8-byte integer type
  integer, parameter :: i8 = selected_int_kind(18)

  ! Define double precision real type
  integer, parameter :: dp = kind(0.0d0)

  ! Create a random number generator
  type(RNG_t) :: rng

  ! Set the initial seed for the generator
  call rng%set_seed([89732_i8, 1892342989_i8])

  ! Print some random numbers
  print *, "Uniform random number: ", rng%unif_01()
  print *, "8-byte random integer: ", rng%int_8()
  print *, "4-byte random integer: ", rng%int_4()
  print *, "Two normal numbers:    ", rng%two_normals()
  print *, "Poisson(10.0) deviate: ", rng%poisson(10.0_dp)
  print *, "Point on unit circle:  ", rng%circle(1.0_dp)
  print *, "Point on unit sphere:  ", rng%sphere(1.0_dp)
end program example
