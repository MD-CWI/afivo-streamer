# rng_fortran

A module for generating (pseudo) random numbers. Internally, the
xoroshiro128plus generator is used. The following types of random numbers are
currently supported:

* 8 byte random integers
* 4 byte random integers
* (0,1] uniform random numbers
* Normal random numbers (in pairs of two)
* Poisson-distributed random numbers
* Random points on a circle
* Random points on a sphere

A usage example is given in `example.f90`.
