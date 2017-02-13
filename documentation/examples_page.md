# Examples {#examples-page}

The examples listed can be run from the `examples` folder. They are
automatically compiled after typing `make`.

# Basic

* @ref computational_domain_2d.f90 "2D", @ref computational_domain_3d.f90 "3D":
  How to define the computational domain
* @ref random_refinement_2d.f90 "2D", @ref random_refinement_3d.f90 "3D": Examples
  that show how to create an AMR tree, perform refinement at random, write
  output files, and how to fill ghost cells
* @ref boundary_conditions_2d.f90 "2D", @ref boundary_conditions_3d.f90 "3D":
  How to write a routine for boundary conditions

# Drift-diffusion

* @ref drift_diffusion_2d.f90 "2D", @ref drift_diffusion_3d.f90 "3D": A
  drift-diffusion example

# Multigrid

* @ref poisson_basic_2d.f90 "2D", @ref poisson_basic_3d.f90 "3D": Example showing how to use multigrid to solve Poisson's equation, and compare with an analytic solution
* @ref poisson_benchmark_2d.f90 "2D", @ref poisson_benchmark_3d.f90 "3D":
  Programs to benchmark the multigrid routines
* @ref poisson_cyl.f90 "2D": Solving Poisson's equation in cylindrical coordinates
* @ref poisson_cyl_dielectric.f90 "2D": Solving Poisson's equation in
  cylindrical coordinates with a jump in \f$\varepsilon\f$

# Implicit diffusion

* @ref implicit_diffusion_2d.f90 "2D", @ref implicit_diffusion_3d.f90 "3D": Using
  the multigrid methods to solve the diffusion equation with a backward Euler
  scheme

# Other

* @ref simple_streamer_2d.f90 "2D": A simplified model for streamers in 2D

