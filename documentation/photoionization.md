# Photoionization

[TOC]

# Introduction

TODO

# Source term

TODO

# Helmholtz approach

TODO

# Monte Carlo approach

In our Monte Carlo method for photoionization, we assume that photon scattering
can be neglected, and that the direction of ionizing photons is isotropic. The
method can be divided in three parts:
	-  Determine the coordinates at which photons are produced, and
store these coordinates in a list \f$L_{src}\f$.
	- Determine the coordinates at which photons are absorbed, and
store these coordinates in a list \f$L_{dst}\f$.
	- Compute the resulting photoionization profile on a mesh.

The implementation of these steps is described below.

## The source of ionizing photons

In a plasma fluid simulation, the computational domain is divided into cells. We
assume that for each cell the production of ionizing photons \f$I\f$ is known within
a given time step \f$\Delta t\f$. In our Monte Carlo method, this information is converted
to a list  \f$L_{src}\f$ of approximately \f$N\f$ discrete photons in the following way.
1. Determine the total photon production \f$I_{\Sigma}\f$ within the time step \f$\Delta t\f$,
by summing over all the cells.
2. For each cell \f$n_{\gamma} = N I / I_{\Sigma}\f$ photons should be produced.
To convert \f$n_{\gamma}\f$ to an integer, first draw a uniform random number
\f$U (0 , 1)\f$.
If \f$n_{\gamma} - \lfloor{ n_{\gamma}} \rfloor > U\f$ round up,
else round down.
3. For each produced photon, add the coordinate of the cell center to the list
\f$L_{src}\f$ of photons.

Some additional remarks: In principle the number of photons \f$N\f$ can be chosen
adaptively, for example to create discrete ‘super-photons’ with a given weight.
It is also possible to assign different weights to different photons by modifying
the above procedure, see section 11.3.4.
Instead of rounding \f$n_{\gamma}\f$ to an integer as described above, it can be more
realistic to sample from the Poisson distribution with parameter \f$n_{\gamma}\f$.
For \f$n_{\gamma} \ll 1\f$ the result would be almost identical, but for larger
values there are differences.
However, most of the random fluctuations in our method are due to the stochastic
photon direction and absorption, as discussed in the next section. Therefore our
rounding with uniform random numbers should be fine for most applications.
When grid cells are large compared to typical photoionization length scales,
one could determine a ‘subgrid’ production location in the cell, instead of using
the cell center.

## Absorption of ionizing photons

Now that we know where photons are produced, we need to determine where
they are absorbed. This is done in the following way. First, we determine
the absorption distance \f$r\f$ for a photon. Given an absorption function \f$f(r)\f$,
the cumulative absorption function \f$F(r) =\int_0^r f(r')dr'\f$ can be computed, either
numerically or analytically. Then a so-called <b>lookup table</b> is constructed with
columns \f$F(r)\f$ and \f$r\f$. Now we can do inversion sampling: given a random number
\f$U(0, 1)\f$, the corresponding distance \f$r\f$ is obtained by linear interpolation from the
table.
A random orientation \f$\hat{r}\f$ for the photon is determined, using a procedure
for picking a point on a sphere proposed by Marsaglia [205]:
1. Get two random numbers \f$U_1(-1,1)\f$ and \f$U_2(-1,1)\f$. If
\f$U_1^2 + U_2^2 \leq 1\f$ accept them, otherwise draw new numbers.
2. Let \f$a = U_1^2 + U_2^2\f$ and \f$b = 2 \sqrt{1-a}\f$.
The isotropic random direction is
\f$\hat{r} = (b U_1, b U_2, 1 - 2a)\f$ in Cartesian coordinates.
In the special case of a 2D Cartesian \f$(x, y)\f$ coordinate system, we should not
pick points on a sphere but points on a circle. Then the second step is replaced
by
\f[
\hat{r} = ( \frac{U_1^2-U_2^2}{U_1^2+U_2^2},\frac{2 U_1 U_2}{U_1^2+U_2^2}) .
\f]
Given the direction and the distance, the location of absorption
\f$\mathbf{r} = r \hat{\mathbf{r}}\f$ is known,
which is added to the list \f$L_{dst}\f$.
Sometimes, the typical absorption length of photons is much larger than the
domain size. Since most photons will not be absorbed within the domain, the
Monte Carlo approximation becomes less accurate. This can be resolved by
changing the lookup table. Suppose that all photons with absorption distances
\f$r > r_{max}\f$ can be ignored, then one can construct a lookup table up to
\f$r_{max}\f$ and scale the corresponding \f$F(r)\f$ values to the range (0; 1).
Each photon that is produced now gets a smaller weight \f$F(r_{max})\f$.

A Lookup table with the absorption locations of the photons has been created.
These coordinates need to be mapped to a mesh to get the photoionization profile
Two options can be considered: the photons can be absorbed on a single mesh with a constant grid spacing, or they can be absorbed at different levels in a refined mesh. 
In both cases the nearest grid point (NGP) is used to map the photons to densities,
which means that photons are absorbed within a single cell.
With linear (and higher order) interpolation the resulting density profiles are smoother,
but it becomes harder to handle refinement boundaries.
