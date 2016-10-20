# Geometric multigrid

[TOC]

\todo Add references

# Introduction {#mg-intro}

Multigrid can be seen as a technique to improve the convergence of a relaxation
method, by using a hierarchy of grids.
Afivo comes with a built-in geometric multigrid solver, to solve problems of the form
\f[
  A(u) = \rho,
\f]
where \f$A\f$ is a (nearly) elliptic operator, \f$\rho\f$ the right-hand side and
\f$u\f$ the solution to be computed.
In discretized form, we write the above equation as 
\f[
  A_h(u_h) = \rho_h,
\f]
where \f$h\f$ denotes the mesh spacing at which the equation is discretized.

There already exists numerous sources on the foundations of multigrid, the
different cycles and relaxation methods, convergence behavior and other aspects.
Here we will not provide a general introduction to multigrid. Instead we briefly
summarize the main ingredients, and focus on one particular topic: how to
implement multigrid on an adaptively refined quadtree or octree mesh. On such a
mesh, the solution has to be specified on all levels. Therefore we use FAS
multigrid, which stands for Full Approximation Scheme. Below, the implementation
of the various multigrid components in Afivo are described.

# The V-cycle {#mg-v-cycle}

Suppose there are levels \f$l = l_\mathrm{min}, l_\mathrm{min}+1, \ldots, l_\mathrm{max}\f$, then the FAS
V-cycle can be described as

 1. For \f$l\f$ from \f$l_\mathrm{max}\f$ to \f$l_\mathrm{min}+1\f$, perform
    \f$N_\mathrm{down}\f$ relaxation steps on level \f$l\f$, then update level
    \f$l-1\f$ (see below)
 2. Perform \f$N_\mathrm{base}\f$ relaxation steps on level
    \f$l_\mathrm{min}\f$, or apply a direct solver
 3. For \f$l\f$ from \f$l_\mathrm{min}+1\f$ to \f$l_\mathrm{max}\f$, perform a
	correction using the data from level \f$l-1\f$ <a name="eq_coarse-corr"
	>(3)</a> \f[ u_h \leftarrow u_h + I_H^h(v_H - v'_H), \f] then perform
	\f$N_\mathrm{up}\f$ relaxation steps on level \f$l\f$. (See below for the
	notation)

The first two steps require some extra explanation. Let us denote the level
\f$l-1\f$ grid by \f$H\f$ and the level \f$l\f$ grid by \f$h\f$, and let \f$v\f$
denote the current approximation to the solution \f$u\f$. Furthermore, let
\f$I_H^h\f$ be an interpolation operator to go from coarse to fine and
\f$I_h^H\f$ a restriction operator to go from fine to coarse. For these
operators, the schemes described in section \ref main-interp-restrict are used.

In the first step, the coarse grid is updated in the following way

 1. Set \f$v_H \leftarrow I_h^H v_h\f$, and store a copy \f$v'_H\f$ of \f$v_H\f$
 2. Compute the fine grid residual \f$r_h = \rho_h - A_h(v_h)\f$
 3. Update the coarse grid right-hand side
	\f[
	  \rho_H \leftarrow I_h^H r_h + A_H(v_H).
	\f]

This last equation can also be written as
\f[
  \rho_H \leftarrow I_h^H \rho_h + \tau_h^H,
\f]
where \f$\tau_h^H\f$ is given by
\f[
  \tau_h^H = A_H(I_h^H v_h) - I_h^H A_h(v_h).
\f]
This term can be seen as a correction to \f$\rho\f$ on the coarse grid.
When a solution \f$u_h\f$ is found such that \f$A_h(u_h) = \rho_h\f$, then
\f$u_H = I_h^H u^h\f$ will be a solution to \f$A_H(u_H) = \rho_H\f$.

In the second step, relaxation takes place on the coarsest grid. In order to
quickly converge to the solution with a relaxation method, this grid should
contain very few points (e.g., \f$2\times 2\f$ or \f$4\times 4\f$ in 2D). Alternatively,
a direct solver can be used on the coarsest grid, in which case it can contain more points.
Such a direct method has not yet been built into Afivo, although this is
planned for the future. As a temporary solution, additional coarse grids can be
constructed below the coarsest quadtree/octree level. For example, if a quadtree
has boxes of \f$16\times 16\f$ cells, then three levels can be added below it
(\f$8\times 8\f$, \f$4\times 4\f$ and \f$2\times 2\f$), which can then be used in the
multigrid routines.

# The FMG-cycle {#mg-fmg-cycle}

The full multigrid (FMG) cycle that is implemented in Afivo works in the
following way:

1. If there is no approximation of the solution yet, then set the initial guess
   to zero on all levels, and restrict \f$\rho\f$ down to the coarsest grid
   using \f$I_h^H\f$. If there is already an approximation \f$v\f$ to the
   solution, then restrict \f$v\f$ down to the coarsest level. Use <a
   href="#eq_coarse-rhs"> equation(4)</a> to set \f$\rho\f$ on coarse grids.
2. For \f$l = l_\mathrm{min}, l_\mathrm{min}+1, \ldots, l_\mathrm{max}\f$
    * Store the current approximation \f$v_h\f$ as \f$v`_h\f$
    * If \f$l > l_\mathrm{min}\f$, perform a coarse grid correction using
      <a href="#eq_coarse-corr" >equation(3)</a>
    * Perform a V-cycle starting at level \f$l\f$, as described in the previous
      section

# Gauss-Seidel red-black {#main-gsrb}

In Afivo, we have implemented Gauss-Seidel red-black or GS-RB as a relaxation
method.
This method is probably described in almost all textbooks on multigrid, such as
Xcite{Brandt_2011,trottenberg2000multigrid,Briggs_2000,Hackbusch_1985}, so we
just give a very brief description.

The *red-black* refers to the fact that points are relaxed in an
alternating manner, using a checkerboard-like pattern. For example, in two
dimensions with indices \f$(i,j)\f$ points can be labeled *red* when \f$i+j\f$ is
even and *black* when \f$i+j\f$ is odd. Now consider 
<a href="#eq_mg-discr-equation">equation(2)</a>, which typically relates a value \f$u_h^{(i,j)}\f$ to
neighboring values and the source term \f$\rho\f$. If we keep the values of the
neighbors fixed, then we can determine the value \f$u_h^{(i,j)}\f$ that locally
solves the linear equation. This is precisely what is done in GS-RB: the linear
equations are solved for all the red points while keeping the old black values,
and then vice versa.


# Conservative filling of ghost cells {#mg-ghost-cells}

The finer levels will typically not cover the complete grid in Afivo, so that
ghost cells have to be used near refinement boundaries.
These ghost cells can be filled in multiple ways, which will affect the
multigrid solution and convergence behavior.
Here we consider *conservative* schemes for filling ghost cells
Xcite{Bai_1987,trottenberg2000multigrid}.
A conservative scheme ensures that the coarse flux across a refinement boundary
equals the average of the fine fluxes, see <a href="#fig_mg-ref-bound" >Figure 8</a>.

Ensuring consistent fluxes near refinement boundaries helps in obtaining a
*consistent* solution. For example, if we consider a general equation of
the form \f$\nabla \cdot \vec{F} = \rho\f$, then the divergence theorem gives
\f[
  \int_V \rho \, dV = \int_V \nabla \cdot \vec{F} \, dV = \int \vec{F} \cdot
  \vec{n} \, dS,
\f]
where the last integral runs over the surface of the volume \f$V\f$, and \f$\vec{n}\f$
is the normal vector to this surface.
This means that when fine and coarse fluxes are consistent, the integral over
\f$\rho\f$ will be same on the fine and the coarse grid.

![Two coarse cells, of which the right one is refined. The cell centers are indicated by dots. There are two ghost values (red dots) on the left of the refinement boundary. Fluxes across the refinement boundary are indicated by arrows.](mg_refinement_boundary.png)

The construction of a conservative scheme for filling ghost cells is perhaps
best explained with an example.
Consider a 2D Poisson problem
\f[
  \nabla^2 u = \nabla \cdot (\nabla u) = \rho,
\f]
with a standard 5-point stencil for the Laplace operator
<a name="eq_5-point-stencil" >(6)</a>
\f[
  L_h = h^{-2}
  \begin{bmatrix}
    & \;1 &\\
    1 & -4 & 1\\
    & \;1 &
  \end{bmatrix}.
\f]
With this stencil, the coarse flux \f$f_H\f$ across the refinement boundary in
 <a href="#fig_mg-ref-bound" >Figure 8</a> is given by
\f[
  f_H = [u_H^{(2,1)} - u_H^{(1,1)}]/H,
\f]
and on the fine grid, the two fluxes are given by
\begin{align}
  f_{h,1} &= [u_h^{(3,1)} - g_h^{(2,1)}]/h,\\
  f_{h,2} &= [u_h^{(3,2)} - g_h^{(2,2)}]/h.
\end{align}
The task is now to fill the ghost cells \f$g_h^{(2,1)}\f$ and \f$g_h^{(2,2)}\f$ in such
a way that the coarse flux equals the average of the fine fluxes, i.e., such
that
<a name="eq_mg-flux-cons" >(7)</a>
\f[
  f_H = (f_{h,1} + f_{h,2})/2
\f]
To relate \f$u_H^{(2,1)}\f$ to the refined values \f$u_h\f$, the restriction operator
\f$I_h^H\f$ needs to be specified.
In our implementation, this operator does averaging over the children, which can
be represented as
\f[
  I_h^H = \frac{1}{4}
  \begin{bmatrix}
    1 & 1\\
    1 & 1
  \end{bmatrix}.
\f]
The constraint from <a href="#eq_mg-flux-cons">equation(7)</a> can then be written as
\f[
  g_h^{(2,1)} + g_h^{(2,2)} = u_H^{(1,1)} + \frac{3}{4} \left(u_h^{(3,1)} + u_h^{(3,2)}\right)
  - \frac{1}{4} \left(u_h^{(4,1)} + u_h^{(4,2)}\right).
\f]
Any scheme for the ghost cells that satisfies this constraint will be a
conservative discretization.

Bilinear *extrapolation* (similar to standard bilinear interpolation) gives
the following scheme for \f$g_h^{(2,1)}\f$
\f[
  g_h^{(2,1)} = \frac{1}{2} u_H^{(1,1)} + \frac{9}{8} u_h^{(3,1)} -
  \frac{3}{8} \left (u_h^{(3,2)} + u_h^{(4,1)} \right)
  + \frac{1}{8} u_h^{(4,2)}.
\f]
(The scheme for \f$g_h^{(2,2)}\f$ should then be obvious.)
Another option is to use only the closest two neighbors for the extrapolation,
which gives the following expression for \f$g_h^{(2,1)}\f$
<a name="eq_ghost-cell-standard-2d" >(8)</a>
\f[
  g_h^{(2,1)} = \frac{1}{2} u_H^{(1,1)} + u_h^{(3,1)} -
  \frac{1}{4} \left (u_h^{(3,2)} + u_h^{(4,1)} \right).
\f]
This last scheme is how refinement-boundary ghost cells are filled by default in
Afivo.


## Three-dimensional case {#main-3d_case}

In three spatial dimensions, the 5-point stencil of
<a href="#eq_5-point-stencil">equation(6)</a> becomes a 7-point stencil with \f$-6\f$ at the center, and
the restriction operator has eight entries of 1/8.
The analog of <a href="#eq_ghost-cell-standard-2d">equation(8)</a> then becomes
<a name="eq_ghost-cell-standard-3d" >(9)</a>
\f[
  g_h^{(2,1,1)} = \frac{1}{2} u_H^{(1,1,1)} + \frac{5}{4} u_h^{(3,1,1)} -
  \frac{1}{4} \left (u_h^{(4,1,1)} + u_h^{(3,2,1)} + u_h^{(3,1,2)}\right).
\f]


## Change in coefficient at cell face {#mg-varepsilon}

For the more general equation with a coefficient \f$\varepsilon\f$
\f[
  \nabla \cdot (\varepsilon \nabla u) = \rho,
\f]
we consider a special case: \f$\varepsilon\f$ jumps from \f$\varepsilon_1\f$ to
\f$\varepsilon_2\f$ at a cell face.
Local reconstruction of the solution shows that a gradient \f$(\phi_{i+1} -
\phi_i) / h\f$ has to be replaced by
\f[
  \frac{2 \, \varepsilon_1 \varepsilon_2}{\varepsilon_1 \varepsilon_2} \, \frac{\phi_{i+1} -
    \phi_i} {h},
\f]
or in other words, the gradient is multiplied by the harmonic mean of the
\f$\varepsilon\f$`s (see for example chapter 7.7 of
Xcite{trottenberg2000multigrid}).
The 5-point stencil for the Laplacian can be modified accordingly.

When a jump in \f$\varepsilon\f$ occurs on a coarse cell face, it will also be
located on a fine cell face, see <a href="#fig_mg-ref-bound" >Figure 8</a>.
In this case, the ghost cell schemes described above for constant \f$\varepsilon\f$
still ensure flux conservation.
The reason is that the coarse and fine flux are both weighted by a factor
\f$2 \, \varepsilon_1 \varepsilon_2 / (\varepsilon_1 \varepsilon_2)\f$.


## Cylindrical case {#mg-cyl}

In cylindrical coordinates, the Laplace operator can be written as
\f[
  \nabla^2 u = \frac{1}{r} \partial_r (r \partial_r u) + \partial^2_z u
  = \partial_r^2 u + \frac{1}{r} \partial_r u + \partial_z^2 u,
\f]
where we have assumed cylindrical symmetry (no \f$\phi\f$ dependence).
At a radius \f$r \neq 0\f$, the 5-point stencil is
<a name="eq_5-point-stencil-cyl" >(10)</a>
\f[
  L_h = h^{-2}
  \begin{bmatrix}
    & \;1 &\\
    1-\frac{h}{2 r} & -4 & 1+\frac{h}{2 r}\\
    & \;1 &
  \end{bmatrix}.
\f]
With the cell-centered grids in Afivo, radial grid points are located at
\f$(i - \frac{1}{2}) h\f$ for \f$i = 1, 2, 3, \ldots\f$, which means we do not have to
consider the special case \f$r = 0\f$.
For this type of grid indexing, the 5-point stencil can also be written as
\f[
  L_h = h^{-2}
  \begin{bmatrix}
    & \;1 &\\
    \frac{2i-2}{2i-1} & -4 & \frac{2i}{2i-1}\\
    & \;1 &
  \end{bmatrix}.
\f]

If we do not modify the restriction operator, then the ghost cells can still be
filled with the schemes from <a href="#eq_ghost-cell-standard-2d">equations(8)</a> and
<a href="#eq_ghost-cell-standard-3d">(9)</a>.
One way to interpret this is that fluxes are computed in the same way in
cylindrical coordinates, although their divergence is weighted by the radius:
\f[
  \nabla \cdot \vec{F} = \frac{1}{r} \partial_r (r F_r) + \ldots
\f]
From <a href="#fig_mg-ref-bound">Figure 8</a>, we can see that for refinement in the
\f$r\f$-direction, the coarse and fine flux are `weighted` by the same radius.
For the fluxes in the \f$z\f$-direction, the computations are the same as for the
Cartesian case.

\verbatim
Note that when the restriction operator is changed to include radial weighting,
these arguments are no longer valid.
\endverbatim

# Multigrid test problems {#mg-examples}

In this section we present several Poisson test problems
to demonstrate the multigrid behavior on a partially refined mesh.
We use the **method of manufactured solutions**: from an analytic solution the right-hand side and boundary
conditions are computed. Two test problems are considered, a
constant-coefficient
<a class="code" href="poisson__basic__2d_8f90.html#aef02b53cac21b72a47afc1c34286c443">two-</a> and
<a class="code" href="poisson__basic__3d_8f90.html#a1b0ee6adabd66bb15df1120e11776ead">three-dimensional</a>
Poisson equation shown in <a href="examples.html"><span>Examples</span></a>
<a name="eq_mg-example-lpl-1" >(11)</a>
\f[
  \nabla^2 u = \nabla \cdot (\nabla u) = \rho
\f]
and also a 
<a class="code" href="poisson__cyl_8f90.html#aaa62a5d2b70da45c7561014a92ccb9f0">cylindrical</a>
version with a coefficient \f$\varepsilon\f$
<a name="eq_mg-example-lpl-2" >(12)</a>
\f[
  \frac{1}{r} \partial_r (r \varepsilon \partial_r u) + \partial_z (\varepsilon \partial_z u) = \rho,
\f]
both on a two-dimensional rectangular domain \f$[0,1] \times [0,1]\f$.
For the cylindrical case, \f$\varepsilon\f$ has a value of \f$100\f$ in the lower left
quadrant \f$[0,0.25] \times [0,0.25]\f$, and a value \f$1\f$ in the rest of the domain.
In all cases, we pick the following solution for \f$u\f$
\f[
  u(r) = \exp(|{\vec{r}-\vec{r}_1}|/\sigma) + \exp(|{\vec{r}-\vec{r}_2}|/\sigma),
\f]
where \f$\vec{r_1} = (0.25, 0.25)\f$, \f$\vec{r_2} = (0.75, 0.75)\f$ and
\f$\sigma = 0.04\f$.
An analytic expression for the right-hand side \f$\rho\f$ is obtained by plugging
the solution in <a href="#eq_mg-example-lpl-1">equation(11)</a> and
<a href="#eq_mg-example-lpl-2">equation(12)</a>.
Note that jumps in \f$\varepsilon\f$ also contribute to the source term \f$\rho\f$,
and the solution itself is used to set boundary conditions.

The two different problems can now be solved numerically.
For the cylindrical case with the varying \f$\varepsilon\f$, a modified Laplacian
operator is used, as described in section \ref mg-ghost-cells.
The Gauss-Seidel red-black relaxation methods are also modified, because they
depend on the applied operator, see section \ref main-gsrb.
For these examples, we have used
\f$N_\mathrm{down} = N_\mathrm{up} = N_\mathrm{base} = 2\f$ (number of down/up/base
smoothing steps), and a coarsest grid of \f$2\times 2\f$ cells.

It is possible to do adaptive mesh refinement in multigrid, for example by using
an estimate of the local truncation error based on <a href="#eq_mg-tau">equation(5)</a>
(see also chapter 9 of Xcite{trottenberg2000multigrid}).
Such a technique is not used here, instead the refinement criterion is based on
the right-hand side: refine if \f$\Delta x^2 |\rho| > 5.0 \times 10^{-4}\f$.
The resulting mesh spacing is shown in <a href="#fig_mg-ex1">Figure 9a</a>.

In both cases, one FMG (full multigrid) cycle is enough to achieve convergence
up to the discretization error, which was approximately \f$10^{-4}\f$ for the mesh
of <a href="#fig_mg-ex1">Figure 9a</a>.
Consecutive FMG cycles have a negligible effect on the absolute error, although
they do reduce the residual \f$r = \rho - \nabla^2 u\f$.
The maximum value of \f$|r|\f$ is shown versus iteration number in
<a href="#fig_mg-ex1">Figure 9b</a>.
The convergence behavior is similar for both cases, with each iteration
reducing the residual by a factor of about \f$0.056\f$.
The offset between the lines is caused by the \f$\varepsilon = 100\f$ region, which
locally amplifies the source term by a factor of 100.

<a name="fig_mg-ex1" />
a) Mesh spacing | b) Maximum residual versus FMG
------------- | -------------
![a](mesh_mg_ex_1.png) | b todo

a) The left figure shows the mesh spacing used for the multigrid examples, in a
\f$[0,1] \times [0,1]\f$ domain. Each step in color is a factor two in
refinement, with red indicating \f$\Delta x = 2^{-5}\f$ and the darkest blue
indicating \f$\Delta x = 2^{-12}\f$. b) The maximum residual versus FMG
iteration
