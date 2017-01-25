# Implementing a plasma fluid model {#main-example-fluid}


To illustrate how Afivo can be used, we describe the implementation of a simple
2D/3D *plasma fluid model* for streamer discharges below. For simplicity,
photoionization is not included in this example. A review of fluid models for
streamer discharges can be found in Xcite{Luque_2012}.


## Model formulation {#main-model-formulation}

We use the so-called drift-diffusion-reaction approximation:
<a name="eq_fluid-model" >(13)</a>
\f[
\begin{align}
  \partial_t n_e &= -\nabla \cdot \vec{j}_e + \bar{\alpha} |{\vec{j}_e}|,\\
  \partial_t n_i &= \bar{\alpha} |{\vec{j}_e}|,\\
  \vec{j}_e &= -\mu_e n_e \vec{E} - D_e \nabla n_e,
\end{align}
\f]

where \f$n_e\f$ is the electron density, \f$n_i\f$ the positive ion density, \f$\vec{j}_e\f$
the electron flux, \f$\bar{\alpha}\f$ the effective ionization coefficient, \f$\mu_e\f$
the electron mobility, \f$D_e\f$ the electron diffusion coefficient and \f$\vec{E}\f$
the electric field. The above equations are coupled to the electric field, which
we compute in the electrostatic approximation:
\f[
\begin{align}
  \vec{E} &= -\nabla \phi,\\
  \nabla^2 \phi &= -\rho/\varepsilon_0\\
  \rho &= e (n_i - n_e),
\end{align}
\f]
where \f$\phi\f$ is the electric potential, \f$\varepsilon_0\f$ the permittivity of
vacuum and \f$e\f$ the elementary charge. The electric potential is computed with
the multigrid routines described in section \ref main-afivo-multigrid.

We make use of the *local field approximation* Xcite{Li_2007}, so that
\f$\mu_e\f$, \f$D_e\f$ and \f$\bar{\alpha}\f$ are all functions of the local electric field
strength \f$E = |{\vec{E}}|\f$. These coefficients can be obtained
experimentally, or they can be computed with a Boltzmann solver
Xcite{Bolsighagelaar011,Dujko_2011} or particle swarms Xcite{Li_hybrid_i_2010}.


## Flux calculation and time stepping {#main-flux-calc-time-stepping}

The electron flux is computed as in Montijn \cite Montijn_2006. For the diffusive part,
we use central differences. The advective part is computed using the Koren
limiter Xcite{koren_limiter}. The Koren limiter was not designed to include
refinement boundaries, and we use linear interpolation to obtain fine-grid ghost
values. These ghost cells lie inside a coarse-grid neighbor cell, and we limit
them to twice the coarse values to preserve positivity. (We would like to
improve this in the future.)

Time stepping is also performed as in Montijn \cite Montijn_2006, using the explicit
trapezoidal rule, also known as the modified Euler`s method. The time step is
determined by a CFL condition for the electron flux and the dielectric
relaxation time, as in Montijn \cite Montijn_2006.

## Refinement criterion {#main-refinement-crit}

Our refinement criterion contains two components: a *curvature monitor*
\f$c_\phi\f$ for the electric potential and a monitor \f$\bar{\alpha} \Delta x\f$ which
gives information on how well the ionization length (\f$1/\bar{\alpha}\f$) is
resolved. For both, we use the maximum value found in a box in order to decide
whether to (de)refine it.

Since \f$\nabla^2 \phi = -\rho / \varepsilon_0\f$, the curvature monitor can be
computed as \f$c_\phi = \Delta x^2 |\rho| / \varepsilon_0\f$. The quantity
\f$\bar{\alpha} \Delta x\f$ is computed by locating the highest electric field in
the box, and looking up the corresponding value of \f$\bar{\alpha}\f$. The combined
refinement criterion is then as follows, where later rules can override earlier
ones:
	- If \f$\bar{\alpha} \Delta x < 0.1\f$ and \f$\Delta x < 25 \, \mu\textrm{m}\f$,
	derefine.
	- If \f$t < 2.5 \, \textrm{ns}\f$, ensure that there is enough refinement
	around the initial seed to resolve it.
	- If \f$\bar{\alpha} \Delta x > 1.0\f$ and \f$c_\phi > 0.1 \, \textrm{Volt}\f$,
	refine.

## Simulation conditions and results {#main-results}

<a name="fig_ex-streamer-seed" />
<img src="../../documentation/figures/streamer_seed.png" width=300 />
**Figure 10**. Cross section through the center of the three-dimensional simulation
domain. The ionized seeds with a density of \f$10^{20} \, \textrm{m}^{-3}\f$
electrons and ions are indicated in red. There is a background density of
\f$5 \times 10^{15} \, \textrm{m}^{-3}\f$ electrons and ions, and the background
electric field points down with a magnitude \f$E_0 = 2.5 \, \textrm{MV/m}\f$.

A cross section through the computational domain of \f$(32 \, \textrm{mm})^3\f$ is
shown in <a href="#fig_ex-streamer-seed" >Figure 10</a>. The background field points down,
with a magnitude \f$E_0 = 2.5 \, \textrm{MV/m}\f$, which is about
\f$5/6\f$\textsuperscript{th} of the critical field. The background field is applied
by grounding the bottom boundary of the domain, and applying a voltage at the
top. At the other sides of the domain we use Neumann boundary conditions for the
potential. We use transport coefficients (e.g., \f$\bar{\alpha}\f$, \f$\mu_e\f$) for
atmospheric air, but for simplicity photoionization has not been included.
Instead a background density of \f$5 \times 10^{15} \, \textrm{m}^{-3}\f$ electrons
and positive ions is present.

Two seeds of electrons and ions locally enhance the background electric field,
see <a href="#fig_ex-streamer-seed" >Figure 10</a>. These seeds have a density of
\f$10^{20} \, \textrm{m}^{-3}\f$, a width of about \f$0.3 \, \textrm{mm}\f$ and a length
of \f$1.6 \, \textrm{mm}\f$. The electrons from these seeds will drift upwards,
enhancing the field at the bottom of the seed where a positive streamer can
form.

In <a href="#fig_ex-streamer-dns">Figure 11</a>, the time evolution of the electron
density is shown, and in <a href="#fig_ex-streamer-fld">Figure 12</a>. the electric field
is shown. Two positive streamers grow downwards from the ionized seeds. The
upper one is attracted to the negatively charged end of the lower one, and
connects to it at around \f$9.5 \, \textrm{ns}\f$. The three-dimensional simulation
took about 3.5 hours on a 16-core machine, and eventually used about
\f$1.3 \times 10^7\f$ grid cells.
