# Transport data

[TOC]

# Introduction {#td-intro}

This page describes how to obtain electron transport data for use in plasma fluid simulations. Typical transport parameters are:

* The [electron mobility](https://en.wikipedia.org/wiki/Electron_mobility) \f$ \mu_e \f$
* The electron diffusion coefficient \f$ D_e \f$. Simulations are usually not very sensitive to this coefficient, so for simplicity a scalar diffusion coefficient is typically used.
* Ionization coefficient(s) \f$ \alpha \f$, which describe the mean number of ionization events per electron per unit distance
* Attachment coefficient(s) \f$ \eta \f$, which describe the mean number of attachment events per electron per unit distance
* Electron-neutral [rate coefficients](https://en.wikipedia.org/wiki/Reaction_rate_constant)

# Cross sections {#td-crosssec}

Transport coefficients are typically computed from electron-neutral [cross
sections](https://en.wikipedia.org/wiki/Cross_section_(physics)) using a
so-called Boltzmann solver, see @ref td-boltzmann. The first step in computing
transport coefficients is therefore to select a set of cross sections for a
given gas mixture. The most convenient source for such cross sections is
[lxcat](lxcat.net), which aims to
> address, at least in part, the well-recognized needs for the community to
> organize the means of collecting, evaluating and sharing data both for
> modeling and for interpretation of experiments.
Electron-neutral cross sections are usually given in tabulated form for a range of electron energies.

To download electron-neutral cross sections from lxcat, take the following steps:

1. Go to **data center**, **browse and download**
2. Highlight **scattering cross sections** and **electrons**
3. Next highlight all databases (unless you already know what you are looking for)
4. Highlight the ground states present in your gas mixture, for example **N2**,
   **O2** and **Ar**
5. Next, the types of cross sections has to be selected. Highlight all physical
   processes: **attachment**, **ionization**, **rotational** and **excitation**.
   Furthermore, a choice has to be made for the [elastic
   collision](https://en.wikipedia.org/wiki/Elastic_collision) cross sections,
   choosing between **elastic** and **effective**. This choice is discussed
   below in @ref td-elastic-cs, but you can simply select both here.
6. A list of cross section databases is now shown. It is important to use a
   so-called *complete set* of cross sections, which means that all relevant
   electron energy loss processes are accounted for. For each gas component,
   select such a set using the `+` button.
7. Next, you can see a picture of the included cross sections and download the
   data in **TXT** format. The beginning of the text file describes the format
   that is used. Furthermore, each database has a description with often
   important information, and instructions on how to cite the database, which is
   also important.

Next, the downloaded cross sections can be used in a Boltzmann solver, see @ref td-boltzmann.

## Elastic collision cross sections {#td-elastic-cs}

There exist various types of elastic collision cross sections on
[lxcat](lxcat.net), categorized under *elastic* and *effective*. Before, there
were also *momentum* cross sections, but a [post from Sergey
Pancheshnyi](https://us.lxcat.net/discussion_board.php?place=topic%2Flxcat%2F3qPvcBhvmhg%2Fdiscussion)
indicates that these have been renamed to *effective*. According to lxcat
documentation:

> "elastic" is used to denote the elastic momentum transfer cross section and
> where "effective" denotes the total momentum transfer cross section (sum of
> elastic momentum transfer and total inelastic cross sections). The latter is
> useful for solving the Boltzmann equation in the 2-term approximation.

Some further documentation is present in a [post from Leanne Pitchford](https://us.lxcat.net/discussion_board.php?place=topic%2Flxcat%2FrF5drzG_WjQ%2Fdiscussion):

> The cross sections for each species include EITHER elastic cross
> sections OR momentum transfer cross sections.
> - Elastic cross sections : these correspond to elastic momentum transfer cross sections.
> - Momentum transfer cross sections : these correspond to a total
>   momentum transfer or "effective momentum transfer" cross sections.
>   These include the effects of inelastic collisions as is appropriate
>   for use in the two-term spherical harmonic  expansion.  See, for
>   example, Baraff and Buchsbaum, Phys. Rev. 130, 1007 (1963) and Sec.
>   IIB of Pitchford and Phelps, Phys. Rev. A 25, 540 (1982).  Where data
>   is available, the effective Qm is set equal to the sum of the
>   inelastic cross sections plus the elastic
>   momentum transfer cross section. This is an approximate relation.

In other words, the *effective* cross sections include the momentum loss of all
types of electron-neutral collisions, not only elastic processes. It depends on
the type of Boltzmann solver (see @ref td-boltzmann) which type of cross section
can be used. *Effective* cross sections should probably only be used with a
two-term Boltzmann solver.

A further complication is that different types of *elastic* cross sections can
be found, for example the *elastic momentum transfer cross section* and the
*elastic total cross section*. The main idea is the following: when an electron
elastically collides with a gas molecule, the scattering angle for the electron
is generally anisotropic. In principle, such collisions should be described with
[differential cross
sections](https://en.wikipedia.org/wiki/Cross_section_(physics)#Differential_cross_section),
that give the probability of scattering per unit angle \f$ \frac{d\sigma}{d\Omega} \f$.
These differential cross sections can be integrated to give the total cross section
\f[
\sigma_{\mathrm{total}} = \int_{0}^{2\pi} \int_{0}^{\pi} \frac{d\sigma}{d\Omega} \, \sin(\theta) d\theta d\phi
\f]
which describes the probability of an elastic collision, regardless of scattering angle. Note that there is usually no dependence on the \f$ \phi \f$ direction, so that this just gives a factor \f$ 2 \pi \f$. In practice measuring or computing differential cross sections is difficult. The main effect these collisions is to reduce the momentum of electrons. How much momentum is lost is described by the *momentum transfer cross section*, given by
\f[
\sigma_{\mathrm{mt}} = \int_{0}^{2\pi} \int_{0}^{\pi} [1- \cos(\theta)] \frac{d\sigma}{d\Omega} \, \sin(\theta) d\theta d\phi
\f]
The factor \f$ 1- \cos(\theta) \f$ here weighs the cross section by the fraction of momentum that is lost. In a model it is convenient to assume that electrons scatter isotropically, meaning that all angles are equally probably. The differential cross section then simplifies to \f[ \frac{d\sigma}{d\Omega} = \frac{\sigma_{\mathrm{isotropic}}}{4 \pi} \f] With this assumption, the momentum transfer cross section is simply given by \f[ \sigma_{\mathrm{mt}} = \sigma_{\mathrm{isotropic}} \f] In other words, **one can directly use the momentum transfer cross section in an isotropic scattering model, and have the correct amount of momentum loss**.

## Three-body processes {#td-3body}

One has to be careful with three-body processes, such as three-body attachment
to oxygen: `e + O2 + O2 -> O2- + O2`. There are two options. The first is to
multiply the corresponding cross section with the density of the third body (in
this case O2) in the *correct units* before applying a Boltzmann solver. The
second is to essentially have a negligible cross section for this process in the
Boltzmann solver, but to scale it back to the desired amplitude afterwards in
the simulation model. This second option is convenient not fully consistent. It
can lead to problems when the three-body process has a significant effect on the
shape of the electron distribution function, for example at high pressure.

# Boltzmann solvers {#td-boltzmann}

Transport data can be computed from cross sections with a so-called Boltzmann solver. The idea is to solve [Boltzmann's equation](https://en.wikipedia.org/wiki/Boltzmann_equation) for electrons under simplified conditions, for example by assuming that the electron distribution function \f$ f(\vec{x}, \vec{v}, t) \f$ is homogeneous in space, and that the electric field is constant and homogeneous. A few Boltzmann solvers are listed below.

* [Bolsig+](http://www.bolsig.laplace.univ-tlse.fr) (Popular two-term Boltzmann
  solver, closed source)
* [Bolos](https://github.com/aluque/bolos) (Open-source two-term Boltzmann
  solver, implemented in Python)
* [Magboltz](http://consult.cern.ch/writeup/magboltz/) (Particle / Monte Carlo
  method implemented in Matlab, open source)
* [Particle_swarm](https://gitlab.com/MD-CWI-NL/particle_swarm) (Particle /
  Monte Carlo method, open source, implemented in Fortran)

## Bolsig+ {#td-bolsigplus}

An easy way to use Bolsig+ is through the online interface at [lxcat](lxcat.net)
via the **online calculations** button. (The online version has some downside however, see below!) After selecting cross sections from the desired database(s), you could for example use the following options

option | value | description
---|---|---
EEDF | non-Maxwellian | The EEDF is generally non-Maxwellian
E/N range | 1-1200 Td | A reasonable range for streamer simulations
`# of points` | 100 | This seems to be the limit for online calculations
E/N profile | exponential | To have more points at lower fields, quadratic could also work
`Tgas` | 300 K | Gas temperature
Super-elastic collisions | include | Probably not significant at low gas temperature and high E/N

After pressing **run calculations**, you can download the swarm parameters in a text file. This file can be converted to input for the afivo-streamer code with the script `tools/bolsig_convert.py`, for example like:

    ./bolsig_convert.py swarm.txt my_new_file.txt

where `swarm.txt` is the output of Bolsig+. Execute `bolsig_convert.py` with the `-h` options to see its full usage. The output should now contain the following:

* `Mean energy (eV)`
* `Mobility *N (1/m/V/s)`
* `Diffusion coefficient *N (1/m/s)`
* `Townsend ioniz. coef. alpha/N (m2)`
* `Townsend attach. coef. eta/N (m2)`
* Data for specific collisions, either selected by type or by name. By default ionization and attachment processes are included.

One still has to define the reactions before the file can be used as input for the afivo-streamer code, see @ref td-format

A more advanced way to compute transport data is to download an offline version of a Boltzmann solver and compute the data locally, for example using [BOLSIG+](http://www.bolsig.laplace.univ-tlse.fr/download.html) (Windows), [BOLSIG-](http://www.bolsig.laplace.univ-tlse.fr/download.html) (linux), or [particle_swarm](https://gitlab.com/MD-CWI-NL/particle_swarm). When input data is linearly interpolated, it is best to compute a large number of data points to reduce interpolation errors, in particular in the rate coefficients, which are often far from linear. Another option is to enable cubic interpolation of the input data, see @ref td-interpolation.

# Transport data format {#td-format}

Finally, a list of reactions has to be added by hand, since it is hard to deduce
the species appearing on the right-hand side from the output of a Boltzmann
solver. The format of these reactions is illustrated in the files in the
`transport_data` folder. A standard file for `O2` could look like this

    reaction_list
    -----------------------
    # Ionization
    e + O2 -> e + e + O2+,field_table,C43 O2 Ionization 12.06 eV
    # Attachment
    e + O2 + O2 -> O2-,field_table,C27 O2 Attachment
    e + O2 -> O-,field_table,C28 O2 Attachment
    -----------------------

In this example, the reactions are defined by a `field_table`, which tabulates
the reaction rates versus reduced electric field. The strings *C43 O2 Ionization
12.06 eV* etc. refer to the names of the corresponding tables. Note that the
three-body reaction rate constant here will be multiplied by the density of O2
squared. The cross section for this process should therefore not have been
pre-multiplied by the O2 density when it was used in the Boltzmann solver.

Further details about the format of chemical reactions can be found in the
[chemistry documentation](documentation/chemistry.md).

# Transport data interpolation {#td-interpolation}

To speed up computations, transport and reaction coefficients are stored in lookup tables, see [lookup_table_fortran](https://github.com/jannisteunissen/lookup_table_fortran). These tables have a fixed size and spacing, so that values can quickly be obtained without searching through a list. These tables are constructed and used as follows:

**Step 1** The input data is read in from a text file

**Step 2** The input data is interpolated to have a regular spacing in `E/N`. This interpolation is controlled by the following parameters:

    table_data%min_townsend = 0. (minimum E/N)
    table_data%max_townsend = 1000. (maximum E/N)
    table_data%size = 1000 (number of points in the table)
    table_data%input_interpolation = linear (or cubic_spline)

Note there are two options for the interpolation of input data: linear (the default) and cubic spline interpolation. The interpolated data is then stored in a lookup table.

**Step 3** When a value is required at a certain `E/N`, the values in the lookup table are *linearly interpolated*.
