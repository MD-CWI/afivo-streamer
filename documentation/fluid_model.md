# Fluid model equations

[TOC]

# Classical fluid model

By default, the "classical" fluid model is used, which is a drift-diffusion-reaction model with the local field approximation. The equations solved in this model are:

\f[
\partial_t n_{e} + \nabla \cdot \boldsymbol{\Gamma}_e = S_\mathrm{chem} + S_\mathrm{ph},
\f]
with the electron flux given by
\f[
\boldsymbol{\Gamma}_e = -\mu_e \boldsymbol{E} n_e - D_e \nabla n_e.
\f]
In these equations, the symbols have the following meaning:
* \f$\mu_e\f$ is the electron mobility
* \f$D_e\f$ the scalar electron diffusion coefficient
* \f$S_\mathrm{chem}\f$ is a source term due to chemical reactions (including ionization and attachment), see @ref documentation/chemistry.md
* \f$S_\mathrm{ph}\f$ is a source term due to photoionization, see @ref documentation/photoionization.md

With the local field approximation (LFA), electron transport coefficients are assumed to depend on the local electric field strength. In the code, so-called reduced transport coefficients are used. For example, the input data for the electron mobility is a table with rows

\f[\mathrm{E/N}, N \, \mu_e(\mathrm{E/N}),\f]

where \f$N\f$ is the gas number density and \f$\mathrm{E/N}\f$ is the reduced electric field in Townsend. The code then interpolates this table (see @ref td-interpolation) and multiplies with the inverse of \f$N\f$ to obtain the actual electron mobility.

## Ions

The equation for ions looks similar to that of electrons,
\f[
\partial_t n_{j} + \nabla \cdot \boldsymbol{\Gamma}_j = S_\mathrm{chem} + S_\mathrm{ph},
\f]
but there are a couple of differences:

* By default, the ion flux is assumed to be zero
* If a non-zero ion mobility has been specified using the parameters `input_data%mobile_ions` and `input_data%ion_mobilities`, then the ion flux is computed as \f$\boldsymbol{\Gamma}_j = \pm \mu_j \boldsymbol{E} n_e\f$, with the sign depending on the ion charge
* The photoionization source term only contributes to the production of one particular ion species

## Neutrals

Neutrals evolve only due to chemical reactions
\f[
\partial_t n_{j} = S_\mathrm{chem}.
\f]

# Fluid model with energy equation

An experimental feature is the use of a fluid model with an energy equation. In this case, transport and reaction data is interpolated based on the local mean electron energy. The code can automatically convert input data to this format.

The use of a different models is controlled by the following code:

\snippet m_model.f90 model_types


