# Time integration

# Time integration method

Time integration is performed explicitly, using the Afivo time integration
module, which is described in the [Afivo online
documentation](https://teunissen.net/afivo/md_documentation_time_integration.html).

Several time integrators can be used, which can be set using the `time_integrator` variable. The following options are currently supported:

* `forward_euler` (first order)
* `heuns_method` (second order, SSP)
* `midpoint_method` (second order)
* `ssprk3` (third order, SSP)

## Time step restrictions

Due to the explicit time integration, there are a number of time step restrictions:

* A [CFL condition](https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition)
* A condition due to explicit diffusion (usually not important, that's why diffusion is handled explicitly)
* A time step restriction due to the plasma chemistry
* A time step restriction due to the dielectric relaxation time

For more details about these time step restrictions, see e.g. the [paper about afivo-streamer](https://doi.org/10.1088/1361-6463/aa8faf). More on the dielectric relaxation time, and ways to avoid it, can be found [in this paper](https://doi.org/10.1088/1361-6595/ab6757).

The plasma chemistry time step restriction is of the form

    dt < max(n, n0) / (d/dt n)'

where `d/dt n` is the time derivative due to chemical reactions. When a density is
(almost) zero, the parameter `n0` avoids a very small time step. This parameter
can be set in the m_dt module.

## Relevant parameters

In the m_dt module, there are a couple of relevant parameters for the time integration:

\snippet m_dt.f90 relevant_parameters

## Abort on small time steps

When the time step becomes smaller than `dt_min`, the simulation will abort, assuming an instability has occurred. This can be prevented by:

* Decreasing the minimal time step `dt_min`, see m_dt
* Increasing the `n0` for the chemistry time step (see above)



