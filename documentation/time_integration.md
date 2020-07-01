# Time integration

The `m_af_advance` module can be used to perform time integration. 

## Methods

Currently, the following explicit methods are included:

* [Forward Euler](https://en.wikipedia.org/wiki/Euler_method) (`af_forward_euler`)
* [Midpoint method](https://en.wikipedia.org/wiki/Midpoint_method) (`af_midpoint_method`)
* [Heun's method](https://en.wikipedia.org/wiki/Heun%27s_method) (`af_heuns_method`)

New integrators can be added relatively easily by modifying the `af_advance` routine.

## How to use

Time integration can be performed with a call to `m_af_advance::af_advance`.

## User-supplied forward Euler method

To use the built-in time integration, a `forward_euler` routine has to be
provided, see `m_af_advance::subr_feuler` for details. This routine will then
be used to construct the various time integration schemes.

For higher-order schemes, the idea is that multiple copies will be stored for variables that will be integrated, which correspond to different temporal states. The number of copies can be specified in `m_af_core::af_add_cc_variable`.

The user-provided forward Euler method then takes three extra arguments that
indicate which temporal states to operate on. For a simple equation of the form

    y' = f(y)

the method should provide a solution

    y_out = y_prev + dt * f(y_deriv)

More details can be found at `m_af_advance::subr_feuler`, and an example from @ref euler_gas_dynamics.f90 is shown below

\snippet euler_gas_dynamics.f90 forward_euler_gasd

