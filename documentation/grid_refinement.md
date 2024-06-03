# Grid refinement {#grid-refinement}

## Parameter refine_adx
Grid refinement is controlled by several parameters. Usually the most import one is

    # Refine if alpha*dx is larger than this value:
    refine_adx = 1.0

The idea behind this parameter is that the electron density will locally approximately increase like \f$ n_e(x) \sim e^{\alpha x} \f$, where \f$ \alpha \f$ is the ionization coefficient. To accurately capture this growth, one should resolve the length scale \f$ 1/\alpha \f$. The code will therefore mark cells where `alpha(E) * dx > refine_adx` for refinement. To get an idea of the convergence of a model, one can run it several times, using for example:

    refine_adx = 2.0
    refine_adx = 1.0
    refine_adx = 0.5

## Other refinement parameters

Some of the other parameters that control refinement are shown below:

\snippet m_refine.f90 basic_refinement_settings

## Using a static grid for a specific region in the domain

The user can specify a specific region of the domain to have fixed grid spacing. The following flags are utilized for that purpose:

* `refine_regions_dr`: the minimum grid spacing in the region
* `refine_regions_rmin`: the starting boundary of the region
* `refine_regions_rmax`: the closing boundary of the region

The flag `refine_regions_rmin` should have two values if working in a two dimensional domain and three values if in a three dimensional domain. The same is the case for `refine_regions_rmax`. The values are separated by spaces. For example, having `refine_regions_rmin = 0.0 0.0` and  `refine_regions_rmax = 6.0e-4 5.0e-2` in a cylindrically symmetric simulation corresponds to the region bounded by the lines `r = 0.0`, `z = 0.0`, `r = 6.0e-4`, and `z = 5.0e-2`.


