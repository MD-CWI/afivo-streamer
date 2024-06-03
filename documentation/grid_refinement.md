# Grid refinement {#grid-refinement}

## Parameter refine_adx
Grid refinement is controlled by several parameters. Usually the most import one is

    # Refine if alpha*dx is larger than this value:
    refine_adx = 1.0

The idea behind this parameter is that the electron density will locally approximately increase like \f$ n_e(x) \sim e^{\alpha x} \f$, where \f$ \alpha \f$ is the ionization coefficient. To accurately capture this growth, one should resolve the length scale \f$ 1/\alpha \f$. The code will therefore mark cells where `refine_adx * alpha(E) > dx` for refinement. To get an idea of the convergence of a model, one can run it several times, using for example:

    refine_adx = 2.0
    refine_adx = 1.0
    refine_adx = 0.5

## Other refinement parameters

Some of the other parameters that control refinement are shown below:

\snippet m_refine.f90 basic_refinement_settings

