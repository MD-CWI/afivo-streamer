# Python tools for input, output and analysis {#tools}

[TOC]

# Overview of included tools {#tools-overview}

The tools described below can be found in the `tools` folder. For all the Python-based tools, execute them with a `-h` argument to see how they can be used.

## afivo/tools/plot_raw_data.py

This script can be used to visualize and convert Silo files from Python, see @ref documentation/output_and_visualization.md

## bolsig_convert.py

This script can be used to convert the output from [Bolsig+](https://us.lxcat.net/solvers/BolsigPlus/) to a format that `afivo-streamer` can work with. The script also works with output generated from [lxcat](https://lxcat.net).

## plot_log_file.py

This script can be used to get a quick look at streamer properties (position, velocity, maximal field) from one or more log files. It also shows how to conveniently parse log files.

## plot_log_xy.py

This script can be used to plot two quantities from a log file against each other, for example the maximal field as a function of time.

## chemistry_visualize_rates.py

This script can be used to visualize (and study) the most important chemical reactions during a simulation.

## absorption_function.py

This script can be used to generate coefficients for Helmholtz-based photoionization approximations, which can be important when using gases other than air.

## Visit scripts

There are multiple scripts that can help to get data from Silo files. These can be executed in the following way:

      visit -nowin -cli -s <path to script> <script arguments>

* `visit_integrate_conditional.py`: volume integration, but with a conditional
* `visit_integrate_region.py`: volume integration over a region
* `visit_integrate_volume_overTime.py`: volume integration over time
* `visit_lineout.py`: perform a lineout
* `visit_plot_dataOverTime.py`: get data at a point over time
