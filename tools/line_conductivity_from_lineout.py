#!/usr/bin/env python3

# Compute line conductivity from a radial lineout

import numpy as np
import argparse
from scipy.integrate import trapezoid

# Get and parse the command line arguments
pr = argparse.ArgumentParser(
    description='''Compute line conductivity from a radial lineout''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
pr.add_argument('lineout', type=str,
                help='File with radial lineout of conductivity')
pr.add_argument('-conductivity_threshold', type=float, default=10e-2,
                help='Threshold for conductivity for determining radius')
args = pr.parse_args()

lineout = np.genfromtxt(args.lineout).T

# Define radius as location where conductivity is 1% of its maximum
i = np.argmax(lineout[1] < args.conductivity_threshold * lineout[1].max())
radius = lineout[0][i]

sigma_integral = trapezoid(2 * np.pi * lineout[0] * lineout[1], lineout[0])
sigma_radial_average = sigma_integral/(np.pi * radius**2)
print(radius, sigma_integral, sigma_radial_average)
