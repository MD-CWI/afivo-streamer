#!/usr/bin/env python3

# This script can be used to compute the radius of a streamer discharge from the
# on-axis electric field profile E(z)
#
# Author: Jannis Teunissen

import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Determine radius from on-axis electric field profile')
parser.add_argument('infile', type=str, help='Input file')
parser.add_argument('-z_column', type=int, default=0,
                    help='Index of column with z-coordinate')
parser.add_argument('-E_column', type=int, default=1,
                    help='Index of column with electric field E(z)')
parser.add_argument('-E_bg', type=float,
                    help='Background electric field')
parser.add_argument('-factor', type=float, default=0.5,
                    help='Fit until value is below max(E) * factor')
parser.add_argument('-sep', type=str, default=r'\s+',
                    help='Delimiter for data')
parser.add_argument('-skiprows', type=int, default=0,
                    help='Skip this many rows when reading the data')
parser.add_argument('-charge_layer_width', type=float,
                    help='Manually set width of the charge layer')
parser.add_argument('-compact_output', action='store_true',
                    help='Compact output for processing with other scripts')
args = parser.parse_args()

# Load data
df = pd.read_csv(args.infile, sep=r'\s+',
                 skiprows=args.skiprows, dtype=float, header=0,
                 names=['z', 'E'], usecols=[args.z_column, args.E_column])

z = df['z'].values
E = df['E'].values

# Locate maximum and ensure it is positive
i_max = np.argmax(np.abs(E))

if E[i_max] < 0:
    E = -E

E_max = E[i_max]

# Get background field
E_bg = args.E_bg

if E_bg is None:
    E_bg = np.median(E)
    if not args.compact_output:
        print(f'Estimated background field: {E_bg:.3e}')

# Locate distance where value has dropped to factor * max(E)
distance_pos = np.argmax(E[i_max:] < args.factor * E_max)
distance_neg = np.argmax(np.flip(E[:i_max+1]) < args.factor * E_max)

direction = np.sign(distance_pos - distance_neg)

if direction > 0:
    z = z[i_max:i_max+distance_pos+1] - z[i_max]
    E = E[i_max:i_max+distance_pos+1]
else:
    z = z[i_max] - np.flip(z[i_max-distance_neg:i_max+1])
    E = np.flip(E[i_max-distance_neg:i_max+1])


def fit_func(z, R, E_max):
    return E_bg + (E_max - E_bg) * (z/R + 1)**-2


R_guess = (args.factor + args.factor**0.5)/(1 - args.factor) * z[-1]

if args.charge_layer_width is None:
    # Start fit from location where the derivative starts to decrease
    dEdz = np.abs(np.gradient(E))
    n_skip = np.argmax(np.diff(dEdz) < 0)
    if not args.compact_output:
        print(f'Estimated charge layer width: {z[n_skip]:.3e}')
else:
    n_skip = np.argmax(z - args.charge_layer_width >= 0)

popt, pcov = curve_fit(fit_func, z[n_skip:], E[n_skip:],
                       p0=[R_guess, E[n_skip]])

if args.compact_output:
    print(f'{popt[0]:.3e}')
else:
    print(f'Fitted radius: {popt[0]:.3e}')
    print(f'Fitted E_max:  {popt[1]:.3e}')

    fig, ax = plt.subplots()

    ax.plot(z, E, label='data')
    ax.plot(z[n_skip:], E[n_skip:], ls='--', label='fit range')
    ax.plot(z, fit_func(z, *popt), label='fit')
    ax.legend()

    plt.show()
