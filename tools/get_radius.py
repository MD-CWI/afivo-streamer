#!/usr/bin/env python3

# This script can be used to compute the full width half maximum (FWHM) of the
# time-integrated optical emission of a streamer simulation.
# It can be used for 3D and axisymmetric simulations.
#
# Author: Jannis Teunissen

import pandas as pd
import numpy as np
import argparse
import sys
sys.path.append('../afivo/tools')

from plot_raw_data import get_uniform_data, plot_uniform_data
from raw_reader import load_file

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Determine radius from log file and Silo data')
parser.add_argument('log', type=str, help='Log file')
parser.add_argument('silo_file', type=str,
                    help='Silo file containing optical emission')
parser.add_argument('-csv', type=str, default='fwhm.csv',
                    help='Output csv file')
parser.add_argument('-ix_range', type=int, nargs=2, help='Start and end index')
parser.add_argument('-max_radius', type=float, default=5e-3,
                    help='Maximal radius')
parser.add_argument('-axisymmetric', action='store_true',
                    help='Whether the simulation is axisymmetric')
parser.add_argument('-variable', type=str, default='N2_B3',
                    help='Variable corresponding to integrated light emission')
parser.add_argument('-min_pixels', type=int, default=4096,
                    help='Min. pixels if full domain was on a uniform grid')
parser.add_argument('-threshold', type=float, default=1e7,
                    help='Threshold for projected intensity to compute FWHM')
args = parser.parse_args()

df = pd.read_csv(args.log, delim_whitespace=True)

if args.ix_range is not None:
    df = df.iloc[args.ix_range[0]:args.ix_range[1]]

if args.axisymmetric:
    coords_names = ['x', 'y']
    project_dims = None
    abel_transform = True
else:
    coords_names = ['x', 'y', 'z']
    project_dims = [1]
    abel_transform = False

r_min = np.array([df[c].min() for c in coords_names]) - args.max_radius
r_max = np.array([df[c].max() for c in coords_names]) + args.max_radius

if args.axisymmetric:
    r_min[0] = 0.
else:
    # Drop projected dimension
    r_min = np.delete(r_min, project_dims)
    r_max = np.delete(r_max, project_dims)

grids, domain = load_file(args.silo_file, project_dims, args.variable)
values, coords = get_uniform_data(grids, domain, args.min_pixels,
                                  rmin=r_min, rmax=r_max,
                                  abel_transform=abel_transform)

if values.shape[0] < 100:
    raise ValueError('Fewer than 100 radial points, increase min_pixels')


def get_fwhm_3d(x, y):
    if y.max() < args.threshold:
        return 0.

    # Substract half of maximum
    hm = 0.5 * y.max()
    tmp = y.copy() - hm

    # Determine first and last index where tmp >= 0
    ix_valid = np.where(tmp >= 0)[0]
    i0, i1 = ix_valid[0], ix_valid[-1]

    # Interpolate linearly
    a, b = tmp[i0], tmp[i0-1]
    f0 = -b/(a-b)
    x0 = f0 * x[i0] + (1 - f0) * x[i0-1]

    c, d = tmp[i1], tmp[i1+1]
    f1 = -d/(c-d)
    x1 = f1 * x[i1] + (1 - f1) * x[i1+1]

    return x1 - x0


def get_fwhm_cyl(x, y):
    if y.max() < args.threshold:
        return 0.

    # Substract half of maximum
    hm = 0.5 * y.max()
    tmp = y.copy() - hm

    if tmp[0] < 0.:
        return 0.

    # Determine first index where tmp < 0
    i0 = np.where(tmp < 0)[0][0]

    # Interpolate linearly
    a, b = tmp[i0], tmp[i0-1]
    f0 = -b/(a-b)
    x0 = f0 * x[i0] + (1 - f0) * x[i0-1]

    return 2 * x0


if args.axisymmetric:
    fwhm = np.array([get_fwhm_cyl(coords[0], values[:, i])
                     for i in range(values.shape[1])])
else:
    fwhm = np.array([get_fwhm_3d(coords[0], values[:, i])
                     for i in range(values.shape[1])])

last_dim = coords_names[-1]
df['fwhm'] = np.interp(df[last_dim], coords[1], fwhm)

columns = ['time', 'fwhm', last_dim]
df[columns].to_csv(args.csv, index=False)
print('Saved ' + args.csv)
