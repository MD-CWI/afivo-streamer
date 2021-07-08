#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

p.add_argument("log_file", type=str, nargs='+', help="Input log file(s)")
p.add_argument("-savgol_width", type=int, default=5, help="Width of savgol filter")
p.add_argument("-savgol_order", type=int, default=2, help="Order of savgol filter")
args = p.parse_args()

logs = [pd.read_csv(f, delim_whitespace=True) for f in args.log_file]
numbered_files = [f'{i}: {f}' for i, f in enumerate(args.log_file)]

fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle('\n'.join(numbered_files))

for i, log in enumerate(logs):
    # Determine dt for computing derivatives
    dt = log['time'][1] - log['time'][0]

    # Check whether the data is 1d, 2d or 3d and compute the derivatives of x,
    # y, z coordinates
    coords = ['vx']
    log['vx'] = np.abs(np.gradient(log['x'], log['time']))
    log['vx_savgol'] = np.abs(
        savgol_filter(log['x'], args.savgol_width, args.savgol_order,
                      deriv=1, delta=dt))

    if 'y' in log.columns:
        log['vy'] = np.abs(np.gradient(log['y'], log['time']))
        log['vy_savgol'] = np.abs(
            savgol_filter(log['y'], args.savgol_width, args.savgol_order,
                          deriv=1, delta=dt))
        coords += ['vy']
    if 'z' in log.columns:
        log['vz'] = np.abs(np.gradient(log['z'], log['time']))
        log['vz_savgol'] = np.abs(
            savgol_filter(log['z'], args.savgol_width, args.savgol_order,
                          deriv=1, delta=dt))
        coords += ['vz']

    # Compute the norm of the x, y, z velocities
    log['v_gradient'] = np.linalg.norm(log[coords], axis=1)
    coords_savgol = [c + '_savgol' for c in coords]
    log['v_savgol'] = np.linalg.norm(log[coords_savgol], axis=1)

    log.plot('time', 'v', ax=axes, label=f'v_log-{i}')
    log.plot('time', 'v_gradient', ax=axes, label=f'v_gradient-{i}')
    log.plot('time', 'v_savgol', ax=axes, label=f'v_savgol-{i}')

plt.legend()
plt.show()
