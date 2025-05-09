#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
p.add_argument("log_file", type=str, nargs='+', help="Input log file(s)")
p.add_argument("-velocity_window", type=int,
               help="Compute average velocity over this window size")
args = p.parse_args()

logs = [pd.read_csv(f, sep=r'\s+') for f in args.log_file]
numbered_files = [f'{i}: {f}' for i, f in enumerate(args.log_file)]

fig, axes = plt.subplots(2, 2, constrained_layout=True)
fig.suptitle('\n'.join(numbered_files))

coord = 'x'
if 'y' in logs[0].columns:
    coord = 'y'
if 'z' in logs[0].columns:
    coord = 'z'

for i, log in enumerate(logs):
    # Compute smoothed velocity
    log['velocity'] = np.gradient(log[coord], log['time'])

    # Remove outliers in velocity
    v_low = log["velocity"].abs().quantile(0.1)
    v_hi = log["velocity"].abs().quantile(0.9)
    mask = np.logical_or(log['velocity'].abs() < 0.5 * v_low,
                         log['velocity'].abs() > 2 * v_hi)
    log['velocity'] = log['velocity'].mask(mask, np.nan)

    if args.velocity_window is not None:
        # Compute average
        w = args.velocity_window
        log['velocity'] = log['velocity'].rolling(window=w).mean()

    # Remove outliers in radius
    v_low = log['x.2'].abs().quantile(0.1)
    v_hi = log['x.2'].abs().quantile(0.9)
    mask = np.logical_or(log['x.2'].abs() < 0.5 * v_low,
                         log['x.2'].abs() > 2 * v_hi)
    log['x.2'] = log['x.2'].mask(mask, np.nan)

    log.plot('time', coord, ax=axes[0, 0], label=f'y-{i}')
    log.plot('time', 'velocity', ax=axes[0, 1], label=f'velocity-{i}')
    log.plot('time', 'max(E)', ax=axes[1, 0], label=f'max(E)-{i}')
    log.plot('time', 'x.2', ax=axes[1, 1], label=f'radius-{i}')
plt.legend()
plt.show()
