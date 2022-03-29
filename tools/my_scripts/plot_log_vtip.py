#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

p.add_argument("log_file", type=str, nargs='+', help="Input log file(s)")
args = p.parse_args()

logs = [pd.read_csv(f, delim_whitespace=True) for f in args.log_file]
numbered_files = [f'{i}: {f}' for i, f in enumerate(args.log_file)]

fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle('\n'.join(numbered_files))

for i, log in enumerate(logs):
    # Determine dt for computing derivatives
    dt = log['time'][1] - log['time'][0]
    
    # Compute smoothed velocity
    log['velocity'] = -savgol_filter(log['y.3'], 11, 2, deriv=1, delta=dt)

    # Remove outliers in velocity
    v_low = log["velocity"].abs().quantile(0.1)
    v_hi = log["velocity"].abs().quantile(0.9)
    mask = np.logical_or(log['velocity'].abs() < 0.5 * v_low,
                         log['velocity'].abs() > 2 * v_hi)
    log['velocity'].mask(mask, np.nan, inplace=True)
    
    log.plot('y.3', 'velocity', ax=axes, label=f'velocity-{i}')

plt.xlabel('y.3')
plt.ylabel('Velocity')
plt.grid(linestyle='-.',alpha=0.7)
plt.legend()
plt.show()
