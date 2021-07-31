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
    # Compute smoothed delta_phi
    log['delta_phi_smooth'] = savgol_filter(log['delta_phi'], 21, 2)

    # Remove outliers in velocity
    phi_low = log['delta_phi'].abs().quantile(0.1)
    phi_hi = log['delta_phi'].abs().quantile(0.9)
    mask = np.logical_or(log['delta_phi'].abs() < 0.5 * phi_low,
                         log['delta_phi'].abs() > 2 * phi_hi)
    log['delta_phi'].mask(mask, np.nan, inplace=True)

    log.plot('y', 'delta_phi', ax=axes, label=f'delta_phi-{i}')
    log.plot('y', 'delta_phi_smooth', ax=axes, linestyle='--', label=f'delta_phi_smooth-{i}')
    
plt.xlabel('y')
plt.ylabel('delta_phi')
plt.grid(linestyle='-.',alpha=0.7)
plt.legend()
plt.show()
