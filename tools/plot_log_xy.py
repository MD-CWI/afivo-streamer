#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

p.add_argument("log_files", type=str, nargs='+', help="Input log file(s)")
p.add_argument("-x", type=str, default='time', help="Name of x variable")
p.add_argument("-y", type=str, nargs='+', default='[max(E)]',
               help="Name of y variables")
args = p.parse_args()

logs = [pd.read_csv(f, delim_whitespace=True) for f in args.log_files]
numbered_files = [f'{i}: {f}' for i, f in enumerate(args.log_files)]

fig, axes = plt.subplots(1, 1, constrained_layout=True)
fig.suptitle('\n'.join(numbered_files))

for i, log in enumerate(logs):
    for y in args.y:
        log.plot(args.x, y, ax=axes, label=f'{y}-{i}')

plt.xlabel(args.x)
plt.legend()
plt.show()
