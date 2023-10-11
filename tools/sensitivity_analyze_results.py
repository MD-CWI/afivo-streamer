#!/usr/bin/env python3

# TODO:
# - Add time argument (compare at specific time in log files)
# - Have option to read in chemistry output (species amounts)
# - Sort reactions by their effect on all quantities of interest
# - Option to show e.g. box plots
# - Show information about derivatives

import numpy as np
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='''Script to analyze results from a sensitivity study
    Authors: Hemaditya Malla, Jannis Teunissen''')
parser.add_argument('logs', type=str, nargs='+',
                    help='Log files')
parser.add_argument('-y', type=str, nargs='+', default=["sum(n_e)"],
                    help='y variable of the files to analyse')
args = parser.parse_args()

logs = sorted(args.logs)
logs_df = [pd.read_csv(f, delim_whitespace=True) for f in args.logs]

all_cases = {}

for log, df in zip(logs, logs_df):
    # Extract reaction index and factor from file name
    tmp = log.split('_')
    ix, fac = int(tmp[-3][2:]), float(tmp[-2][3:])

    # Store all dataframes for a reaction index in a list
    if ix in all_cases:
        all_cases[ix].append(df)
    else:
        all_cases[ix] = [df]

if 0 not in all_cases:
    raise ValueError('Base case not found (..._ix0000_...)')

base_case = all_cases[0][0]
reaction_ix = [ix for ix in all_cases.keys() if ix != 0]

# Include base case in all cases
for ix in reaction_ix:
    all_cases[ix].append(base_case)

print(f'{"  ix":4} {"quantity":15} {"     std/mean":15}')

for y in args.y:
    for ix in reaction_ix:
        values = np.array([df[y] for df in all_cases[ix]])
        means = np.mean(values[:, -1], axis=0)
        std = np.std(values[:, -1], axis=0)
        print(f'{ix:4} {y:15} {std/means:15.8f}')
