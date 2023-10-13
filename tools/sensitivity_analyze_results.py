#!/usr/bin/env python3

# TODO:
# - Have option to read in chemistry output (species amounts)

import numpy as np
import matplotlib.pyplot as plt
import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='''Script to analyze results from a sensitivity study
    Authors: Hemaditya Malla, Jannis Teunissen''')
parser.add_argument('logs', type=str, nargs='+',
                    help='Log files')
parser.add_argument('-y', type=str, nargs='+', default=["sum(n_e)"],
                    help='y variable of the files to analyse')
parser.add_argument('-time_index', type=int, default=-1,
                    help='Which time index to consider')
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
        all_cases[ix].append((fac, df))
    else:
        all_cases[ix] = [(fac, df)]

if 0 not in all_cases:
    raise ValueError('Base case not found (..._ix0000_...)')

base_case = all_cases[0][0][1]
times = np.array(base_case['time'])
reaction_ix = [ix for ix in all_cases.keys() if ix != 0]

effect_magnitudes = np.zeros(len(reaction_ix))

print(f'Using data at time t = {times[args.time_index]}')

# Here mu indicates the mean derivative, mustar the mean absolute derivative,
# and sigma the standard deviation in the derivative. All derivatives w.r.t. a
# factor f.
print(f'R{"#":<4} {"variable":15} {"mu":>15} {"mustar":>15} {"sigma":>15}')

for i, ix in enumerate(reaction_ix):
    values = np.array([df[args.y] for _, df in all_cases[ix]])
    factors = np.array([f for f, _ in all_cases[ix]])

    # Get values at time index
    values = values[:, args.time_index]

    # Get base values
    base_values = np.array(base_case[args.y])[args.time_index]
    base_values = np.tile(base_values, [len(factors), 1])

    # Compute dg/df ~ (g(f) - g(1))/(f-1), where f is the factor centered
    # around 1. E.g., if a factor is 1.1, the delta is 0.1
    deltas = np.tile(factors - 1, [len(args.y), 1]).T
    derivs = (values - base_values) / deltas

    derivs_normalized = derivs / base_values
    derivs_mean = np.mean(derivs_normalized, axis=0)
    derivs_meanabs = np.mean(np.abs(derivs_normalized), axis=0)
    derivs_sigma = np.std(derivs_normalized, axis=0)

    for name, mu, mustar, sigma in zip(
            args.y, derivs_mean, derivs_meanabs, derivs_sigma):
        print(f'R{ix:<4} {name:15} {mu:15.8f} {mustar:15.8f} {sigma:15.8f}')

    effect_magnitudes[i] = derivs_meanabs.max()

print('\nReactions sorted by their overall importance:')

# Load reaction names
base_name = logs[0].replace('_log.txt', '')
with open(base_name + '_reactions.txt', 'r') as f:
    reactions_list = [x.strip() for x in f.readlines() if x.strip()]

ix_sort = np.argsort(effect_magnitudes)[::-1]
for i in ix_sort:
    ix = reaction_ix[i]
    print(f'R{ix:<4} {reactions_list[ix-1]:40} {effect_magnitudes[i]:<15.8f}')
