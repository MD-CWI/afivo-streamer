#!/usr/bin/env python3

# TODO:
# - Remove zero effect reactions while making the bar plots
# - Make the bar plots prettier

import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='''Script to analyze results from a sensitivity study
    Authors: Hemaditya Malla, Jannis Teunissen''')
parser.add_argument('logs', type=str, nargs='+',
                    help='Log/amounts files')
parser.add_argument('-y', type=str, nargs='+', default=["sum(n_e)"],
                    help='Variables in the log files to compare')
parser.add_argument('-time_index', type=int, default=-1,
                    help='Which time index in the log files to consider')
parser.add_argument('-bar_plot', action='store_true',
               help='Make a bar plot for each variable y')
args = parser.parse_args()

logs = sorted(args.logs)

# Determine whether we analyse log files(pandas dataframe) or species amounts(txt files)
if not all([x.endswith('amounts.txt') for x in logs]):
    logs_df = [pd.read_csv(f, delim_whitespace=True) for f in args.logs]
    base_name = logs[0].replace('_log.txt', '')
else:
    # Can we use the below variable elsewhere below?
    analyse_chemistry = True
    # Make sure that the default argument is changed
    if args.y[0] == "sum(n_e)":
        args.y = ["e"]

    # Loading the species list and appending the time column to it
    base_name = logs[0].replace('_amounts.txt', '')
    with open(base_name + '_species.txt', 'r') as f:
        species_list = [x.strip() for x in f.readlines() if x.strip()]
    species_list.insert(0, "time")

    # Load amounts of species and create a dataframe of them so that the below code wont break
    logs_df = [pd.DataFrame(np.loadtxt(f), columns=species_list) for f in args.logs]


log_sizes = np.array([len(df) for df in logs_df])
max_size, min_size = log_sizes.max(), log_sizes.min()

if max_size > min_size:
    print(f'Warning: logs have different size, truncating to {min_size} rows')
    logs_df = [df.head(min_size) for df in logs_df]

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

print(f"Reaction test: {len(reaction_ix)}")
effect_magnitudes = np.zeros(len(reaction_ix))

print(times)
print(f'Using data at time t = {times[args.time_index]}\n')

# Here mu indicates the mean derivative, mustar the mean absolute derivative,
# and sigma the standard deviation in the derivative. All derivatives w.r.t. a
# factor f.
print(f'R{"#":<4} {"variable":15} {"mu":>15} {"mustar":>15} {"sigma":>15}')

stats = np.zeros((len(reaction_ix), len(args.y), 3))
for i, ix in enumerate(reaction_ix):
    # print("test: ", args.y)
    # test = [df[args.y] for _, df in all_cases[ix]]
    values = np.array([df[args.y].to_numpy() for _, df in all_cases[ix]])
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

    stats[i] = np.array([derivs_mean, derivs_meanabs, derivs_sigma]).T
    for name, mu, mustar, sigma in zip(
            args.y, derivs_mean, derivs_meanabs, derivs_sigma):
        print(f'R{ix:<4} {name:15} {mu:15.8f} {mustar:15.8f} {sigma:15.8f}')

    effect_magnitudes[i] = derivs_meanabs.max()

print('\nReactions sorted by their overall importance:')
print(f'{"rank":<6} R{"#":<6} {"reaction_list":40} {"max(mustar)":15}')

# Load reaction names
# base_name = logs[0].replace('_log.txt', '')
with open(base_name + '_reactions.txt', 'r') as f:
    reactions_list = [x.strip() for x in f.readlines() if x.strip()]

ix_sort = np.argsort(effect_magnitudes)[::-1]
used_reactions_list = []
for n, i in enumerate(ix_sort):
    ix = reaction_ix[i]
    used_reactions_list.append(reactions_list[ix-1])
    print(f'{n+1:<6} R{ix:<6} {reactions_list[ix-1]:40} ' +
          f'{effect_magnitudes[i]:<15.8f}')


# Plotting stuff
if args.bar_plot:
    for i, iy in enumerate(args.y):
        fig, ax = plt.subplots(1,1, figsize=(5, 20))
        fig.suptitle(f"Sensitivity mean and sigma for $\Delta$ {iy}")
        ax.barh(np.arange(stats.shape[0]), stats[:, i, 0],
                xerr=stats[:,i,2], align='center', 
                tick_label=used_reactions_list)
        plt.tight_layout()
    plt.show()
