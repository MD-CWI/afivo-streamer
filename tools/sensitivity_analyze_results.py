#!/usr/bin/env python3

import numpy as np
import argparse
import pandas as pd
import re

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
parser.add_argument('-num_bar_plot', type=int, default=0,
                    help='If >0, show N most important reactions for y[0]')
parser.add_argument('-bar_plot_xlim', type=float, nargs=2,
                    help='x-range for bar plots')
parser.add_argument('-figname', type=str,
                    help='Name of figure to save')
parser.add_argument('-label_replacement', nargs='+', type=str,
                    help='Pairs of (text to replace) and (replacement)')
parser.add_argument('-font', type=str,
                    help='Name of font to use')
args = parser.parse_args()

if args.num_bar_plot > 0 and len(args.y) > 1:
    raise ValueError('For bar plot, specify only one y variable')

logs = sorted(args.logs)

# Determine whether we analyse log files (pandas dataframes) or species
# amounts (text files)
if not all([x.endswith('amounts.txt') for x in logs]):
    logs_df = [pd.read_csv(f, sep=r'\s+') for f in args.logs]
    base_name = logs[0].replace('_log.txt', '')
else:
    # Make sure that the default argument is changed
    if args.y[0] == "sum(n_e)":
        args.y = ["e"]

    # Loading the species list and appending the time column to it
    base_name = logs[0].replace('_amounts.txt', '')
    with open(base_name + '_species.txt', 'r') as f:
        species_list = [x.strip() for x in f.readlines() if x.strip()]
    species_list.insert(0, "time")

    # Load amounts of species and create a dataframe of them so that the below
    # code wont break
    logs_df = [pd.DataFrame(np.loadtxt(f), columns=species_list)
               for f in args.logs]

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
reaction_ix = np.array([ix for ix in all_cases.keys() if ix != 0])

effect_magnitudes = np.zeros(len(reaction_ix))
derivatives_mean = np.zeros((len(reaction_ix), len(args.y)))
derivatives_sigma = np.zeros((len(reaction_ix), len(args.y)))

print(f'Using data at time t = {times[args.time_index]}\n')

# Here mu indicates the mean derivative, mustar the mean absolute derivative,
# and sigma the standard deviation in the derivative. All derivatives w.r.t. a
# factor f.
print(f'R{"#":<4} {"variable":15} {"mu":>15} {"mustar":>15} {"sigma":>15}')

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
    derivs_sigma = np.std(derivs_normalized, axis=0, ddof=1)

    for name, mu, mustar, sigma in zip(
            args.y, derivs_mean, derivs_meanabs, derivs_sigma):
        print(f'R{ix:<4} {name:15} {mu:15.8f} {mustar:15.8f} {sigma:15.8f}')

    effect_magnitudes[i] = derivs_meanabs.max()
    derivatives_mean[i] = derivs_mean
    derivatives_sigma[i] = derivs_sigma

print('\nReactions sorted by their overall importance:')
print(f'{"rank":<6} R{"#":<6} {"reaction_list":40} {"max(mustar)":15}')

# Load reaction names
with open(base_name + '_reactions.txt', 'r') as f:
    reactions_list = [x.strip() for x in f.readlines() if x.strip()]

ix_sort = np.argsort(effect_magnitudes)[::-1]
for n, i in enumerate(ix_sort):
    ix = reaction_ix[i]
    print(f'{n+1:<6} R{ix:<6} {reactions_list[ix-1]:40} ' +
          f'{effect_magnitudes[i]:<15.8f}')


def reaction_to_latex_format(label):
    # Convert a reaction to latex format

    tex = label

    if args.label_replacement is not None:
        for orig, repl in zip(args.label_replacement[::2],
                              args.label_replacement[1::2]):
            tex = tex.replace(orig, repl)

    # Replace +, - with superscript version ($^+$ or $^-$)
    pattern_plusmin = r'\b([+-])'
    tex = re.sub(pattern_plusmin, r'$^{\1}$', tex)

    # Match number following a letter (case of `letter+number`) and
    # replace with $_N$
    pattern_letter_number = r'([a-zA-Z])(\d+)'
    tex = re.sub(pattern_letter_number, r'\1$_{\2}$', tex)

    # Match number following '(' and replace with $^N$
    pattern_paren_number = r'\((\d+)'
    tex = re.sub(pattern_paren_number, r'($^{\1}$', tex)

    # Nicer arrows
    tex = tex.replace('->', r'$\to$')
    # Remove double $$
    tex = tex.replace('$$', '')

    return tex


if args.num_bar_plot > 0:
    # Make bar plot for args.y[0]
    import matplotlib.pyplot as plt

    if args.font is not None:
        # Set specific font
        plt.rcParams["font.family"].insert(0, args.font)

    N = args.num_bar_plot
    ixs = ix_sort[:N]
    r_ixs = reaction_ix[ixs]
    labels = [reaction_to_latex_format(reactions_list[i-1]) for i in r_ixs]
    colors = ['green' if x > 0 else 'red' for x in derivatives_mean[ixs, 0]]
    fig, ax = plt.subplots(1, 1, figsize=(4.5, 3.6), layout='constrained')
    thickness = 0.8
    x_values = np.arange(N, 0, -1)

    bars = ax.barh(x_values, np.abs(derivatives_mean[ixs, 0]),
                   tick_label=labels, color=colors, height=thickness)

    ax.set_xlabel(r'relative sensitivity')

    sigma_labels = [r'$\pm$ ' + f'{s:.1f}' for s in derivatives_sigma[ixs, 0]]
    ax.bar_label(bars, sigma_labels, padding=5, color='black')

    if args.bar_plot_xlim:
        ax.set_xlim(args.bar_plot_xlim)
    if args.figname:
        plt.savefig(args.figname, dpi=200, bbox_inches='tight')
        print(f'Saved {args.figname}')
    else:
        plt.show()
