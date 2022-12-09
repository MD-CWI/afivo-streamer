#!/usr/bin/env python3

# Authors: Hemaditya Malla, Jannis Teunissen

import numpy as np
import matplotlib.pyplot as plt
import argparse

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
p.add_argument('rates_file', type=str, help='File with reaction rates')
p.add_argument('-soi', type=str, default=[], nargs='+',
               help='Species of interest')
p.add_argument('-list_species', action='store_true', help='List species')
p.add_argument('-list_reactions', action='store_true', help='List reactions')
p.add_argument('-plot_all', action='store_true',
               help='Plot all reaction rates together')
p.add_argument('-time_interval', nargs=2, type=float,
               help='Time interval over which to analyse the reactions (s)')
args = p.parse_args()

# Assume the other files are in the same folder
base_name = args.rates_file.replace('_rates.txt', '')

with open(base_name + '_species.txt', 'r') as f:
    species_list = [x.strip() for x in f.readlines() if x.strip()]
with open(base_name + '_reactions.txt', 'r') as f:
    reactions_list = [x.strip() for x in f.readlines() if x.strip()]

stoich_matrix = np.loadtxt(base_name + '_stoich_matrix.txt')
n_species, n_reactions = stoich_matrix.shape

tmp = np.loadtxt(base_name + '_rates.txt')
time = tmp[:, 0]
rates = tmp[:, 1:]

if args.time_interval is not None:
    # Only consider chemistry within the given interval
    t1_idx = np.where(time >= args.time_interval[0])[0][0]
    # +1 for array slicing, but do not go beyond array size
    t2_idx = np.where(time <= args.time_interval[1])[0][-1] + 1
    t2_idx = min(t2_idx, len(time))

    time = time[t1_idx:t2_idx]
    rates = rates[t1_idx:t2_idx]

if args.list_species:
    for i, name in enumerate(species_list):
        print(f'{i:4} {name}')

if args.list_reactions:
    for i, name in enumerate(reactions_list):
        print(f'{i:4} {name}')

if args.plot_all:
    # Visualize all the rates vs. time
    plt.figure(figsize=(8, 8))

    ix = np.argsort(rates[-1, :])[::-1]
    sum_of_rates = rates[-1, :].sum()

    for i in ix:
        plt.plot(time, rates[:, i], label=reactions_list[i] +
                 f' ({100*rates[-1, i]/sum_of_rates:.2f}%)')
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")

for soi in args.soi:
    # Visualize the source and sink reactions for a given specie
    sidx = species_list.index(soi)
    fig, ax = plt.subplots(3, figsize=(8, 12), sharex=True)

    srce_idx = np.where(stoich_matrix[sidx, :] > 0)[0]
    sink_idx = np.where(stoich_matrix[sidx, :] < 0)[0]
    titles = ['Source', 'Sink']

    for i, (ix, text) in enumerate(zip([srce_idx, sink_idx], titles)):
        amount = stoich_matrix[sidx, ix] * rates[:, ix]
        frac = amount[-1]/amount[-1].sum()

        for j, idx in enumerate(ix):
            ax[i].plot(time, amount[:, j], label=reactions_list[idx] +
                       f' ({100*frac[j]:.2f}%)')

        ax[i].set_title(text + ' reactions')
        ax[i].set_xlabel('Time (s)')
        ax[i].set_ylabel('Production (#)')
        ax[i].legend()

    gross_prod = np.dot(rates[:, srce_idx], stoich_matrix[sidx, srce_idx])
    net_prod = np.dot(rates, stoich_matrix[sidx])
    ax[2].plot(time, gross_prod, label='gross production')
    ax[2].plot(time, net_prod, label='net production')
    ax[2].set_xlabel('Time (s)')
    ax[2].set_ylabel('Production (#)')
    ax[2].legend()
    fig.suptitle(f'{len(srce_idx)+len(sink_idx)} of {n_reactions}'
                 f' influence {soi}')

if plt.get_fignums():
    plt.tight_layout()
    plt.show()
