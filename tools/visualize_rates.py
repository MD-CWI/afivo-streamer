#!/usr/bin/env python3

#Author: Hemaditya Malla

import numpy as np
import matplotlib.pyplot as plt
import argparse

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
p.add_argument('rates_file', type=str, help='File with reaction rates')
p.add_argument('-soi', type=str, help='Species of interest')
p.add_argument('-list_species', action='store_true', help='List species')
p.add_argument('-list_reactions', action='store_true', help='List reactions')
p.add_argument('-plot_all', action='store_true',
               help='Plot all reaction rates together')
args = p.parse_args()

# Assume the other files are in the same folder
base_name = args.rates_file.replace('_rates.txt', '')

with open(base_name + '_species.txt', 'r') as f:
    species_list = [x.strip() for x in f.readlines() if x.strip()]
with open(base_name + '_reactions.txt', 'r') as f:
    reactions_list = [x.strip() for x in f.readlines() if x.strip()]

rate_matrix = np.loadtxt(base_name + '_rate_matrix.txt')
n_reactions, n_species = rate_matrix.shape

tmp = np.loadtxt(base_name + '_rates.txt')
time = tmp[:, 0]
rates = tmp[:, 1:]

if args.list_species:
    for i, name in enumerate(species_list):
        print(f'{i:4} {name}')

if args.list_reactions:
    for i, name in enumerate(reactions_list):
        print(f'{i:4} {name}')

if args.plot_all:
    # Visualize all the rates vs. time
    plt.figure()
    for i, rxn in enumerate(reactions_list):
        plt.plot(time, rates[:, i], label=rxn)
    plt.legend()
    plt.tight_layout()
    plt.show()

if args.soi:
    # Visualizing the source and sink reaction rates for a given specie
    sidx = species_list.index(args.soi)
    fig, ax = plt.subplots(2, figsize=(8, 8), sharex=True)

    srce_idx = np.where(rate_matrix[:, sidx] > 0)[0]
    sink_idx = np.where(rate_matrix[:, sidx] < 0)[0]
    titles = ['Source', 'Sink']

    for i, (ix, text) in enumerate(zip([srce_idx, sink_idx], titles)):
        amount = rate_matrix[ix, sidx] * rates[:, ix]
        frac = amount[-1]/amount[-1].sum()

        for j, idx in enumerate(ix):
            ax[i].plot(time, amount[:, j], label=reactions_list[idx] +
                       f' ({100*frac[j]:.2f}%)')

        ax[i].set_title(text + ' reactions')
        ax[i].set_xlabel('Time (s)')
        ax[i].set_ylabel('Production (#)')
        ax[i].legend()

    fig.suptitle(f'{len(srce_idx)+len(sink_idx)} of {n_reactions}'
                 f' influence {args.soi}')

    plt.tight_layout()
    plt.show()
