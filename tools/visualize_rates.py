#!/usr/bin/env python3

#Author: Hemaditya Malla

import numpy as np
import matplotlib.pyplot as plt
import argparse

p = argparse.ArgumentParser()
p.add_argument('base_name', type=str, help='Base of input filename')
args = p.parse_args()

rmfile = args.base_name + "_rate_matrix.txt"
rates = args.base_name + "_rates.txt"

with open(args.base_name + "_species.txt", 'r') as f:
    species_list = [x.strip() for x in f.readlines()]
with open(args.base_name + "_reactions.txt", 'r') as f:
    reactions_list = [x.strip() for x in f.readlines()]

#Reading the matrix
# dim 0 -- reactions
# dim 1 -- species
rate_matrix = np.loadtxt(rmfile)
#Readting the rates
rate_data = np.loadtxt(rates)
time = rate_data[:,0]


#Visualizing all the rates vs. time
plt.figure(1)
for i,rxn in enumerate(reactions_list):
        plt.plot(time, rate_data[:,i],label=rxn)
plt.legend()
#Visualizing the source and sink reaction rates for a given specie
soi = "O4_plus" #Taking the electron as an examplee
sidx = species_list.index(soi)
fig, ax = plt.subplots(2, sharex = True)
specie_reaction_idx = np.nonzero(rate_matrix[:,sidx])[0] + 1
print(specie_reaction_idx)
for idx in specie_reaction_idx:
    multiplicity = rate_matrix[idx-1, sidx]
    print(multiplicity)
    if multiplicity < 0: #Loss reactions
        ax[0].plot(time, multiplicity*rate_data[:,idx], label=reactions_list[idx-1])
    elif multiplicity > 0: #Source reactions
        ax[1].plot(time, multiplicity*rate_data[:,idx], label=reactions_list[idx-1])
    else:
        continue
fig.suptitle(str(len(specie_reaction_idx))+" of "+ str(len(reactions_list))+" influence " + soi)
ax[0].title.set_text("Sink reactions")
ax[1].title.set_text("Source reactions")
ax[0].legend()
ax[1].legend()
plt.show()

