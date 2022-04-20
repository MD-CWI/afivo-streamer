#!/usr/bin/env python3

#Author: Hemaditya Malla

import numpy as np
import matplotlib.pyplot as plt
import re

#TODO: Add an argparse framework to read the filename, specie/reaction of interest



fname = "./output/streamer_cyl_electrode_"
suff = ".txt"
rmfile = fname+"rate_matrix"+suff
chem_deetsfile = fname+"chemistry_details"+suff
rates = fname+"rates"+suff


#Making a list of species and reactions
species_list = list()
reactions_list = list()
with open(chem_deetsfile, "r") as f:
    fr = f.read()

fr = re.sub(r'\r\n', '\n', fr)  # Fix newlines
lines = fr.splitlines()
for il, line in enumerate(lines[2:]):
    if line == "":
        break
    species_list.append(line.strip())
print(species_list, len(species_list))
for line in lines[il+5:]:
    if line == " ":
        break
    reactions_list.append(line.strip())
print(reactions_list, len(reactions_list))


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

