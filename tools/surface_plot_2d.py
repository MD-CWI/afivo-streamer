#!/usr/bin/env python3

# Plot surface variables from a 2D simulation
#
# Authors: Xiaoran Li, Jannis Teunissen

import numpy as np
import matplotlib.pyplot as plt
import argparse

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
p.add_argument('npz', type=str, help='Surface npz file')
args = p.parse_args()

dimnames = ['x', 'y']

X = np.load(args.npz)
n = X['n_surfaces']
nc = X['n_cell']

# Assume the surface is flat, so that all normal_dim are the same
normal_dim = X['normal_dim'][0] - 1
dim = 1 - normal_dim

r = X['r'][dim]

# For convenience, determine the grid spacing for each cell
dr = np.repeat(X['dr'][0], nc)

# Get indices that sort values along the surface coordinate
ix = np.argsort(r)

fig, ax = plt.subplots(3)
ax[0].plot(r[ix], X['surf_dens'][ix])
ax[0].set_xlabel(dimnames[dim] + ' (m)')
ax[0].set_ylabel('surface charge')

ax[1].plot(r[ix], X['photon_flux'][ix])
ax[1].set_xlabel(dimnames[dim] + ' (m)')
ax[1].set_ylabel('photon flux')

ax[2].plot(r[ix], dr[ix])
ax[2].set_xlabel(dimnames[dim] + ' (m)')
ax[2].set_ylabel('grid spacing')
plt.show()
