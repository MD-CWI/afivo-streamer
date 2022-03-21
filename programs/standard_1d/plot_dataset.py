#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import argparse

p = argparse.ArgumentParser()
p.add_argument('-dataset', type=str, default='dataset.npz',
               help='Dataset filename')
args = p.parse_args()

f = np.load(args.dataset)
X = f['X']
ix_list = f['ix_list']
dx_list = f['dx_list']

fig, ax = plt.subplots(3, sharex=True)

for x, ix in zip(X, ix_list):
    ax[0].plot(x[0])
    ax[1].plot(x[1])
    ax[2].plot(x[2])

plt.show()
