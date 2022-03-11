#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse

p = argparse.ArgumentParser()
p.add_argument('files', type=str, nargs='+',
               help='Input files (e.g. *.npz)')
p.add_argument('-w', type=int, default=128,
               help='Width of samples (in cells)')
p.add_argument('-min_rhs', type=float, default=1e9,
               help='Required amplitude of rhs')
p.add_argument('-output', type=str, default='dataset.npz',
               help='Where to store the output')
args = p.parse_args()

files = sorted(args.files)

variables = ['e', 'electric_fld', 'rhs']
X = np.zeros((len(files), len(variables), args.w))
dx_list = np.zeros(len(files))
ix_list = np.zeros(len(files), dtype=int)
n_samples = 0


for fname in files:
    f = np.load(fname)

    # Locate maximum of the rhs, which is -rho/eps0
    i_origin = np.argmax(f['rhs'])

    if f['rhs'][i_origin] < args.min_rhs:
        continue

    dx_list[n_samples] = f['dr'][0]
    ix_list[n_samples] = int(fname[-10:-4])

    i0 = i_origin - args.w//2
    i1 = i0 + args.w

    for n, name in enumerate(variables):
        X[n_samples, n, :] = f[name][i0:i1]
    n_samples += 1

np.savez(args.output,
         X=X[:n_samples],
         dx_list=dx_list[:n_samples],
         ix_list=ix_list[:n_samples])

# x_min = f['r_min'][0]
# x_max = f['r_max'][0]
# nx = f['nx'][0]
# x = np.linspace(x_min+0.5*dx, x_max-0.5*dx, nx)
# plt.plot(x, f['rhs'], '--')
# plt.plot(x[i0:i1], f['rhs'][i0:i1])
# plt.show()
