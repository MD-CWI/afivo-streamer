#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from raw_reader import get_raw_data, map_grid_data_to
import argparse
import os
import tempfile
import subprocess

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
p.add_argument('input_file', type=str, help='Input silo or raw file')
p.add_argument('-variable', type=str, help='Which variable to plot')
p.add_argument('-min_pixels', type=int, default=512,
               help='Min. pixels for any dimension')
p.add_argument('-project_dims', type=int, nargs='+', choices=[0, 1, 2],
               help='Project (integrate) along dimension(s)')
p.add_argument('-save_plot', type=str,
               help='Save the plot into this file')
p.add_argument('-axisymmetric', action='store_true',
               help='Assume axisymmetric geometry when averaging')
p.add_argument('-abel_transform', action='store_true',
               help='Perform forward Abel transform (implies axisymmetric)')
p.add_argument('-vmin', type=float, help='Minimum intensity in plot')
p.add_argument('-vmax', type=float, help='Maximum intensity in plot')
p.add_argument('-xlim', type=float, nargs=2,
               help='Plot range in x direction')
p.add_argument('-ylim', type=float, nargs=2,
               help='Plot range in y direction')
p.add_argument('-hide_axes', action='store_true',
               help='Hide axes etc.')
p.add_argument('-save_npz', type=str,
               help='Save the data into npz file and do not plot')
p.add_argument('-cmap', type=str, default='plasma',
               help='Use this colormap')
p.add_argument('-silo_to_raw', type=str, default='./silo_to_raw',
               help='Path to silo_to_raw converter')
args = p.parse_args()

if args.abel_transform:
    args.axisymmetric = True

extension = os.path.splitext(args.input_file)[1]

if extension == ".silo":
    if not args.variable:
        raise ValueError('Specify variable to plot using -variable')

    if not os.path.exists(args.silo_to_raw):
        raise ValueError('Could not find silo_to_raw, specify -silo_to_raw')

    # Convert silo to raw
    with tempfile.NamedTemporaryFile() as fp:
        res = subprocess.call([args.silo_to_raw, args.input_file,
                               args.variable, fp.name])
        grids, domain = get_raw_data(fp.name, args.project_dims)
else:
    grids, domain = get_raw_data(args.input_file, args.project_dims)

# Grid size nx should be of form 2^k * domain['n_cells_coarse']
ratio = args.min_pixels / domain['n_cells_coarse'].min()
ratio = 2**(np.ceil(np.log2(ratio)).astype(int))
nx = domain['n_cells_coarse'] * ratio

# Grid spacing on uniform grid
dr = (domain['r_max'] - domain['r_min']) / nx

print(f'Resolution:   {nx}')
print(f'Grid spacing: {dr}')

# List of coordinates
x = [np.linspace(a, b, n) for a, b, n in
     zip(domain['r_min']+0.5*dr, domain['r_max']-0.5*dr, nx)]

uniform_data = np.zeros(nx)

# Map each grid to the uniformly spaced output
for g in grids:
    grid_data, ix_lo, ix_hi = map_grid_data_to(g, domain['r_min'],
                                               dr, args.axisymmetric)

    # Index range on the output array
    g_ix = tuple([np.s_[i:j] for (i, j) in zip(ix_lo, ix_hi)])
    uniform_data[g_ix] += grid_data


if args.abel_transform:
    from abel.hansenlaw import hansenlaw_transform
    if args.xlim and args.xlim[0] > 0:
        raise ValueError('Abel transform requires xlim[0] to be 0.')
    tmp = hansenlaw_transform(uniform_data.T, dr[0], direction='forward')
    uniform_data = tmp.T

if args.save_npz:
    # Save data
    np.savez(args.save_npz, uniform_data=uniform_data, *x)
else:
    # Plot data
    fig, ax = plt.subplots()
    ndim = uniform_data.ndim

    if ndim == 0:
        ax.scatter(0, uniform_data)
    elif ndim == 1:
        ax.plot(x[0], uniform_data)
    elif ndim == 2:
        pos = ax.imshow(uniform_data.T, origin='lower',
                        extent=[x[0][0], x[0][-1], x[1][0], x[1][-1]],
                        vmin=args.vmin, vmax=args.vmax, cmap=args.cmap)
        if not args.hide_axes:
            fig.colorbar(pos, ax=ax)
    else:
        raise NotImplementedError(f'{ndim}D plotting not implemented')

    if args.hide_axes:
        plt.axis('off')
    if args.xlim:
        ax.set_xlim(*args.xlim)
    if args.ylim:
        ax.set_ylim(*args.ylim)

    if args.save_plot:
        print(f'Saving to {args.save_plot}')
        plt.savefig(args.save_plot, dpi=200, bbox_inches='tight')
    else:
        plt.show()
