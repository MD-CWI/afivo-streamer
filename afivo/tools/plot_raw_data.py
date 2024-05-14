#!/usr/bin/env python3

# Tool to map Silo data to a uniform grid and plot it with python
#
# Author: Jannis Teunissen
#
# External requirements:
# pyabel (for Abel transform)

import numpy as np
import matplotlib.pyplot as plt
from raw_reader import load_file, map_grid_data_to
import argparse

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
p.add_argument('input_file', type=str, help='Input silo or raw file')
p.add_argument('-variable', type=str, help='Which variable to plot')
p.add_argument('-min_pixels', type=int, default=512,
               help='Min. pixels for any dimension')
p.add_argument('-project_dims', type=int, nargs='+', choices=[0, 1, 2],
               help='Project (integrate) along dimension(s)')
p.add_argument('-r_min', type=float, nargs='+',
               help='Only consider data above this coordinate')
p.add_argument('-r_max', type=float, nargs='+',
               help='Only consider data below this coordinate')
p.add_argument('-save_plot', type=str,
               help='Save the plot into this file')
p.add_argument('-axisymmetric', action='store_true',
               help='Assume axisymmetric geometry when averaging')
p.add_argument('-abel_transform', action='store_true',
               help='Perform forward Abel transform (implies axisymmetric)')
p.add_argument('-interpolation', default='linear',
               choices=['linear', 'nearest'],
               help='Interpolation method. Linear requires valid ghost cells'
               ' and is not conservative')
p.add_argument('-vmin', type=float, help='Minimum intensity in plot')
p.add_argument('-vmax', type=float, help='Maximum intensity in plot')
p.add_argument('-xlim', type=float, nargs=2, help='Plot range in x direction')
p.add_argument('-ylim', type=float, nargs=2, help='Plot range in y direction')
p.add_argument('-hide_axes', action='store_true',
               help='Hide axes etc.')
p.add_argument('-save_npz', type=str,
               help='Save the data into npz file and do not plot')
p.add_argument('-cmap', type=str, default='plasma',
               help='Use this colormap')
p.add_argument('-silo_to_raw', type=str, default=None,
               help='Path to silo_to_raw converter')
p.add_argument('-q', action='store_true',
               help='Print no information to stdout')


def get_uniform_data(grids, domain, min_pixels, interpolation='linear',
                     rmin=None, rmax=None, axisymmetric=False,
                     abel_transform=False):
    # Grid size nx should be of form 2^k * domain['n_cells_coarse']
    ratio = min_pixels / domain['n_cells_coarse'].min()
    ratio = 2**(np.ceil(np.log2(ratio)).astype(int))
    nx = domain['n_cells_coarse'] * ratio

    # Grid spacing on uniform grid
    dr = (domain['r_max'] - domain['r_min']) / nx

    r_min, r_max = domain['r_min'], domain['r_max']

    # Get coordinates of user-defined box on coarse grid
    if rmin is not None:
        r_min = domain['r_min'] + \
            np.floor((rmin-domain['r_min'])/dr) * dr
    if rmax is not None:
        r_max = domain['r_max'] - \
            np.floor((domain['r_max']-rmax)/dr) * dr

    nx = np.round((r_max - r_min)/dr).astype(int)

    if not args.q:
        print(f'Resolution:   {nx}')
        print(f'Grid spacing: {dr}')

    # List of coordinates
    x = [np.linspace(a, b, n) for a, b, n in
         zip(r_min+0.5*dr, r_max-0.5*dr, nx)]

    uniform_data = np.zeros(nx)

    # Map each grid to the uniformly spaced output
    for g in grids:
        axisymmetric = axisymmetric or abel_transform
        grid_data, ix_lo, ix_hi = map_grid_data_to(
            g, r_min, r_max, dr, axisymmetric, interpolation)

        g_ix = tuple([np.s_[i:j] for (i, j) in zip(ix_lo, ix_hi)])
        uniform_data[g_ix] += grid_data

    if abel_transform:
        from abel.hansenlaw import hansenlaw_transform
        if rmin is not None and rmin[0] > 0:
            raise ValueError('Abel transform requires r_min[0] to be 0.')
        tmp = hansenlaw_transform(uniform_data.T, dr[0], direction='forward')
        uniform_data = tmp.T

    return uniform_data, x


def plot_uniform_data(uniform_data, x, time, vmin=None, vmax=None, cmap=None,
                      xlim=None, ylim=None, hide_axes=False, save_plot=None):
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
                        vmin=vmin, vmax=vmax, cmap=cmap)
        if not hide_axes:
            fig.colorbar(pos, ax=ax)
    else:
        raise NotImplementedError(f'{ndim}D plotting not implemented')

    if hide_axes:
        plt.axis('off')
    if xlim:
        ax.set_xlim(*xlim)
    if ylim:
        ax.set_ylim(*ylim)

    ax.set_title(f't = {time:.3e}')

    if save_plot:
        print(f'Saving to {save_plot}')
        plt.savefig(save_plot, dpi=200, bbox_inches='tight')
    else:
        plt.show()


if __name__ == '__main__':
    args = p.parse_args()

    grids, domain = load_file(args.input_file, args.project_dims,
                              args.axisymmetric, args.variable,
                              args.silo_to_raw)

    # No longer axisymmetric if projected
    axisymmetric = args.axisymmetric
    if args.project_dims is not None and 0 in args.project_dims:
        axisymmetric = False

    if domain['n_dims'] > 0:
        values, coords = get_uniform_data(grids, domain, args.min_pixels,
                                          args.interpolation,
                                          args.r_min, args.r_max,
                                          axisymmetric,
                                          args.abel_transform)

        if args.save_npz:
            # Save data
            np.savez(args.save_npz, uniform_data=values,
                     cycle=domain['cycle'], time=domain['time'], *coords)
        else:
            plot_uniform_data(values, coords, domain['time'],
                              vmin=args.vmin, vmax=args.vmax,
                              cmap=args.cmap, xlim=args.xlim, ylim=args.ylim,
                              hide_axes=args.hide_axes,
                              save_plot=args.save_plot)
    else:
        # All spatial dimensions are projected, only print time and sum
        grid_values = np.array([g['values'] for g in grids])
        print(f'{domain["time"]:<16.8e} {grid_values.sum():<16.8e}')
