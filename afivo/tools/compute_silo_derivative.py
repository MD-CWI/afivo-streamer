#!/usr/bin/env python3

# Tool to compute time derivative of Silo data on a changing mesh. The data
# from time t1 is projected onto the mesh at time t2.
#
# Author: Jannis Teunissen

import numpy as np
from raw_reader import map_grid_data_to, write_raw_data
from plot_raw_data import load_file
import argparse
import copy

p = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
p.add_argument('input_files', type=str, nargs='+',
               help='Input Silo or raw files at times t_1, t_2, ..., t_n')
p.add_argument('-deriv_type', type=str, default='1st_central',
               choices=['1st_central', '2nd_central'],
               help='Type of derivative to compute')
p.add_argument('-output', type=str, required=True,
               help='Output raw file with the time derivative')
p.add_argument('-variable', type=str, help='Which variable to consider')
p.add_argument('-project_dims', type=int, nargs='+', choices=[0, 1, 2],
               help='Project (integrate) along dimension(s)')
p.add_argument('-axisymmetric', action='store_true',
               help='Assume axisymmetric geometry when averaging')
p.add_argument('-interpolation', default='nearest',
               choices=['linear', 'nearest'],
               help='Interpolation method. Linear requires valid ghost cells'
               ' and is not conservative')


def map_grid_a_to_b(grids_a, grids_b, axisymmetric=False,
                    interpolation='linear'):
    """Map grids_a to mesh of grids_b. All ghost cells will be set to zero.

    :param grids_a: will be projected onto grids_b
    :param grids_b: target grid
    :param axisymmetric: whether the data is axisymmetric
    :param interpolation: interpolation method

    """
    grids_c = copy.deepcopy(grids_b)

    for gc in grids_c:
        gc['values'][:] = 0.    # Erase all values

    for gc in grids_c:
        for ga in grids_a:
            # Determine which grids from grid_b overlap
            overlap = np.all(np.minimum(gc['r_max'], ga['r_max']) -
                             np.maximum(gc['r_min'], ga['r_min']) > 1e-10)

            if overlap:
                grid_data, ix_lo, ix_hi = map_grid_data_to(
                    ga, gc['r_min'], gc['r_max'], gc['dr'],
                    args.axisymmetric, args.interpolation)

                # Account for indexing offset in 'values'
                g_ix = tuple([np.s_[i:j] for (i, j) in
                              zip(ix_lo+gc['ilo'], ix_hi+gc['ilo'])])
                gc['values'][g_ix] = grid_data

    return grids_c


if __name__ == '__main__':
    args = p.parse_args()

    all_gd = [load_file(f, args.project_dims, args.variable)
              for f in args.input_files]

    grids = [gd[0] for gd in all_gd]
    domains = [gd[1] for gd in all_gd]
    times = np.array([d['time'] for d in domains])
    cycles = np.array([d['cycle'] for d in domains])

    if args.deriv_type == '1st_central':
        if len(grids) != 2:
            raise ValueError('1st_central requires len(grids) == 2')

        # Map onto last grid
        grids[0] = map_grid_a_to_b(grids[0], grids[1],
                                   args.axisymmetric, args.interpolation)

        # Compute 1st order derivative
        dt = times[1] - times[0]
        w = np.array([-1, 1]) / dt

        # Store results in grids[0]
        for i in range(len(grids[0])):
            grids[0][i]['values'] = w[0] * grids[0][i]['values'] + \
                w[1] * grids[1][i]['values']

        domains[0]['time'] = 0.5 * (times[0] + times[1])
        domains[0]['cycle'] = int((cycles[0] + cycles[1])//2)

    elif args.deriv_type == '2nd_central':
        if len(grids) != 3:
            raise ValueError('2nd_central requires len(grids) == 3')

        # Map onto last grid
        grids[0] = map_grid_a_to_b(grids[0], grids[2],
                                   args.axisymmetric, args.interpolation)
        grids[1] = map_grid_a_to_b(grids[1], grids[2],
                                   args.axisymmetric, args.interpolation)
        t_grid = times[1]

        # Compute 2nd order derivative, allowing for variable dt
        dt_a = times[1] - times[0]
        dt_b = times[2] - times[1]
        w = np.zeros(3)

        w[0] = dt_b/dt_a * 2 / (dt_b**2 + dt_a * dt_b)
        w[2] = 2 / (dt_b**2 + dt_a * dt_b)
        w[1] = -(w[0] + w[2])

        # Store results in grids[0]
        for i in range(len(grids[0])):
            grids[0][i]['values'] = w[0] * grids[0][i]['values'] + \
                w[1] * grids[1][i]['values'] + \
                w[2] * grids[2][i]['values']

        domains[0]['time'] = times[1]
        domains[0]['cycle'] = cycles[1]

    write_raw_data(args.output, grids[0], domains[0])
    print(f'Written {args.output}')
