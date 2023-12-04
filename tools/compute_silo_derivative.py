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
p.add_argument('input_files', type=str, nargs=2,
               help='Input Silo or raw files at times t1 and t2')
p.add_argument('raw_output', type=str,
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
p.add_argument('-silo_to_raw', type=str, default='./silo_to_raw',
               help='Path to silo_to_raw converter')


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

    grids_a, domain_a = load_file(args.input_files[0], args.project_dims,
                                  args.variable, args.silo_to_raw)
    grids_b, domain_b = load_file(args.input_files[1], args.project_dims,
                                  args.variable, args.silo_to_raw)
    grids_c = map_grid_a_to_b(grids_a, grids_b, args.axisymmetric,
                              args.interpolation)

    dt = domain_b['time'] - domain_a['time']
    inv_dt = 1/dt

    t_avg = 0.5 * (domain_b['time'] + domain_a['time'])
    cycle_avg = (domain_b['cycle'] + domain_a['cycle'])//2

    for gb, gc in zip(grids_b, grids_c):
        gc['values'] = (gb['values'] - gc['values']) * inv_dt

    domain_c = domain_b.copy()
    domain_c['time'] = t_avg
    domain_c['cycle'] = cycle_avg
    write_raw_data(args.raw_output, grids_c, domain_c)
    print(f'Time derivative has been written to {args.raw_output}')
