#!/usr/bin/env python3

import pickle
import numpy as np
import argparse
from subprocess import run


def get_path_x(path, t):
    return path['x0'] + path['v'] * t + \
        0.5 * path['a'] * t**2


def get_path_v(path, t):
    return path['v'] + path['a'] * t


def get_args():
    pr = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pr.add_argument('silo_file', type=str,
                    help='silo file (e.g. sim.silo)')
    pr.add_argument('-paths', type=str, default='paths.pkl',
                    help='File with paths')
    pr.add_argument('-visit', type=str, default='visit',
                    help='Link to visit executable')
    return pr.parse_args()

args = get_args()

with open(args.paths, 'rb') as f:
    paths = pickle.load(f)

for i, path in enumerate(paths):
    # Define some location along the path
    t = 0.5 * (path['t0'] + path['t1'])
    origin = get_path_x(path, t)
    normal = get_path_v(path, t)

    # Use ' number' to allow for negative numbers as options
    visit_args = ['-nowin', '-cli', '-s', 'visit_get_cross.py'] + \
        [args.silo_file] + \
        ['-origin'] + [f' {str(x)}' for x in origin] + \
        ['-normal'] + [f' {str(x)}' for x in normal]

    # Add this to save figure output
    # + ['-save_fig', 'test.png']

    output = run([args.visit] + visit_args, capture_output=True, text=True)

    # First line of output contains area
    area = float(output.stdout.splitlines()[0])

    print(i, area)


