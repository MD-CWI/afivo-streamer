#!/usr/bin/env python3

import h5py
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='TODO')
parser.add_argument('h5file', type=str, help='HDF5 dataset')
parser.add_argument('var', type=str, help='variable')
parser.add_argument('index', type=int, help='index')
args = parser.parse_args()

with h5py.File(args.h5file, 'r') as f:
    X = f[args.var]
    run_number = f['run_number']
    output_index = f['output_index']
    print(run_number[args.index], output_index[args.index])
    plt.imshow(X[args.index].T, origin='lower')
    plt.show()
