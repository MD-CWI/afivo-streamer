#!/usr/bin/env python3

import numpy as np
import h5py
import glob
import argparse
import re
from tqdm import tqdm
from itertools import compress

# Get and parse the command line arguments
parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Construct dataset from collection of npz files')
parser.add_argument('glob', type=str, help='Glob pattern for npz files')
parser.add_argument('h5file', type=str, help='HDF5 output file')
args = parser.parse_args()

# Read NPZ files
files = sorted(glob.glob(args.glob))

# For each file, extract run number, output index and variable name
# The files are assumed to be formed like ..._run1_000001_phi.npz
pattern = re.compile(r'.*_run([\d]+)_([\d]+)_(.*)\.npz$')

matches = [re.search(pattern, f) for f in files]
output_index = np.array([int(m.group(2)) for m in matches])
run_number = np.array([int(m.group(1)) for m in matches])
variable = [m.group(3) for m in matches]

# Remove files with output_index 0
mask = output_index != 0

N = sum(mask)
if N == 0:
    raise ValueError('No input files')

files = list(compress(files, mask))
variable = list(compress(variable, mask))
output_index = output_index[mask]
run_number = run_number[mask]

if len(set(variable)) != 1:
    raise ValueError('Multiple variables present, adjust glob')

# Get data for first file
tmp = np.load(files[0])
shape = tmp['uniform_data'].shape
r = tmp['arr_0']
z = tmp['arr_1']

# Append all the data to hdf5 file
with h5py.File(args.h5file, 'a') as f:
    if variable[0] not in f:
        X = f.create_dataset(variable[0], [N, *shape], dtype=float)
    else:
        X = f[variable[0]]

    if 'r' not in f:
        f.create_dataset('r', data=r)
        f.create_dataset('z', data=z)

    if 'run_number' not in f:
        f.create_dataset('run_number', data=run_number)
        f.create_dataset('output_index', data=output_index)

    for i in tqdm(range(N)):
        tmp = np.load(files[i])
        X[i] = tmp['uniform_data']
