#!/usr/bin/env python3

import argparse
import sys
import numpy as np
from pathlib import Path


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument('first_file', type=str,
                   help='First input file, e.g. path/sim_Emax_000001.txt')
    p.add_argument('-n', type=int, default=1000,
                   help='Search up to this index')
    p.add_argument('-max_points', type=int, default=10000,
                   help='Max total number of points')
    p.add_argument('-Emin', type=float, default=8e6,
                   help='Electric field threshold')
    p.add_argument('-dmax', type=float, default=1.5e-3,
                   help='Maximum distance between points on a path')
    p.add_argument('-min_points', type=int, default=5,
                   help='Minimum number of points on a path')
    p.add_argument('-lmin', type=float, default=1.5e-3,
                   help='Minimum length of a branch')
    return p.parse_args()


args = get_args()

base_name = args.first_file[:-10]
start_index = int(args.first_file[-10:-4])

points = np.zeros((args.max_points, 4))
times = np.zeros((args.max_points), dtype=int)
parent = np.full((args.max_points), -1, dtype=int)
n_children = np.zeros((args.max_points), dtype=int)
path_ix = np.full((args.max_points), 0, dtype=int)

t_i0 = np.zeros((args.n), dtype=int)
t_num = np.zeros((args.n), dtype=int)
n_points = 0

for i in range(start_index, args.n):
    fname = "{}{:06d}.txt".format(base_name, i)
    exists = Path(fname).is_file()
    if not exists:
        break
    end_index = i

    # Load data and convert single rows to 2D arrays
    d = np.genfromtxt(fname)
    if d.ndim == 1:
        d = d[np.newaxis, :]
    # Filter out points with a too low field
    d = d[d[:, 3] > args.Emin, :]

    # Store points
    t_i0[i] = t_i0[i-1] + t_num[i-1]
    t_num[i] = d.shape[0]
    points[t_i0[i]:t_i0[i]+t_num[i], :] = d
    times[t_i0[i]:t_i0[i]+t_num[i]] = i
    n_points += t_num[i]

if n_points == 0:
    print("No input files found")
    sys.exit(1)


def distance(a, b):
    return np.linalg.norm(a[0:3] - b[0:3])


# Set parents
for t in range(end_index, start_index, -1):
    for i in range(t_i0[t], t_i0[t] + t_num[t]):
        # Search for closest parent point
        dmin = 1e10
        for p in range(t_i0[t-1], t_i0[t-1] + t_num[t-1]):
            d = distance(points[i], points[p])
            if d < dmin:
                parent[i] = p
                dmin = d
        # Update number of children
        if dmin < args.dmax:
            n_children[parent[i]] += 1
        else:
            parent[i] = -1


def number_paths():
    i_path = -1
    for i in range(n_points):
        if path_ix[i] == -1:
            continue
        elif parent[i] == -1:
            i_path += 1
            path_ix[i] = i_path
        elif n_children[parent[i]] > 1:
            i_path += 1
            path_ix[i] = i_path
        else:
            path_ix[i] = path_ix[parent[i]]


number_paths()

while True:
    n_paths = path_ix.max() + 1
    print(n_paths)
    path_count = np.zeros((n_paths))
    path_length = np.zeros((n_paths))
    path_branches = np.zeros((n_paths), dtype=int)

    for i in range(n_paths):
        ixs = (path_ix[:n_points] == i).nonzero()[0]
        path_length[i] = np.linalg.norm(
            points[ixs[0], 0:3] - points[ixs[-1], 0:3])

    for i in range(n_points):
        if path_ix[i] != -1:
            path_count[path_ix[i]] += 1
            if n_children[i] > 1:
                path_branches[path_ix[i]] = 1

    # Remove short paths without branches
    remove_path = (path_count < args.min_points)
    remove_path |= (path_length < args.lmin)
    remove_path &= (path_branches == 0)
    # print(path_count)
    print(path_length)
    print(path_branches)
    print(remove_path)

    if np.sum(remove_path) == 0:
        break

    for i in range(n_points):
        if remove_path[path_ix[i]]:
            if parent[i] != -1:
                n_children[parent[i]] -= 1
            path_ix[i] = -1
    number_paths()

# sys.exit(1)
n_paths = path_ix.max()+1
paths = [{}]
for i in range(n_paths):
    pt = {}
    ixs = (path_ix[:n_points] == i).nonzero()[0]
    pt['points'] = points[ixs, :]
    pt['ix'] = i
    pt['t0'] = times[ixs[0]]
    pt['children'] = []
    if parent[ixs[0]] != -1:
        pt['parent'] = path_ix[parent[ixs[0]]]
        # paths[pt['parent']]['children'] += [i]
    else:
        pt['parent'] = -1
    paths.append(pt)

# for i in range(1, n_paths+1):
#     print(paths[i])

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect('equal')
for i in range(1, n_paths+1):
    pt = paths[i]
    ax.plot(pt['points'][:, 0], pt['points'][:, 1], pt['points'][:, 2],
            '-', label=str(pt['t0']))

ax.legend()
plt.show()
# for i in range(n_points):
#     print(i, times[i], parent[i], n_children[i], path_ix[i])
