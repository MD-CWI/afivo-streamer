#!/usr/bin/env python3

import argparse
import sys
import numpy as np
from pathlib import Path
from sklearn.linear_model import LinearRegression, Lasso, HuberRegressor


def get_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('first_file', type=str,
                   help='First input file, e.g. path/sim_Emax_000001.txt')
    p.add_argument('-n', type=int, default=1000,
                   help='Search up to this index')
    p.add_argument('-dt', type=float, default=0.1e-9,
                   help='Time step for input files')
    p.add_argument('-max_points', type=int, default=10000,
                   help='Max total number of points')
    p.add_argument('-Emin', type=float, default=10e6,
                   help='Filter points below this threshold')
    p.add_argument('-Efac', type=float, default=0.8,
                   help='Minimum field compared to a parent point')
    p.add_argument('-dmax', type=float, default=0.3e-3,
                   help='Maximum distance between points on a path')
    p.add_argument('-min_points', type=int, default=10,
                   help='Minimum number of points on a path')
    p.add_argument('-lmin', type=float, default=0.5e-3,
                   help='Minimum length of a branch')
    p.add_argument('-branch_dt', type=float, default=1.0e-9,
                   help='Maximal difference in branch start time')
    p.add_argument('-branch_dmax', type=float, default=0.5e-3,
                   help='Maximal distance between branches')
    p.add_argument('-show_plot', action='store_true',
                   help='Show plot of the data')
    return p.parse_args()


# Distance between points
def distance(a, b):
    return np.linalg.norm(a[0:3] - b[0:3])


# Compute angle (in radians) between two vectors
def get_angle(va, vb):
    return np.arccos(np.dot(va, vb) / (np.linalg.norm(va) *
                                       np.linalg.norm(vb)))


# Set parents
def set_parents():
    n_children[:] = 0
    parent[:] = missing

    for t in range(end_index, start_index, -1):
        for i in range(t_i0[t], t_i0[t] + t_num[t]):
            # Search for closest parent point
            dmin = 1e10
            for p in range(t_i0[t-1], t_i0[t-1] + t_num[t-1]):
                if path_ix[p] == removed:
                    continue
                d = distance(points[i], points[p])
                if d < dmin:
                    dmin = d
                    # Field should not be much lower than parent
                    if points[i, 3] > args.Efac * points[p, 3]:
                        parent[i] = p

            # Update number of children
            if parent[i] >= 0:
                n_children[parent[i]] += 1


# Number the paths from 0 - N-1
def number_paths():
    i_path = -1
    for i in range(n_points):
        if path_ix[i] == removed:
            continue
        elif parent[i] == missing:
            i_path += 1
            path_ix[i] = i_path
        elif n_children[parent[i]] > 1:
            # Give branches a new index
            i_path += 1
            path_ix[i] = i_path
        elif distance(points[i], points[parent[i]]) > args.dmax:
            # Assume this is a new path
            i_path += 1
            path_ix[i] = i_path
        elif path_ix[parent[i]] < 0:
            i_path += 1
            path_ix[i] = i_path
        else:
            path_ix[i] = path_ix[parent[i]]


# Make sure paths are numbered from 0 - N-1
def renumber_paths():
    i_path = -1
    n_paths = path_ix.max() + 1
    new_ixs = np.full((n_paths), -1, dtype=int)

    for i in range(n_points):
        if path_ix[i] >= 0 and new_ixs[path_ix[i]] == -1:
            i_path += 1
            new_ixs[path_ix[i]] = i_path

    for i in range(n_points):
        if path_ix[i] >= 0:
            path_ix[i] = new_ixs[path_ix[i]]


def get_paths():
    n_paths = path_ix.max() + 1
    paths = [{'children': []} for i in range(n_paths)]

    for i in range(n_paths):
        ixs = (path_ix[:n_points] == i).nonzero()[0]
        paths[i]['ix'] = i
        paths[i]['n_points'] = ixs.size
        paths[i]['points'] = points[ixs, :]
        paths[i]['E0'] = points[ixs, 3].mean()
        paths[i]['t0'] = times[ixs[0]]
        paths[i]['t1'] = times[ixs[-1]]
        paths[i]['times'] = times[ixs]

        x0, v, a = fit_path(times[ixs], points[ixs, 0:3])
        paths[i]['x0'] = x0
        paths[i]['v'] = v
        paths[i]['a'] = a

        if parent[ixs[0]] != missing:
            paths[i]['parent'] = path_ix[parent[ixs[0]]]
            paths[paths[i]['parent']]['children'] += [i]
        else:
            paths[i]['parent'] = missing

    return paths


# Fit a curve x0 + v0 * t + 0.5 * a0 * t**2 through a path
def fit_path(t, xyz):
    # Doing this fit component-wise makes it easier to try other methods (for
    # example to deal with outliers)
    lr = LinearRegression()
    x0 = np.zeros((3))
    v = np.zeros((3))
    a = np.zeros((3))

    for dim in range(3):
        t_t2 = np.vstack([t, 0.5 * t**2]).T
        lr.fit(t_t2, xyz[:, dim])
        x0[dim] = lr.intercept_
        v[dim] = lr.coef_[0]
        a[dim] = lr.coef_[1]

    return x0, v, a


def get_path_x(path, t):
    return path['x0'] + path['v'] * t + 0.5 * path['a'] * t**2


def get_path_v(path, t):
    return path['v'] + path['a'] * t


def find_closest_point(xa, va, xb, vb):
    x = xb - xa
    v = vb - va
    t = -np.dot(x, v) / (v**2).sum()
    x1 = xa + t * va
    x2 = xb + t * vb
    return t, np.linalg.norm(x1-x2)


# Compute overlap between two intervals
def get_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


# Constants
missing = -1
removed = -2

# Get arguments
args = get_args()

# Name of input files
base_name = args.first_file[:-10]
start_index = int(args.first_file[-10:-4])

# Store points in arrays
points = np.zeros((args.max_points, 4))
times = np.zeros((args.max_points), dtype=int)
parent = np.full((args.max_points), missing, dtype=int)
n_children = np.zeros((args.max_points), dtype=int)
path_ix = np.full((args.max_points), 0, dtype=int)

t_i0 = np.zeros((args.n), dtype=int)
t_num = np.zeros((args.n), dtype=int)
n_points = 0

# Load all the data into the arrays
for i in range(start_index, args.n):
    fname = "{}{:06d}.txt".format(base_name, i)
    exists = Path(fname).is_file()
    if not exists:
        break
    end_index = i

    # Load data and convert single rows to 2D arrays
    d = np.genfromtxt(fname)

    if d.size == 0:
        # Create empty list with right size
        d = np.zeros((0, 4))
    elif d.ndim == 1:
        # Convert 1D array to 2D array
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

set_parents()
number_paths()

# Remove short paths
while True:
    n_paths = path_ix.max() + 1
    path_count = np.zeros((n_paths), dtype=int)
    path_length = np.zeros((n_paths))
    path_children = np.zeros((n_paths), dtype=int)

    for i in range(n_paths):
        ixs = (path_ix[:n_points] == i).nonzero()[0]
        path_length[i] = np.linalg.norm(
            points[ixs[0], 0:3] - points[ixs[-1], 0:3])

    for i in range(n_points):
        if path_ix[i] >= 0:
            path_count[path_ix[i]] += 1
            path_children[path_ix[i]] += n_children[i]

    # Remove short paths without branches
    remove_path = (path_count < args.min_points)
    remove_path |= (path_length < args.lmin)
    remove_path &= (path_children < path_count)

    if np.sum(remove_path) == 0:
        break

    for i in range(n_points):
        ip = path_ix[i]
        if ip != removed and remove_path[ip]:
            path_ix[i] = removed
            if parent[i] != missing:
                n_children[parent[i]] -= 1

    number_paths()

paths = get_paths()

# Merge short paths with adjacent ones
for path in paths:
    if path['n_points'] < args.min_points:
        tmean = 0.5 * (path['t0'] + path['t1'])
        xmean = get_path_x(path, tmean)

        nearby_paths = [path['parent']] + path['children']
        distances = np.zeros((len(nearby_paths)))
        for n, j in enumerate(nearby_paths):
            x = get_path_x(paths[j], tmean)
            distances[n] = np.linalg.norm(xmean-x)

        if distances.min() < args.dmax:
            i_closest = np.argmin(distances)
            path_ix[path_ix == path['ix']] = nearby_paths[i_closest]

# Remove remaining short paths
for path in paths:
    if path['n_points'] < args.min_points:
        path_ix[path_ix == path['ix']] = removed

renumber_paths()
paths = get_paths()
# TODO: compute temporal overlap (in symmetric way) between branches?

branchings = []

# Loop over pairs of path i, j
for j in range(len(paths)):
    # Find neighboring path
    for i in range(j):
        # The neighbor should start at around the same time
        if abs(paths[i]['t0'] - paths[j]['t0']) * args.dt > args.branch_dt:
            continue

        n_min = min(paths[i]['n_points'], paths[j]['n_points'])

        # The paths should be sufficiently long
        # if n_min < 2 * args.min_points:
        #     continue

        # The paths should have sufficient temporal overlap
        overlap = get_overlap([paths[i]['t0'], paths[i]['t1']],
                              [paths[j]['t0'], paths[j]['t1']])

        if overlap / n_min < 0.75:
            continue

        # Get position and velocity at average start time
        t0 = 0.5 * (paths[i]['t0'] + paths[j]['t0'])
        xa = get_path_x(paths[i], t0)
        va = get_path_v(paths[i], t0)
        xb = get_path_x(paths[j], t0)
        vb = get_path_v(paths[j], t0)

        t, d = find_closest_point(xa, va, xb, vb)
        t += t0

        if d < args.branch_dmax:
            x1 = get_path_x(paths[i], t)
            x2 = get_path_x(paths[j], t)
            xmean = 0.5 * (x1 + x2)

            # Merge with existing branching event?
            merged = False
            for b in branchings:
                xb = sum(b['x']) / len(b['x'])
                if np.linalg.norm(xb - xmean) < args.branch_dmax and \
                   abs(t - b['t']) * args.dt < args.branch_dt:
                    if i not in b['ixs']:
                        b['ixs'] += [i]
                        b['v'] += [va]
                        b['x'] += [x1]
                    if j not in b['ixs']:
                        b['ixs'] += [j]
                        b['v'] += [vb]
                        b['x'] += [x2]
                    b['t'] = (b['n'] * b['t'] + t) / (b['n'] + 1)
                    b['n'] += 1
                    merged = True
                    break

            if not merged:
                b = {'t': t,
                     'n': 1,
                     'ixs': [i, j],
                     'v': [va, vb],
                     'x': [x1, x2]}
                branchings.append(b)

for b in branchings:
    print(b['ixs'], b['x'], b['v'], b['t'])

if args.show_plot:
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    for pt in paths:
        ax.plot(pt['points'][:, 0], pt['points'][:, 1], pt['points'][:, 2],
                '.', label='{},{},{},{}'.format(pt['ix'], pt['parent'],
                                                pt['t0'], pt['points'].shape[0]))
        line = pt['x0'] + np.outer(pt['times'], pt['v']) + \
            0.5 * np.outer(pt['times']**2, pt['a'])
        ax.plot(line[:, 0], line[:, 1], line[:, 2], '-')

    # This requires a recent version of matplotlib
    ax.set_box_aspect([ub - lb for lb, ub in
                       (getattr(ax, f'get_{a}lim')() for a in 'xyz')])

    ax.legend()
    plt.show()
