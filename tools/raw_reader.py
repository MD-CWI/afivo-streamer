#!/usr/bin/env python3

# Module to read and convert raw data obtained by transforming Silo files.
#
# Author: Jannis Teunissen (2022)

import numpy as np
import copy
from struct import unpack, calcsize
from scipy.interpolate import RegularGridInterpolator


def read_single_grid(f):
    """Read single grid from binary file"""
    n_dims = unpack('=i', f.read(calcsize('=i')))[0]

    fmt = '=' + str(n_dims) + 'i'
    dims = np.array(unpack(fmt, f.read(calcsize(fmt))))

    # lo and hi index range of non-phony data
    ilo = np.array(unpack(fmt, f.read(calcsize(fmt))))
    ihi = np.array(unpack(fmt, f.read(calcsize(fmt))))

    # Mesh coordinates
    coords = []
    coords_cc = []
    for d in dims:
        fmt = '=' + str(d) + 'd'
        tmp = np.array(unpack(fmt, f.read(calcsize(fmt))))
        coords.append(tmp)
        # Cell-centered coordinates
        coords_cc.append(0.5 * (tmp[1:] + tmp[:-1]))

    # Number of cell centers is one less than number of faces
    n_cells = dims-1
    fmt = '=' + str(np.product(n_cells)) + 'd'
    tmp = unpack(fmt, f.read(calcsize(fmt)))
    vals = np.array(tmp).reshape(n_cells, order='F')

    grid = {
        'n_dims': n_dims,
        'dims': dims,
        'ilo': ilo,
        'ihi': ihi,
        'coords': coords,
        'coords_cc': coords_cc,
        'values': vals
    }

    return grid


def get_raw_data(fname, project_dims=None):
    """Read raw data extracted from a Silo file

    :param fname: filename of raw data
    :param project_dims: int, integrate over these dimensions
    :returns: grids and domain properties

    """
    with open(fname, 'rb') as f:
        cycle = unpack('=i', f.read(calcsize('=i')))[0]
        time = unpack('=d', f.read(calcsize('=d')))[0]

        grids = []
        n_grids = unpack('=i', f.read(calcsize('=i')))[0]

        for i in range(n_grids):
            grid = read_single_grid(f)

            # Add properties
            props = get_grid_properties(grid)
            grid.update(props)

            if project_dims is not None:
                grid = grid_project(grid, project_dims)

            grids.append(grid)

    domain_props = get_domain_properties(grids)
    domain_props['cycle'] = cycle
    domain_props['time'] = time
    return grids, domain_props


def get_grid_properties(g):
    """Return some additional properties for a grid"""
    # Note: ghost cells are excluded for r_min and r_max
    r_min = [c[i] for c, i in zip(g['coords'], g['ilo'])]
    r_max = [c[i] for c, i in zip(g['coords'], g['ihi'])]
    dr = [c[1] - c[0] for c in g['coords']]

    props = {
        'r_min': np.array(r_min),
        'r_max': np.array(r_max),
        'dr': np.array(dr)
    }
    return props


def get_domain_properties(grids):
    """Return properties for a domain"""
    n_dims = grids[0]['n_dims']
    r_min = np.full(n_dims, 1e100)
    r_max = np.full(n_dims, -1e100)
    dr_max = np.zeros(n_dims)
    dr_min = np.full(n_dims, 1e100)

    for g in grids:
        r_min = np.minimum(r_min, g['r_min'])
        r_max = np.maximum(r_max, g['r_max'])
        dr_min = np.minimum(dr_min, g['dr'])
        dr_max = np.maximum(dr_max, g['dr'])

    nx_coarse = np.rint((r_max - r_min) / dr_max)

    props = {
        'n_dims': n_dims,
        'r_min': np.array(r_min),
        'r_max': np.array(r_max),
        'dr_min': np.array(dr_min),
        'dr_max': np.array(dr_max),
        'n_cells_coarse': nx_coarse.astype(int)
    }

    return props


def grid_project(in_grid, project_dims):
    """Integrate grid data along one or more dimensions"""
    g = copy.deepcopy(in_grid)

    pdims = np.sort(project_dims)

    # Index range to use for values below. Along the projected dimensions,
    # exclude ghost cells.
    ilo = 0 * g['dims']
    ihi = 1 * g['dims']
    ilo[pdims] = g['ilo'][pdims]
    ihi[pdims] = g['ihi'][pdims]

    # Get index range corresponding to non-ghost grid cells
    valid_ix = tuple([np.s_[i:j] for (i, j) in zip(ilo, ihi)])
    fac = np.product(g['dr'][pdims])
    g['values'] = g['values'][valid_ix].sum(axis=tuple(pdims)) * fac

    g['n_dims'] -= len(project_dims)
    g['dims'] = np.delete(g['dims'], pdims)
    g['ilo'] = np.delete(g['ilo'], pdims)
    g['ihi'] = np.delete(g['ihi'], pdims)
    g['r_min'] = np.delete(g['r_min'], pdims)
    g['r_max'] = np.delete(g['r_max'], pdims)
    g['dr'] = np.delete(g['dr'], pdims)

    for d in pdims[::-1]:       # Reverse order
        del g['coords'][d]
        del g['coords_cc'][d]

    return g


def map_grid_data_to(g, r_min, r_max, dr, axisymmetric=False,
                     interpolation_method='linear'):
    """Map grid data to another grid with origin r_min and grid spacing dr

    :param g: input grid
    :param r_min: origin of new grid
    :param r_max: upper location of new grid
    :param dr: grid spacing of new grid
    :param axisymmetric: whether to account for an axisymmetric geometry
    :param interpolation_method: how to interpolate data (linear, nearest)
    :returns: data and indices on new grid

    """
    # Check if grids overlap
    if np.any(g['r_min'] > r_max) or np.any(g['r_max'] < r_min):
        Ndim = len(r_min)
        # Return empty index range
        return 0., np.zeros(Ndim, dtype=int), np.zeros(Ndim, dtype=int)

    ratios = dr / g['dr']

    if not np.allclose(ratios, ratios[0], atol=0.):
        raise ValueError('Grid ratio not uniform')

    # Get index range corresponding to non-ghost grid cells
    valid_ix = tuple([np.s_[i:j] for (i, j) in zip(g['ilo'], g['ihi'])])

    # To avoid issues due to numerical round-off errors
    eps = 1e-10

    # Compute coarse grid min and max index
    if ratios[0] > 1 + eps:
        # Fine-to-coarse, first compute index for fine grid
        ix_lo_fine = np.round((g['r_min'] - r_min)/g['dr']).astype(int)
        ix_hi_fine = np.round((g['r_max'] - r_min)/g['dr']).astype(int) - 1

        # Then convert to coarse grid index
        r = np.round(ratios).astype(int)
        ix_lo, ix_hi = ix_lo_fine//r, ix_hi_fine//r
    else:
        # Grid is at same refinement level or coarser, so we can directly
        # compute index with rounding
        ix_lo = np.round((g['r_min'] - r_min)/dr).astype(int)
        ix_hi = np.round((g['r_max'] - r_min)/dr).astype(int) - 1

    nx = ix_hi - ix_lo + 1

    if ratios[0] > 1 + eps:
        # Reduce resolution. Determine coarse indices for each cell center
        cix = []
        coords_fine = []
        coords_coarse = []
        for d in range(g['n_dims']):
            dim_coords = g['coords_cc'][d]
            cc_fine = dim_coords[g['ilo'][d]:g['ihi'][d]]
            ix = np.floor((cc_fine - r_min[d])/dr[d]).astype(int)
            cc_coarse = r_min[d] + (ix + 0.5) * dr[d]

            cix.append(ix - ix_lo[d])
            coords_fine.append(cc_fine)
            coords_coarse.append(cc_coarse)

        # Create meshgrid of coarse indices
        ixs = np.meshgrid(*cix, indexing='ij')

        # Add fine grid values at coarse indices, weighted by relative volume
        cdata = np.zeros(nx)

        if axisymmetric:
            rvolume = np.product(g['dr']/dr) * coords_fine[0]/coords_coarse[0]
            # Broadcast to have volume weight for every grid cell
            rvolume = np.broadcast_to(rvolume[:, None], g['ihi']-g['ilo'])
            # Important to use Fortran order here
            values = rvolume.ravel() * g['values'][valid_ix].ravel()
        else:
            # Cartesian grid, simple averaging
            rvolume = np.product(g['dr']/dr)
            values = rvolume * g['values'][valid_ix].ravel()

        np.add.at(cdata, tuple(map(np.ravel, ixs)), values)
    elif ratios[0] < 1 - eps:
        # To interpolate data, compute coordinates of new cell centers
        # TODO: could maybe include axisymmetric correction here as well
        c_new = [np.linspace(a, b, n) for a, b, n in
                 zip(g['r_min']+0.5*dr, g['r_max']-0.5*dr, nx)]
        mgrid = np.meshgrid(*c_new, indexing='ij')
        new_coords = np.vstack(tuple(map(np.ravel, mgrid))).T

        # Near physical boundaries, extrapolation will be performed when using
        # linear interpolation
        f_interp = RegularGridInterpolator(
            tuple(g['coords_cc']), g['values'], bounds_error=False,
            fill_value=None, method=interpolation_method)

        cdata = f_interp(new_coords).reshape(nx)
    else:
        # Can directly use available data
        cdata = g['values'][valid_ix]

    # Get index range that is valid on uniform grid
    nx_uniform = np.round((r_max - r_min)/dr).astype(int)
    offset_lo = np.maximum(0, -ix_lo)
    offset_hi = np.maximum(0, ix_hi+1 - nx_uniform)

    # Index range on the grid_data
    grid_lo, grid_hi = offset_lo, nx-offset_hi
    d_ix = tuple([np.s_[i:j] for (i, j) in zip(grid_lo, grid_hi)])

    return cdata[d_ix], ix_lo+offset_lo, ix_hi+1-offset_hi
