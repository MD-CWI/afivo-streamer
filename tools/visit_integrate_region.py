#!/usr/bin/env python2

# This script prints the volume integral of a variable, possibly constrained to
# a region, as well as the total integration volume

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import sys


def get_argparser():
    # Get and parse the command line arguments
    pr = argparse.ArgumentParser(
        description='''Script to integrate a density over a region''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog='visit -nowin -cli -s visit_integrate_region.py')
    pr.add_argument('filename', type=str,
                    help='Input file (e.g. sim.silo)')
    pr.add_argument('variable', type=str,
                    help='Name of variable')
    pr.add_argument('-rmin', nargs='+', type=float, default=[-1e10, -1e10, -1e10],
                    help='Minimal coordinates')
    pr.add_argument('-rmax', nargs='+', type=float, default=[1e10, 1e10, 1e10],
                    help='Maximal coordinates')
    pr.add_argument('-cyl', action='store_true',
                    help='Take cylindrical geometry into account')
    return pr


if __name__ == '__main__':
    pr = get_argparser()

    args = pr.parse_args()

    if args.cyl:
        # Radius * 2 * pi * variable
        expr_var = "coord(mesh)[0] * 6.283185307179586 * " + args.variable
        expr_const = "coord(mesh)[0] * 6.283185307179586 * zonal_constant(mesh, 1.0)"
    else:
        expr_var = args.variable
        expr_const = "zonal_constant(1.0)"

    names = ['tmp_var', 'tmp_const']
    expressions = [expr_var, expr_const]

    for name, expr in zip(names, expressions):
        DefineScalarExpression(
            name,
            """if(and(and(
            and(ge(min_coord(mesh, "X"), {}), le(max_coord(mesh, "X"), {})),
            and(ge(min_coord(mesh, "Y"), {}), le(max_coord(mesh, "Y"), {}))),
            and(ge(min_coord(mesh, "Z"), {}), le(max_coord(mesh, "Z"), {}))),
            {}, 0.0)""".format(args.rmin[0], args.rmax[0],
                               args.rmin[1], args.rmax[1],
                               args.rmin[2], args.rmax[2],
                               expr))

    OpenDatabase(args.filename)

    results = []
    for name in names:
        AddPlot("Pseudocolor", name)
        DrawPlots()
        Query("Weighted Variable Sum")
        results.append(GetQueryOutputValue())

    print('volume_integral integration_volume')
    print(results[0], results[1])

    sys.exit()
