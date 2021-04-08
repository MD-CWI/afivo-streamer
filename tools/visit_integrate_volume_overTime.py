#!/usr/bin/env python2

# This script prints the volume integral of a variable, possibly constrained to
# a region, as well as the total integration volume

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import sys
import numpy as np
import tempfile

def get_argparser():
    # Get and parse the command line arguments
    pr = argparse.ArgumentParser(
        description='''Script to integrate a density over a region''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog='visit -nowin -cli -s visit_integrate_region.py')
    pr.add_argument('database', type=str,
                    help='Database name (e.g. sim.silo or "sim_*.silo"')
    pr.add_argument('variable', type=str,
                    help='Name of variable')
    pr.add_argument('-rmin', nargs='+', type=float, default=[-1e10, -1e10, -1e10],
                    help='Minimal coordinates')
    pr.add_argument('-rmax', nargs='+', type=float, default=[1e10, 1e10, 1e10],
                    help='Maximal coordinates')
    pr.add_argument('-i0', type=int,
                    help='Start index in database (0 to N-1 for N files)')
    pr.add_argument('-i1', type=int,
                    help='Stop index in database (0 to N-1 for N files)')
    pr.add_argument('-cyl', action='store_true',
                    help='Take cylindrical geometry into account')
    pr.add_argument('-output', type=str, default='point',
                    help='Base path name for output')
    return pr


if __name__ == '__main__':
    pr = get_argparser()

    try: 
        import visit as v
    except ImportError:
        pr.print_help()
        raise

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

    # Open either a database or a single file
    use_database = ('*' in args.database)
    if use_database:
        v.OpenDatabase(args.database + ' database', 0)
    else:
        v.OpenDatabase(args.database)
    if not args.i0:
        args.i0 = 0
    if not args.i1:
        args.i1 = v.TimeSliderGetNStates()

    t_values = np.zeros(args.i1-args.i0)
    var_values = np.zeros(args.i1 - args.i0)
    it = 0
    for i in range(args.i0, args.i1):
        results = []
        print("File index:", i)
        #Do the volume integration
        for name in names:
            v.AddPlot("Pseudocolor", name)
            v.DrawPlots()
            v.Query("Weighted Variable Sum")
            results.append(v.GetQueryOutputValue())
        print('Volume integral: ',results[0], 'Integrated volume: ', results[1])
        v.SetTimeSliderState(i)
        Query(str('Time'))
        t_values[it] = v.GetQueryOutputValue()
        var_values[it] = results[0]
        it = it+1

                

    np.savetxt(args.variable+'_integralOverVolume.txt',np.vstack([t_values, var_values]).T)
    print('Saved {}'.format(args.output))
    sys.exit()
