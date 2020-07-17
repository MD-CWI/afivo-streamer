#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import os
import sys


def get_argparser():
    # Get and parse the command line arguments
    pr = argparse.ArgumentParser(
        description='''Save lineout data using Visit in curve format.
        The output file contains two columns: the path length and the
        variable along the path.''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='''Usage: visit -nowin -cli -s visit_lineout.py
        database variable -r0 x y z -r1 x y z''',
        prog='visit -nowin -cli -s visit_lineout.py')
    pr.add_argument('database', type=str,
                    help='Database name (e.g. sim.silo or "sim_*.silo")')
    pr.add_argument('varname', type=str,
                    help='Extract this (scalar) variable')
    pr.add_argument('-r0', type=float, nargs=3, required=True,
                    metavar=('x0', 'y0', 'z0'),
                    help='Start location of line (z=0 in 2D)')
    pr.add_argument('-r1', type=float, nargs=3, required=True,
                    metavar=('x1', 'y1', 'z1'),
                    help='End location of line (z=0 in 2D)')
    pr.add_argument('-output', type=str, default='line',
                    help='Base path name for output')
    pr.add_argument('-i0', type=int,
                    help='Start index in database (0 to N-1 for N files)')
    pr.add_argument('-i1', type=int,
                    help='Stop index in database (0 to N-1 for N files)')
    return pr


if __name__ == '__main__':
    pr = get_argparser()

    try:
        import visit as v
    except ImportError:
        pr.print_help()
        raise

    args = pr.parse_args()

    # Open either a database or a single file
    use_database = ('*' in args.database)
    if use_database:
        v.OpenDatabase(args.database + ' database', 0)
    else:
        v.OpenDatabase(args.database)

    v.AddPlot("Curve", "operators/Lineout/" + args.varname, 1, 1)

    LineoutAtts = v.LineoutAttributes()
    LineoutAtts.point1 = tuple(args.r0)
    LineoutAtts.point2 = tuple(args.r1)
    LineoutAtts.interactive = 0
    LineoutAtts.ignoreGlobal = 0
    v.SetOperatorOptions(LineoutAtts, 1)

    v.DrawPlots()

    if not args.i0:
        args.i0 = 0
    if not args.i1:
        args.i1 = v.TimeSliderGetNStates()

    for i in range(args.i0, args.i1):
        # Set output options
        s = v.SaveWindowAttributes()
        s.format = s.CURVE
        s.outputDirectory = os.path.dirname(args.output)
        s.fileName = os.path.basename(args.output) + '_' + \
            args.varname + '_' + '{:04d}'.format(i)
        s.outputToCurrentDirectory = 0
        s.family = 0
        v.SetSaveWindowAttributes(s)

        v.TimeSliderSetState(i)
        v.SaveWindow()

    sys.exit()
