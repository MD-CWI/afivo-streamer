#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import os
import sys
import numpy as np
import tempfile


def get_argparser():
    # Get and parse the command line arguments
    pr = argparse.ArgumentParser(
        description='''Save data of a variable over time using Visit in curve format.
        The output file contains two columns: the times and the
        variable at the point.''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='''Usage: visit -nowin -cli -s visit_plot_dataOverTime.py
        database variable -r pt''',
        prog='visit -nowin -cli -s visit_plot_dataOverTime.py')
    pr.add_argument('database', type=str,
                    help='Database name (e.g. sim.silo or "sim_*.silo")')
    pr.add_argument('varname', type=str,
                    help='Extract this (scalar) variable')
    #pr.add_argument('-r', type=float, nargs=3, required=True,
    #                metavar=('x', 'y', 'z'),
    #                help='Point to sample at(z=0 in 2D)')
    pr.add_argument('-r', type=float, required = True,
                    help='Point to be sampled')
    pr.add_argument('-output', type=str, default='point',
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
    timeSteps = v.TimeSliderGetNStates()

    LineoutAtts = v.LineoutAttributes()
    #Add input functionality for this later
    #LineoutAtts.point1 = tuple(args.r0)
    #LineoutAtts.point2 = tuple(args.r1)
    LineoutAtts.point1 = (0.0, 0.0, 0.0)
    LineoutAtts.point2 = (0.0, 16e-3, 0.0)
    LineoutAtts.interactive = 0
    LineoutAtts.ignoreGlobal = 0
    v.SetOperatorOptions(LineoutAtts, 1)

    v.DrawPlots()

    if not args.i0:
        args.i0 = 0
    if not args.i1:
        args.i1 = v.TimeSliderGetNStates()

    t_values = np.zeros(args.i1 - args.i0)
    var_values = np.zeros(args.i1 - args.i0)
    with tempfile.NamedTemporaryFile() as fname:
        it = 0
        for i in range(args.i0, args.i1):
            print(fname.name, i)
            # Set output options
            s = v.SaveWindowAttributes()
            s.format = s.CURVE
            s.outputDirectory = os.path.dirname(fname.name)
            s.fileName = os.path.basename(fname.name)
            s.outputToCurrentDirectory = 0
            s.family = 0
            v.SetSaveWindowAttributes(s)
    
            v.SetTimeSliderState(i)
            Query(str('Time'))
            t_values[it] = v.GetQueryOutputValue()
            v.SaveWindow()

            line_data = np.genfromtxt(fname.name + '.curve')
            var_values[it] = np.interp(args.r, line_data[:,0], line_data[:,1])
            it = it+1
    np.savetxt(args.varname+'_point_vs_time.txt', np.vstack([t_values, var_values]).T)
    print('Saved {}'.format(args.output))
    sys.exit()
