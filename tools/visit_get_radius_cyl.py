#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import os
import tempfile
import sys
import numpy as np


def get_argparser():
    # Get and parse the command line arguments
    pr = argparse.ArgumentParser(
        description='''Script to determine the streamer radius versus
        z-coordinate in axisymmetric simulation data''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='''Usage: visit -nowin -cli -s visit_get_radius_cyl.py
        silo_file -zrange z0 z1''',
        prog='visit -nowin -cli -s visit_get_radius_cyl.py')
    pr.add_argument('silo_file', type=str,
                    help='Silo file (e.g. sim.silo)')
    pr.add_argument('-zrange', nargs=2, type=float, required=True,
                    help='Minimal and maximal z coordinate')
    pr.add_argument('-npoints', type=int, default=50,
                    help='Number of z coordinates to use')
    pr.add_argument('-output', type=str, default='r_vs_z.txt',
                    help='File name for output')
    pr.add_argument('-threshold_e', type=float, default=5e17,
                    help='Threshold for electron density')
    return pr


if __name__ == '__main__':
    pr = get_argparser()

    try:
        import visit as v
    except ImportError:
        pr.print_help()
        raise

    args = pr.parse_args()

    v.OpenDatabase(args.silo_file)
    v.AddPlot("Curve", "operators/Lineout/" + "e", 1, 1)

    v.DrawPlots()

    z_values = np.linspace(args.zrange[0], args.zrange[1], args.npoints)
    r_values = np.zeros(args.npoints)

    with tempfile.NamedTemporaryFile() as fname:
        for i in range(args.npoints):
            print(fname.name, i)

            LineoutAtts = v.LineoutAttributes()
            LineoutAtts.point1 = (0., z_values[i], 0.)
            LineoutAtts.point2 = (1., z_values[i], 0.)
            LineoutAtts.interactive = 0
            LineoutAtts.ignoreGlobal = 0
            v.SetOperatorOptions(LineoutAtts, 1)

            # Set output options
            s = v.SaveWindowAttributes()
            s.format = s.CURVE
            s.outputDirectory = os.path.dirname(fname.name)
            s.fileName = os.path.basename(fname.name)
            s.outputToCurrentDirectory = 0
            s.family = 0
            v.SetSaveWindowAttributes(s)

            v.SaveWindow()

            r_data = np.genfromtxt(fname.name + '.curve')

            max_dens = r_data[:, 1].max()

            if args.threshold_e < max_dens:
                # Check at which index the density has dropped by factor 2
                j = np.where(r_data[:, 1] < args.threshold_e)[0][0]

                if j > 0:
                    # Linear interpolation to get intersection
                    w0 = r_data[j-1, 1]/args.threshold_e
                    w1 = r_data[j, 1]/args.threshold_e
                    c = (1 - w1)/(w0 - w1)
                    r_values[i] = c * r_data[j, 0] + (1-c) * r_data[j-1, 0]
                else:
                    r_values[i] = r_data[j, 0]
            else:
                r_values[i] = 0.

    np.savetxt(args.output, np.vstack([z_values, r_values]).T)
    print('Saved {}'.format(args.output))

    sys.exit()
