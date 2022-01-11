#!/usr/bin/env python

# Returns the potential at the location of max(E) (i.e., the streamer head) for
# two silo files

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import sys
import numpy as np


def get_argparser():
    # Get and parse the command line arguments
    pr = argparse.ArgumentParser(
        description='''Script to determine the streamer radius based on the cross section of each segment''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='''Usage: visit -nowin -cli -s visit_get_radius_cross.py vtk_file''',
        prog='visit -nowin -cli -s visit_get_radius_cross.py')
    pr.add_argument('vtk_file', type=str,
                    help='vtk file (e.g. sim.vtk)')

    return pr


if __name__ == '__main__':
    pr = get_argparser()

    args = pr.parse_args()

    # Open vtk file
    OpenDatabase(args.vtk_file,0)
    AddPlot("Pseudocolor", "N2_B3", 1, 1)
    DrawPlots()

    # Get the amplitude of N2_B3 density
    Query("Max")
    dens_max = GetQueryOutputObject()["max"]
    max_loc = GetQueryOutputObject()["max_coord"]+(0.,)
    print(max_loc)

    # get the effective area based on half amplitude
    AddOperator("Isovolume",1)
    SetActivePlots(0)
    SetActivePlots(0)
    IsovolumeAtts = IsovolumeAttributes()
    IsovolumeAtts.lbound = 0.5*dens_max
    IsovolumeAtts.ubound = dens_max
    IsovolumeAtts.variable = 'default'
    SetOperatorOptions(IsovolumeAtts, 0, 1)
    DrawPlots()

    # get the value of effective area
    Query("2D area")
    area = GetQueryOutputObject()["Surface Area"]
    radius = np.sqrt(area/np.pi)
    print("radius from", args.vtk_file, "=", radius)

    Close()
    sys.exit()
