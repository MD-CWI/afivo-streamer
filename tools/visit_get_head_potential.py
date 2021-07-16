#!/usr/bin/env python2

# Returns the potential at the location of max(E) (i.e., the streamer head) for
# two silo files

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import sys


def get_argparser():
    # Get and parse the command line arguments
    pr = argparse.ArgumentParser(
        description='''Script to get the head potential''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        prog='visit -nowin -cli -s visit_get_head_potential.py')
    pr.add_argument('filename_t0', type=str,
                    help='Silo file for background potential')
    pr.add_argument('filename_t1', type=str,
                    help='Silo file at later time')
    return pr


if __name__ == '__main__':
    pr = get_argparser()

    args = pr.parse_args()

    # Open file at t1
    OpenDatabase(args.filename_t1)
    AddPlot("Pseudocolor", "electric_fld")
    DrawPlots()

    # Get coordinate of max(E)
    Query("Max")
    coordinates = GetQueryOutputObject()['max_coord'] + (0.,)

    # Get potential
    AddPlot("Pseudocolor", "phi")
    DrawPlots()
    pick_t1 = ZonePick(coord=coordinates, vars=("phi"))

    ReplaceDatabase(args.filename_t0)
    pick_t0 = ZonePick(coord=coordinates, vars=("phi"))

    print("delta_phi phi_head phi_background location_head")
    print(pick_t1['phi'] - pick_t0['phi'], pick_t1['phi'], pick_t0['phi'], coordinates)

    Close()
    sys.exit()
