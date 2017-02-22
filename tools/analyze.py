#!/usr/bin/env python2

# This script requires that the visit python modules are in PYTHONPATH
# Use for example:
# VISITPYDIR="SOME_DIR/visit2_10_3.linux-x86_64/2.10.3/linux-x86_64/lib/site-packages"
# export PYTHONPATH=$PYTHONPATH:$VISITPYDIR

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import os
import visit as v

def get_args():
    # Get and parse the command line arguments
    pr = argparse.ArgumentParser(
        description='''Generate Visit frames of 3D volume data''')
    pr.add_argument('database', type=str,
                    help='Database string, e.g. "sim_*.silo" (with quotes)')

    return pr.parse_args()

def myQueryOverTime():
   for i in range(v.TimeSliderGetNStates()):
     v.SetTimeSliderState(i)

     # Activate electric field plot
     v.SetActivePlots(i_field)

     v.Query("Time")
     time = v.GetQueryOutputValue()

     v.Query("Max")
     q = v.GetQueryOutputObject()

     max_field = q['max']
     coord_max_field = q['max_coord']

     # Construct an xyz coordinate by appending zeros
     xyz_coord = coord_max_field + (3-len(coord_max_field)) * (0.0,)

     # As an example, pick electron and ion density at coordinate
     p = v.Pick(coord=xyz_coord, vars=("electron", "pos_ion"))

     print(time, max_field, coord_max_field, p['electron'], p['pos_ion'])

if __name__ == '__main__':
    args = get_args()

    v.LaunchNowin()

    v.OpenDatabase(args.database + ' database')

    # Add plots, and keep track of their index (can be used in SetActivePlots)
    i_field = 0
    v.AddPlot("Pseudocolor", "electric_fld")

    i_electron = 1
    v.AddPlot("Pseudocolor", "electron")

    i_pos_ion = 2
    v.AddPlot("Pseudocolor", "pos_ion")

    v.DrawPlots()

    myQueryOverTime()




