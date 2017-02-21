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
                    help='Database string (can have wildcard)')
    pr.add_argument('varname', type=str,
                    help='Plot this scalar quantity')

    return pr.parse_args()

def myQueryOverTime():
   for i in range(v.TimeSliderGetNStates()):
     v.SetTimeSliderState(i)
     v.Query("Time")
     t = v.GetQueryOutputValue()
     v.Query("Max")
     q = v.GetQueryOutputObject()
     t2 = q['max']
     r = q['max_coord']
     print(t, t2, r)

if __name__ == '__main__':
    args = get_args()

    v.LaunchNowin()

    v.OpenDatabase(args.database + ' database')
    n_states = v.GetDatabaseNStates()

    v.AddPlot("Pseudocolor", args.varname)
    v.DrawPlots()

    myQueryOverTime()



