#!/usr/bin/env python2

# Returns the potential, background potential and their difference at the location of max(E) (i.e., the streamer head) and the coordinates of max(E) for a set of silo files

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import sys
import math
import numpy as np


def get_argparser():
    # Get and parse the command line arguments
    pr = argparse.ArgumentParser(
        description='''Script to get the potential, background potential and their difference at the location of max(E)''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='''Usage: visit -nowin -cli -s visit_get_delta_potential.py
        database -i0 -i1 -k -output''',
        prog='visit -nowin -cli -s visit_get_delta_potential.py')
    pr.add_argument('database', type=str,
                    help='Database name (e.g. "sim_*.silo"')
    pr.add_argument('-i0', type=int,
                    help='Start index in database (0 to N-1 for N files)')
    pr.add_argument('-i1', type=int,
                    help='Stop index in database (0 to N-1 for N files)')
    pr.add_argument('-k', type=int, default=1,
                    help='The number of timesteps for extracting a data in database')
    pr.add_argument('-output', type=str, default='delta_phi_output.txt',
                    help='File name for output')
    return pr


if __name__ == '__main__':
    pr = get_argparser()

    try:
        import visit as v
    except ImportError:
        pr.print_help()
        raise

    args = pr.parse_args()

    # Open the database
    use_database = ("*" in args.database)
    if use_database:
        v.OpenDatabase(args.database + " database", 0)
    else:
        v.OpenDatabase(args.database)
        
    if not args.i0:
        args.i0 = 0
    if not args.i1:
        args.i1 = v.TimeSliderGetNStates()
    
    npoint = int(math.ceil((args.i1-args.i0)/args.k))
    it_values = np.zeros(npoint, dtype=int)    
    time_values = np.zeros(npoint)
    var_values = np.zeros((npoint,6))
    
    print("it time delta_phi phi_current phi_background coord_head")
    j = 0
    for i in range(args.i0, args.i1, args.k):
        
        v.SetTimeSliderState(i)
        
        v.AddPlot("Pseudocolor", "electric_fld")
        v.DrawPlots()
        
        # Get time
        v.Query(str("Time"))
        it_values[j] = i
        time_values[j] = v.GetQueryOutputValue()
        
        # Get coordinate of max(E)
        v.Query("Max") 
        coordinates = v.GetQueryOutputObject()["max_coord"] + (0.,)
        
        # Get potential of streamer head at current time
        v.AddPlot("Pseudocolor", "phi")
        v.DrawPlots()
        pick_curr = v.ZonePick(coord=coordinates, vars=("phi"))
        
        # Get potential of streamer head at beginning time
        v.SetTimeSliderState(0)
        pick_back = v.ZonePick(coord=coordinates, vars=("phi")) 
        v.DeleteAllPlots()
        
        var_values[j, 0] = pick_curr['phi'] - pick_back['phi']
        var_values[j, 1] = pick_curr['phi']
        var_values[j, 2] = pick_back['phi']
        var_values[j, 3] = coordinates[0]
        var_values[j, 4] = coordinates[1]
        var_values[j, 5] = coordinates[2]
        
        print(it_values[j], time_values[j], pick_curr['phi'] - pick_back['phi'], pick_curr['phi'], pick_back['phi'], coordinates)
        j += 1        
    
    np.savetxt(args.output, np.column_stack([it_values, time_values, var_values]), header = 'it time delta_phi phi_current phi_background x y z', comments = '', delimiter = '', fmt = ['%d'] + ['%.8e']*7)
    print('Saved {}'.format(args.output))
    
    sys.exit()
