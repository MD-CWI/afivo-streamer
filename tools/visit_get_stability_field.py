#!/usr/bin/env python2

# Returns the stability field, delta_phi and streamer length

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import sys
import numpy as np


def get_argparser():
    # Get and parse the command line arguments
    pr = argparse.ArgumentParser(
        description='''Script to get the stability field''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='''Usage: visit -nowin -cli -s visit_get_stability_field.py
        filename_t0 filename_t1 -output''',
        prog='visit -nowin -cli -s visit_get_stability_field.py')
    pr.add_argument('filename_t0', type=str,
                    help='Silo file at time=0')
    pr.add_argument('filename_t1', type=str, nargs='+',
                    help='Silo file at stagnating time')
    pr.add_argument('-output', type=str, default='stability_field_output.txt',
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

    # Open file at t0
    v.OpenDatabase(args.filename_t0)
    v.AddPlot("Pseudocolor", "electric_fld")
    v.DrawPlots()
    
    # Get coordinate of streamer tail(electrode tip)
    v.Query("Max")
    coordinates_tail = v.GetQueryOutputObject()['max_coord'] + (0.,)
    
    npoint = len(args.filename_t1)
    it_values = np.zeros(npoint, dtype=int)
    var_values = np.zeros((npoint, 11))
    
    print("it filename stability_field delta_phi phi_head phi_tail streamer_length coord_head coord_tail") 
    for i, f in enumerate(args.filename_t1):
        
        # Replace file at t1
        v.ReplaceDatabase(f)

        # Get coordinate of max(E)
        v.Query("Max")
        coordinates_head = v.GetQueryOutputObject()['max_coord'] + (0.,)

        # Get the potential at streamer tail and head
        pick_tail = v.ZonePick(coord=coordinates_tail, vars=("phi"))
        pick_head = v.ZonePick(coord=coordinates_head, vars=("phi"))
        
        # Get delta_phi
        delta_phi = pick_head['phi'] - pick_tail['phi']

        # Get steamer length
        L = np.abs(coordinates_head[1] - coordinates_tail[1])
    
        # Get stability field
        Est = np.abs(delta_phi/L)
        
        it_values[i] = i
        var_values[i, 0] = Est
        var_values[i, 1] = delta_phi
        var_values[i, 2] = pick_head['phi']
        var_values[i, 3] = pick_tail['phi']
        var_values[i, 4] = L
        var_values[i, 5] = coordinates_head[0]
        var_values[i, 6] = coordinates_head[1]
        var_values[i, 7] = coordinates_head[2]
        var_values[i, 8] = coordinates_tail[0]
        var_values[i, 9] = coordinates_tail[1]
        var_values[i, 10] = coordinates_tail[2]

        print(i, f, Est, delta_phi, pick_head['phi'], pick_tail['phi'], L, coordinates_head, coordinates_tail)
    
    np.savetxt(args.output, np.column_stack([it_values, args.filename_t1, var_values]), header = 'it filename Est delta_phi phi_head phi_tail L head_x head_y head_z tail_x tail_y tail_z', comments = '', delimiter = '  ', fmt='%s')
    #np.savetxt(args.output, np.column_stack([it_values, var_values]), header = 'it Est delta_phi phi_head phi_tail L head_x head_y head_z tail_x tail_y tail_z', comments = '', delimiter = '', fmt=''.join(['%3d'] + ['%17.8e']*11))
    print('Saved {}'.format(args.output))
    
    sys.exit()
