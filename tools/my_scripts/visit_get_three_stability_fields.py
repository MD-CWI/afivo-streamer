#!/usr/bin/env python2

# Returns the three types of stability field, delta_phi and streamer length

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
import sys
import numpy as np


def get_argparser():
    # Get and parse the command line arguments
    pr = argparse.ArgumentParser(
        description='''Script to get three types of stability field''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='''Usage: visit -nowin -cli -s visit_get_three_stability_fields.py
        filename -output''',
        prog='visit -nowin -cli -s visit_get_three_stability_fields.py')
    #pr.add_argument('filename_t0', type=str,
                    #help='Silo file at time=0')
    pr.add_argument('filename', type=str, nargs='+',
                    help='Silo file at t0 and t1')
    pr.add_argument('-output', type=str, default='stability_field.txt',
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
    v.OpenDatabase(args.filename[0])
    v.AddPlot("Pseudocolor", "electric_fld")
    v.DrawPlots()
    
    # Get coordinate of streamer tail(electrode tip)
    v.Query("Max")
    coordinates_tail = v.GetQueryOutputObject()['max_coord'] + (0.,)

    it_values = []
    filename_t1 = []    
    var_values = np.zeros((int(len(args.filename)/2), 15))
    j = 0
    
    print("it filename Est_dV Est_dV0 Est_U delta_phi delta_phi_0 phi_head pick_head_0 phi_tail L coord_head coord_tail") 
    for i, f in enumerate(args.filename):
        if (i % 2) == 0:
            # Update file at t0
            filename_t0 = f
        else:
            # Replace file at t1
            v.ReplaceDatabase(f)

            # Get coordinate of max(E)
            v.Query("Max")
            coordinates_head = v.GetQueryOutputObject()['max_coord'] + (0.,)

            # Get the potential at streamer tail and head
            pick_tail = v.ZonePick(coord=coordinates_tail, vars=("phi"))
            pick_head = v.ZonePick(coord=coordinates_head, vars=("phi"))
            
            # Replace file at t0
            v.ReplaceDatabase(filename_t0)
            # Get the original potential at streamer head 
            pick_head_0 = v.ZonePick(coord=coordinates_head, vars=("phi"))
        
            # Get delta_phi and delta_phi_0
            delta_phi = pick_head['phi'] - pick_tail['phi']
            delta_phi_0 = pick_head_0['phi'] - pick_tail['phi']

            # Get steamer length
            L = np.abs(coordinates_head[1] - coordinates_tail[1])
    
            # Get three types of stability field
            Est_dV = np.abs(delta_phi/L)
            Est_dV0 = np.abs(delta_phi_0/L)
            Est_U = np.abs(pick_tail['phi']/L)
        
            it_values += [j]
            filename_t1 += [f]
            var_values[j, 0] = Est_dV
            var_values[j, 1] = Est_dV0
            var_values[j, 2] = Est_U
            var_values[j, 3] = delta_phi
            var_values[j, 4] = delta_phi_0
            var_values[j, 5] = pick_head['phi']
            var_values[j, 6] = pick_head_0['phi']
            var_values[j, 7] = pick_tail['phi']
            var_values[j, 8] = L
            var_values[j, 9] = coordinates_head[0]
            var_values[j, 10] = coordinates_head[1]
            var_values[j, 11] = coordinates_head[2]
            var_values[j, 12] = coordinates_tail[0]
            var_values[j, 13] = coordinates_tail[1]
            var_values[j, 14] = coordinates_tail[2]

            print(j, f, Est_dV, Est_dV0, Est_U, delta_phi, delta_phi_0, pick_head['phi'],                               pick_head_0['phi'], pick_tail['phi'], L, coordinates_head, coordinates_tail)
            j += 1
         
    np.savetxt(args.output, np.column_stack([it_values, filename_t1, var_values]), header = 'it filename Est_dV Est_dV0 Est_U delta_phi delta_phi_0 phi_head phi_head_0 phi_tail L head_x head_y head_z tail_x tail_y tail_z', comments = '', delimiter = '  ', fmt='%s')
    #np.savetxt(args.output, np.column_stack([it_values, var_values]), header = 'it Est_dV Est_dV0 Est_U delta_phi delta_phi_0 phi_head pick_head_0 phi_tail L head_x head_y head_z tail_x tail_y tail_z', comments = '', delimiter = '', fmt = ['%d'] + ['%.8e']*15)
    print('Saved {}'.format(args.output))
    
    sys.exit()
