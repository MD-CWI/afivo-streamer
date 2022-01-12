#!/usr/bin/env python

import argparse
import sys

def get_args():
    # Get and parse the command line arguments
    pr = argparse.ArgumentParser(
        description='''Script to get the cross section of each segment''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='''Usage: visit -nowin -cli -s visit_get_cross.py silo_file cross_file''',
        prog='visit -nowin -cli -s visit_get_cross.py')
    pr.add_argument('silo_file', type=str,
                    help='silo file (e.g. sim.silo)')
    pr.add_argument('-origin', type=float, nargs=3, required=True,
                    help='Origin for the cross section')
    pr.add_argument('-normal', type=float, nargs=3, required=True,
                    help='Normal vector for the cross section')
    pr.add_argument('-max_width', type=float, default=4e-3,
                    help='Maximal window width to consider around origin')
    pr.add_argument('-factor', type=float, default=0.5,
                    help='Determine area at factor * maximum')
    pr.add_argument('-var', type=str, default='N2_B3',
                    help='Variable to consider')
    pr.add_argument('-save_fig', type=str,
                    help='Save output to this figure')
    return pr.parse_args()


if __name__ == '__main__':
    args = get_args()

    OpenDatabase(args.silo_file)
    AddPlot("Pseudocolor", args.var, 1, 1)

    PseudocolorAtts = PseudocolorAttributes()
    PseudocolorAtts.limitsMode = PseudocolorAtts.CurrentPlot
    SetPlotOptions(PseudocolorAtts)

    AddOperator("Slice", 1)
    SetActivePlots(0)
    SliceAtts = SliceAttributes()

    # Determine the slice using a vector and a point
    SliceAtts.originType = SliceAtts.Point
    SliceAtts.originPoint = tuple(args.origin)

    SliceAtts.normal = tuple(args.normal)
    SliceAtts.axisType = SliceAtts.Arbitrary

    SliceAtts.upAxis = (0, 0, 1)
    SliceAtts.project2d = 1
    SliceAtts.interactive = 1
    SetOperatorOptions(SliceAtts, 0, 1)

    AddOperator("Box", 1)
    BoxAtts = BoxAttributes()
    BoxAtts.amount = BoxAtts.Some  # Some, All
    BoxAtts.minx = -args.max_width
    BoxAtts.maxx = args.max_width
    BoxAtts.miny = -args.max_width
    BoxAtts.maxy = args.max_width
    BoxAtts.minz = 0
    BoxAtts.maxz = 1
    BoxAtts.inverse = 0
    SetOperatorOptions(BoxAtts, 1, 1)

    DrawPlots()

    Query("Max", use_actual_data=1)
    res = GetQueryOutputObject()
    half_maximum = args.factor * res['max']

    AddOperator("Isovolume", 1)
    IsovolumeAtts = IsovolumeAttributes()
    IsovolumeAtts.lbound = half_maximum
    IsovolumeAtts.ubound = 1e+37
    IsovolumeAtts.variable = "default"
    SetOperatorOptions(IsovolumeAtts, 2, 1)

    DrawPlots()

    Query("2D area")
    res = GetQueryOutputObject()
    print(res['Surface Area'])

    if args.save_fig:
        SaveWindowAtts = SaveWindowAttributes()
        SaveWindowAtts.fileName = args.save_fig
        SaveWindowAtts.format = SaveWindowAtts.PNG
        SetSaveWindowAttributes(SaveWindowAtts)
        SaveWindow()

    Close()

    sys.exit()
