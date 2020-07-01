#!/usr/bin/env python2.7

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import os
import visit as v
import math
import sys


def rotateXY(xyz, phi):
    return (math.cos(phi) * xyz[0] - math.sin(phi) * xyz[1],
            math.sin(phi) * xyz[0] + math.cos(phi) * xyz[1], xyz[2])


def get_args():
    # Get and parse the command line arguments
    pr = argparse.ArgumentParser(
        description='''Generate Visit frames of 2D data.
        Usage visit -nowin -cli -s visit_plot_2d.py [options]''')
    pr.add_argument('database', type=str,
                    help='Database string (can have wildcard)')
    pr.add_argument('varname', type=str,
                    help='Plot this scalar quantity')

    pr.add_argument('-outdir', type=str, default=os.getcwd(),
                    help='Output directory')
    pr.add_argument('-width', type=int, default='2560',
                    help='Output width (pixels)')
    pr.add_argument('-height', type=int, default='2560',
                    help='Output height (pixels)')

    pr.add_argument('-ct', type=str,
                    help='Color table name')
    pr.add_argument('-cmin', type=float,
                    help='Minimum value')
    pr.add_argument('-cmax', type=float,
                    help='Maximum value')
    pr.add_argument('-attenuation', type=float, default=1.0,
                    help='Attenuatin (transparency)')
    pr.add_argument('-t0', type=int,
                    help='First frame (counting starts at 0)')
    pr.add_argument('-t1', type=int,
                    help='Last frame (counting starts at 0)')

    pr.add_argument('-contour', type=str,
                    help='Add contour of this variable')
    pr.add_argument('-contour_cmin', type=float,
                    help='Minimum value')
    pr.add_argument('-contour_cmax', type=float,
                    help='Maximum value')
    pr.add_argument('-contour_linewidth', type=int, default=1,
                    help='Line width')
    pr.add_argument('-contour_nlevels', type=int, default=10,
                    help='Number of levels')
    pr.add_argument('-contour_rgb', nargs=3, type=int, default=(192, 192, 192),
                    help='RGB color values for contour lines (0-255)')

    pr.add_argument('-scaling', type=str, default='lin',
                    choices=['lin', 'log'],
                    help='Scaling of data (lin/log)')
    pr.add_argument('-nodal', action='store_true',
                    help='Show nodal output')
    pr.add_argument('-time', action='store_true',
                    help='Show time in output')
    pr.add_argument('-axes', action='store_true',
                    help='Show axes in output')
    pr.add_argument('-mesh', type=str,
                    help='Show mesh with this name')
    pr.add_argument('-triad', action='store_true',
                    help='Show triad')
    pr.add_argument('-legend', action='store_true',
                    help='Show legend in output')
    pr.add_argument('-dbase', action='store_true',
                    help='Show database in output')
    pr.add_argument('-bbox', type=float, nargs=4,
                    help='Restrict domain to (x0,x1, y0,y1)')

    return pr.parse_args()

if __name__ == '__main__':
    args = get_args()

    # v.LaunchNowin()

    v.OpenDatabase(args.database + ' database')

    database_name = args.database.split(r'*')[0]
    database_name = database_name.split(r'/')[-1]

    v.AddPlot('Pseudocolor', args.varname, 1, 1)
    if args.mesh:
        v.AddPlot('Mesh', args.mesh)
        matts = v.MeshAttributes()
        matts.lineWidth = 1
        v.SetPlotOptions(matts)

    if args.contour:
        v.AddPlot("Contour", args.contour, 1, 1)
        catts = v.ContourAttributes()
        catts.contourNLevels = args.contour_nlevels
        if args.contour_cmin:
            catts.minFlag = 1
            catts.min = args.contour_cmin
        if args.contour_cmax:
            catts.maxFlag = 1
            catts.max = args.contour_cmax
        if args.contour_rgb:
            catts.colorType = catts.ColorBySingleColor
            catts.singleColor = (args.contour_rgb[0],
                                 args.contour_rgb[1],
                                 args.contour_rgb[2], 255)
        catts.lineStyle = catts.SOLID
        catts.lineWidth = args.contour_linewidth
        v.SetPlotOptions(catts)

    v.DrawPlots()

    vwa = v.View2DAttributes()
    if args.bbox:
        vwa.SetWindowCoords(args.bbox[0], args.bbox[1],
                            args.bbox[2], args.bbox[3])
        v.SetView2D(vwa)

    # Set rendering type
    vatts = v.PseudocolorAttributes()

    if args.nodal:
        vatts.centering = vatts.Nodal

    # Set scaling
    if args.scaling == 'lin':
        vatts.scaling = vatts.Linear
    elif args.scaling == 'log':
        vatts.scaling = vatts.Log

    if args.ct:
        vatts.colorTableName = args.ct

    if args.cmin is not None:
        vatts.minFlag = 1
        vatts.min = args.cmin
    if args.cmax is not None:
        vatts.maxFlag = 1
        vatts.max = args.cmax
    print(vatts.colorTableName)
    v.SetPlotOptions(vatts)

    # Set annotations
    aatts = v.AnnotationAttributes()
    if args.axes:
        aatts.axes2D.visible = 1
    else:
        aatts.axes2D.visible = 0

    if args.triad:
        aatts.axes3D.triadFlag = 1
    else:
        aatts.axes3D.triadFlag = 0
    aatts.axes3D.bboxFlag = 0
    aatts.userInfoFlag = 0
    aatts.timeInfoFlag = 0
    aatts.legendInfoFlag = 0
    aatts.databaseInfoFlag = 0
    if args.time:
        aatts.timeInfoFlag = 1
    if args.legend:
        aatts.legendInfoFlag = 1
    if args.dbase:
        aatts.databaseInfoFlag = 1

    v.SetAnnotationAttributes(aatts)

    # Set output options
    s = v.SaveWindowAttributes()
    s.format = s.PNG

    s.width = args.width
    s.height = args.height
    s.screenCapture = 0
    s.family = 0
    s.outputDirectory = args.outdir
    s.outputToCurrentDirectory = 0
    v.SetSaveWindowAttributes(s)

    if args.t1 is not None:
        imax = args.t1 + 1
    else:
        imax = v.TimeSliderGetNStates()

    if args.t0:
        imin = args.t0
    else:
        imin = 0

    for i in range(imin, imax):
        v.TimeSliderSetState(i)
        s.fileName = database_name + \
            args.varname + '_' + str(i)
        print(s.fileName)
        v.SetSaveWindowAttributes(s)
        v.SaveWindow()

    sys.exit()
