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
        description='''Generate Visit frames of 3D volume data''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pr.add_argument('database', type=str,
                    help='Database string (can have wildcard)')
    pr.add_argument('varname', type=str,
                    help='Plot this scalar quantity')

    pr.add_argument('-outdir', type=str, default=os.getcwd(),
                    help='Output directory')
    pr.add_argument('-fname', type=str, default='visit',
                    help='Output filename (without extension)')
    pr.add_argument('-fmt', type=str, default='png',
                    choices=['png', 'jpeg', 'bmp'],
                    help='Output format')
    pr.add_argument('-width', type=int, default='800',
                    help='Output width (pixels)')
    pr.add_argument('-height', type=int, default='800',
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
    pr.add_argument('-tstep', type=int, default=1,
                    help='Step between frames')

    pr.add_argument('-rndr', type=str, default='texture',
                    choices=['splatting', 'texture', 'raycast'],
                    help='Use this type of volume rendering')
    pr.add_argument('-scaling', type=str, default='lin',
                    choices=['lin', 'log'],
                    help='Scaling of data (lin/log)')
    pr.add_argument('-nsamples', type=int, default=15*1000*1000,
                    help='Number of samples for texture/splatting')
    pr.add_argument('-time', action='store_true',
                    help='Show time in output')
    pr.add_argument('-mesh', type=str,
                    help='Show mesh with this name')
    pr.add_argument('-triad', action='store_true',
                    help='Show triad')
    pr.add_argument('-legend', action='store_true',
                    help='Show legend in output')
    pr.add_argument('-dbase', action='store_true',
                    help='Show database in output')
    pr.add_argument('-bbox', type=float, nargs=6,
                    help='Restrict domain to (x0,x1, y0,y1, z0,z1)')
    pr.add_argument('-zoom', type=float, default=1.0,
                    help='Zoom factor')
    pr.add_argument('-pan', type=float, nargs=2, default=(0., 0.),
                    help='Image pan (relative x,y shift of window)')
    pr.add_argument('-normal', type=float, nargs=3, default=(0., -1., 0.),
                    help='View normal vector')

    pr.add_argument('-rdeg', type=float, default=6,
                    help='Rotation angle per time step (degrees)')
    pr.add_argument('-rfulldeg', type=float, default=360,
                    help='Full rotation angle (degrees)')
    pr.add_argument('-rfullpause', type=int, default=20,
                    help='Pause this number of frames before full rotation')
    pr.add_argument('-rsteps', type=int, default=2,
                    help='Number of rotation steps per time step')
    pr.add_argument('-rframes', type=int, nargs='+', default=[],
                    help='Perform full rotations at these time steps')
    pr.add_argument('-rend', action='store_true',
                    help='Perform full rotation at end')

    return pr.parse_args()


if __name__ == '__main__':
    args = get_args()

    # Treat databases as time-varying
    v.SetTreatAllDBsAsTimeVarying(1)

    v.OpenDatabase(args.database + ' database')
    v.AddPlot('Volume', args.varname, 1, 1)

    if args.mesh:
        v.AddPlot('Mesh', 'mesh')
        matts = v.MeshAttributes()
        v.SetPlotOptions(matts)
    v.DrawPlots()

    # Set rendering type
    vatts = v.VolumeAttributes()
    if args.rndr == 'splatting':
        vatts.rendererType = vatts.Splatting
    elif args.rndr == 'raycast':
        vatts.rendererType = vatts.RayCasting
    elif args.rndr == 'texture':
        vatts.rendererType = vatts.Texture3D

    # Set scaling
    if args.scaling == 'lin':
        vatts.scaling = vatts.Linear
    elif args.scaling == 'log':
        vatts.scaling = vatts.Log

    if args.ct:
        vatts.SetColorControlPoints(v.GetColorTable(args.ct))

    vatts.resampleFlag = 1
    vatts.resampleTarget = args.nsamples
    vatts.lightingFlag = 0

    if args.cmin is not None:
        vatts.useColorVarMin = 1
        vatts.colorVarMin = args.cmin
    if args.cmax is not None:
        vatts.useColorVarMax = 1
        vatts.colorVarMax = args.cmax

    vatts.opacityAttenuation = args.attenuation

    v.SetPlotOptions(vatts)

    # Set annotations
    aatts = v.AnnotationAttributes()
    aatts.axes3D.visible = 0
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

    # Set initial view
    cc = v.GetView3D()
    cc.viewNormal = tuple(args.normal)  # View x-plane
    cc.viewUp = (0, 0, 1)       # Z-axis points up
    cc.imageZoom = args.zoom
    cc.imagePan = tuple(args.pan)
    cc.perspective = 0
    v.SetView3D(cc)

    # Set box selection
    if args.bbox:
        v.AddOperator('Box')
        batts = v.BoxAttributes()
        batts.amount = batts.All  # Some, All
        batts.minx = args.bbox[0]
        batts.maxx = args.bbox[1]
        batts.miny = args.bbox[2]
        batts.maxy = args.bbox[3]
        batts.minz = args.bbox[4]
        batts.maxz = args.bbox[5]
        batts.inverse = 0
        v.SetOperatorOptions(batts)
        v.DrawPlots()

    # Set output options
    s = v.SaveWindowAttributes()
    if args.fmt == 'png':
        s.format = s.PNG
    elif args.r == 'bmp':
        s.format = s.BMP
    elif args.r == 'jpg':
        s.format = s.JPEG

    s.width = args.width
    s.height = args.height
    s.screenCapture = 0
    s.fileName = args.fname
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

    if args.rdeg == 0:
        args.rsteps = 0          # Don't do zero rotation

    if args.rend:
        args.rframes.append(imax-1)

    dphi = args.rdeg * (math.pi/180) / max(1, args.rsteps)
    # Same rotation speed at the end
    if dphi < 1.0e-3:
        full_steps = 0
        dphi_full = 0.0
    else:
        full_steps = int(round(args.rfulldeg * (math.pi/180) / dphi))
        dphi_full = args.rfulldeg * (math.pi/180) / full_steps

    for i in range(imin, imax, args.tstep):
        v.TimeSliderSetState(i)
        v.SaveWindow()
        for j in range(args.rsteps):
            cc.viewNormal = rotateXY(cc.viewNormal, dphi)
            v.SetView3D(cc)
            if (j > 0):
                v.SaveWindow()
        if i in args.rframes:
            for j in range(args.rfullpause):
                v.SaveWindow()
            for j in range(full_steps):
                cc.viewNormal = rotateXY(cc.viewNormal, dphi_full)
                v.SetView3D(cc)
                v.SaveWindow()

    sys.exit()
