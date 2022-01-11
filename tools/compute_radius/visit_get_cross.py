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
        description='''Script to get the cross section of each segment''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='''Usage: visit -nowin -cli -s visit_get_cross.py silo_file cross_file''',
        prog='visit -nowin -cli -s visit_get_cross.py ')
    pr.add_argument('silo_file', type=str,
                    help='silo file (e.g. sim.silo)')
    pr.add_argument('cross_file', type=str,default='./',
                    help='the output .txt file to get the point and vector of slice')
    pr.add_argument('-output_dir', type=str, default='./',
                    help='the output directory of vtk file')
    pr.add_argument('-vtk_file', type=str, default='vtk_file',
                    help='the output vtk file name')
    pr.add_argument('-single_CS', type=bool, default=False,
                    help='whether there is only one CS in the area')
    return pr


if __name__ == '__main__':
    pr = get_argparser()

    args = pr.parse_args()
   
    points_vectors = np.loadtxt(args.cross_file,skiprows=1)
    n_seg = points_vectors.shape[0]
    print(n_seg)
        
    
    OpenDatabase(args.silo_file)
    AddPlot("Pseudocolor", "N2_B3", 1, 1)
    
    for i in range(n_seg):
        
        print('iiiiiiiiiiiii',i)
        AddOperator("Slice", 1)
        SetActivePlots(0)
        SetActivePlots(0)
        SliceAtts = SliceAttributes()
        
        
        # determine the slice using a vector and a point
        SliceAtts.originType = SliceAtts.Point  
        
        # Point, Intercept, Percent, Zone, Node
        
        #point
        #SliceAtts.originPoint = (points_vectors[4],points_vectors[5],points_vectors[6])
        SliceAtts.originPoint = (points_vectors[i,1],points_vectors[i,2],points_vectors[i,3])
        SliceAtts.originIntercept = 0
        SliceAtts.originPercent = 50
        SliceAtts.originZone = 0
        SliceAtts.originNode = 0
        
        #arbitrary
        SliceAtts.normal = (points_vectors[i,4],points_vectors[i,5],points_vectors[i,6])
        SliceAtts.axisType = SliceAtts.Arbitrary 
        
        # XAxis, YAxis, ZAxis, Arbitrary, ThetaPhi
        SliceAtts.upAxis = (0, 0, 1)
        SliceAtts.project2d = 1
        SliceAtts.interactive = 1
        SliceAtts.flip = 0
        SliceAtts.originZoneDomain = 0
        SliceAtts.originNodeDomain = 0
        SliceAtts.meshName = "blk_mesh"
        SliceAtts.theta = 0
        SliceAtts.phi = 0
        SetOperatorOptions(SliceAtts, 0, 1)
        DrawPlots()
        
        #save the vtk file
        SaveWindowAtts = SaveWindowAttributes()
        SaveWindowAtts.outputToCurrentDirectory = 0
        SaveWindowAtts.outputDirectory = args.output_dir
        index = str(i)
        name= args.vtk_file+index
        SaveWindowAtts.fileName = name
        SaveWindowAtts.family = 0
        SaveWindowAtts.format = SaveWindowAtts.VTK  
        # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY, EXR
        SaveWindowAtts.width = 1024
        SaveWindowAtts.height = 1024
        SaveWindowAtts.screenCapture = 0
        SaveWindowAtts.saveTiled = 0
        SaveWindowAtts.quality = 80
        SaveWindowAtts.progressive = 0
        SaveWindowAtts.binary = 0
        SaveWindowAtts.stereo = 0
        #SaveWindowAtts.compression = SaveWindowAtts.None  
        # None, PackBits, Jpeg, Deflate, LZW
        SaveWindowAtts.forceMerge = 0
        SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions  
        # NoConstraint, EqualWidthHeight, ScreenProportions
        SaveWindowAtts.pixelData = 1
        SaveWindowAtts.advancedMultiWindowSave = 0
        SaveWindowAtts.subWindowAtts.win1.position = (0, 0)
        SaveWindowAtts.subWindowAtts.win1.size = (128, 128)
        SaveWindowAtts.subWindowAtts.win1.layer = 0
        SaveWindowAtts.subWindowAtts.win1.transparency = 0
        SaveWindowAtts.subWindowAtts.win1.omitWindow = 0
        SaveWindowAtts.subWindowAtts.win2.position = (0, 0)
        SaveWindowAtts.subWindowAtts.win2.size = (128, 128)
        SaveWindowAtts.subWindowAtts.win2.layer = 0
        SaveWindowAtts.subWindowAtts.win2.transparency = 0
        SaveWindowAtts.subWindowAtts.win2.omitWindow = 0
        SaveWindowAtts.subWindowAtts.win3.position = (0, 0)
        SaveWindowAtts.subWindowAtts.win3.size = (128, 128)
        SaveWindowAtts.subWindowAtts.win3.layer = 0
        SaveWindowAtts.subWindowAtts.win3.transparency = 0
        SaveWindowAtts.subWindowAtts.win3.omitWindow = 0
        SaveWindowAtts.subWindowAtts.win4.position = (0, 0)
        SaveWindowAtts.subWindowAtts.win4.size = (128, 128)
        SaveWindowAtts.subWindowAtts.win4.layer = 0
        SaveWindowAtts.subWindowAtts.win4.transparency = 0
        SaveWindowAtts.subWindowAtts.win4.omitWindow = 0
        SaveWindowAtts.subWindowAtts.win5.position = (0, 0)
        SaveWindowAtts.subWindowAtts.win5.size = (128, 128)
        SaveWindowAtts.subWindowAtts.win5.layer = 0
        SaveWindowAtts.subWindowAtts.win5.transparency = 0
        SaveWindowAtts.subWindowAtts.win5.omitWindow = 0
        SaveWindowAtts.subWindowAtts.win6.position = (0, 0)
        SaveWindowAtts.subWindowAtts.win6.size = (128, 128)
        SaveWindowAtts.subWindowAtts.win6.layer = 0
        SaveWindowAtts.subWindowAtts.win6.transparency = 0
        SaveWindowAtts.subWindowAtts.win6.omitWindow = 0
        SaveWindowAtts.subWindowAtts.win7.position = (0, 0)
        SaveWindowAtts.subWindowAtts.win7.size = (128, 128)
        SaveWindowAtts.subWindowAtts.win7.layer = 0
        SaveWindowAtts.subWindowAtts.win7.transparency = 0
        SaveWindowAtts.subWindowAtts.win7.omitWindow = 0
        SaveWindowAtts.subWindowAtts.win8.position = (0, 0)
        SaveWindowAtts.subWindowAtts.win8.size = (128, 128)
        SaveWindowAtts.subWindowAtts.win8.layer = 0
        SaveWindowAtts.subWindowAtts.win8.transparency = 0
        SaveWindowAtts.subWindowAtts.win8.omitWindow = 0
        SaveWindowAtts.subWindowAtts.win9.position = (0, 0)
        SaveWindowAtts.subWindowAtts.win9.size = (128, 128)
        SaveWindowAtts.subWindowAtts.win9.layer = 0
        SaveWindowAtts.subWindowAtts.win9.transparency = 0
        SaveWindowAtts.subWindowAtts.win9.omitWindow = 0
        SaveWindowAtts.subWindowAtts.win10.position = (0, 0)
        SaveWindowAtts.subWindowAtts.win10.size = (128, 128)
        SaveWindowAtts.subWindowAtts.win10.layer = 0
        SaveWindowAtts.subWindowAtts.win10.transparency = 0
        SaveWindowAtts.subWindowAtts.win10.omitWindow = 0
        SaveWindowAtts.subWindowAtts.win11.position = (0, 0)
        SaveWindowAtts.subWindowAtts.win11.size = (128, 128)
        SaveWindowAtts.subWindowAtts.win11.layer = 0
        SaveWindowAtts.subWindowAtts.win11.transparency = 0
        SaveWindowAtts.subWindowAtts.win11.omitWindow = 0
        SaveWindowAtts.subWindowAtts.win12.position = (0, 0)
        SaveWindowAtts.subWindowAtts.win12.size = (128, 128)
        SaveWindowAtts.subWindowAtts.win12.layer = 0
        SaveWindowAtts.subWindowAtts.win12.transparency = 0
        SaveWindowAtts.subWindowAtts.win12.omitWindow = 0
        SaveWindowAtts.subWindowAtts.win13.position = (0, 0)
        SaveWindowAtts.subWindowAtts.win13.size = (128, 128)
        SaveWindowAtts.subWindowAtts.win13.layer = 0
        SaveWindowAtts.subWindowAtts.win13.transparency = 0
        SaveWindowAtts.subWindowAtts.win13.omitWindow = 0
        SaveWindowAtts.subWindowAtts.win14.position = (0, 0)
        SaveWindowAtts.subWindowAtts.win14.size = (128, 128)
        SaveWindowAtts.subWindowAtts.win14.layer = 0
        SaveWindowAtts.subWindowAtts.win14.transparency = 0
        SaveWindowAtts.subWindowAtts.win14.omitWindow = 0
        SaveWindowAtts.subWindowAtts.win15.position = (0, 0)
        SaveWindowAtts.subWindowAtts.win15.size = (128, 128)
        SaveWindowAtts.subWindowAtts.win15.layer = 0
        SaveWindowAtts.subWindowAtts.win15.transparency = 0
        SaveWindowAtts.subWindowAtts.win15.omitWindow = 0
        SaveWindowAtts.subWindowAtts.win16.position = (0, 0)
        SaveWindowAtts.subWindowAtts.win16.size = (128, 128)
        SaveWindowAtts.subWindowAtts.win16.layer = 0
        SaveWindowAtts.subWindowAtts.win16.transparency = 0
        SaveWindowAtts.subWindowAtts.win16.omitWindow = 0
        SaveWindowAtts.opts.types = ()
        SaveWindowAtts.opts.help = ""
        SetSaveWindowAttributes(SaveWindowAtts)
        SaveWindow()
        RemoveOperator(0, 1)

    Close()

    sys.exit()
