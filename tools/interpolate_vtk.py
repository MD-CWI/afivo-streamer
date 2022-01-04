#!/usr/bin/env python

# This script can be used to interpolate an unstructured VTK file
#
# Author: Jannis Teunissen

import vtk
from vtk.util import numpy_support
import argparse
import numpy as np


def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("vtu", type=str,
                   help="VTU input file")
    p.add_argument("coords", type=str,
                   help="File with list of coordinates")
    p.add_argument("outfile", type=str,
                   help="File to store output in")
    p.add_argument("-var_name", type=str, default="phi",
                   help="Name of scalar variable in vtu file")
    p.add_argument("-xml", action='store_true',
                   help="Whether to use the XML reader")
    return p.parse_args()


args = get_args()

# Create a PLOT3D reader and load the data
if args.xml:
    reader = vtk.vtkXMLUnstructuredGridReader()
else:
    reader = vtk.vtkUnstructuredGridReader()

reader.SetFileName(args.vtu)
reader.Update()
reader_output = reader.GetOutput()

coords = np.loadtxt(args.coords)

if coords.ndim == 2:
    coords_3d = np.zeros([coords.shape[0], 3])
    coords_3d[:, 0:2] = coords
else:
    coords_3d = coords

points = vtk.vtkPoints()

for point in coords_3d:
    points.InsertNextPoint(point)

polydata = vtk.vtkPolyData()
polydata.SetPoints(points)

probe = vtk.vtkProbeFilter()
probe.SetInputData(polydata)
probe.SetSourceData(reader_output)
probe.Update()

result = probe.GetOutput()
pointdata = result.GetPointData()

phi_vtk = pointdata.GetArray(args.var_name)
phi = numpy_support.vtk_to_numpy(phi_vtk)

np.savetxt(args.outfile, phi)
