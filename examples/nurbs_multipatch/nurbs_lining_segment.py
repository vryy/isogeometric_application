##################################################################
##################################################################
import sys
import os
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()

import geometry_factory

## create rings
rin = 4.8
rout = 5.0
leng = 15.0
sangle = 75.0
eangle = 105.0
ring1_ptr = geometry_factory.CreateSmallRing([0.0, 0.0, 0.0], 'x', rin, rout, sangle, eangle)
ring1 = ring1_ptr.GetReference()
ring2_ptr = geometry_factory.CreateSmallRing([leng, 0.0, 0.0], 'x', rin, rout, sangle, eangle)
ring2 = ring2_ptr.GetReference()

## create segment patch by connect the two rings
segment_patch_ptr = bsplines_patch_util.CreateLoftPatch(ring1, ring2)
segment_patch = segment_patch_ptr.GetReference()
print(segment_patch)
mpatch_export.Export(segment_patch, "segment.m")


