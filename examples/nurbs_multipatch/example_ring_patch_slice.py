##################################################################
######################## include.py   ############################
##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Institute for Structural Mechanics, RUB   #####
##### all rights reserved                                    #####
##################################################################
##################################################################
##################################################################
##################################################################
import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

bsplines_patch_util = BSplinesPatchUtility()

import geometry_factory

## create arc 1
r1 = 1.0
arc1_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r1, 0.0, 90.0) # put 'y' will give interesting ring
arc1 = arc1_ptr.GetReference()
arc1.Id = 1

## create arc 2
r2 = 2.0
arc2_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r2, 0.0, 90.0)
arc2 = arc2_ptr.GetReference()
arc2.Id = 2

## create ring patch by connect the two arcs
ring_patch_ptr = bsplines_patch_util.CreateLoftPatch(arc2, arc1)
ring_patch = ring_patch_ptr.GetReference()
ring_patch.Id = 1
print(ring_patch)

mpatch_export = MultiNURBSPatchMatlabExporter()
mpatch_export.Export(ring_patch, "ring.m")

spatch_ptr = ring_patch.ConstructSlicedPatch(1, 0.3)
spatch = spatch_ptr.GetReference()
spatch.Id = 2
print(spatch)
mpatch_export.Export(spatch, "arc.m")

spatch_ptr = ring_patch.ConstructSlicedPatch(0, 0.3)
spatch = spatch_ptr.GetReference()
spatch.Id = 3
print(spatch)
mpatch_export.Export(spatch, "line.m")
