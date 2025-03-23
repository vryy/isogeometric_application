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

import geometry_factory

mpatch_export3 = MultiNURBSPatchMatlabExporter()
bsplines_patch_util = BSplinesPatchUtility()

def CreatePatch():
    arc1_ptr = geometry_factory.CreateHalfCircle([0.0, 0.0, 0.0], 'z', 1.0, start_angle=0.0)
    arc1 = arc1_ptr.GetReference()

    arc2_ptr = geometry_factory.CreateHalfCircle([0.0, 0.0, 0.0], 'z', 2.0, start_angle=0.0)
    arc2 = arc2_ptr.GetReference()

    # create patch
    patch_ptr = bsplines_patch_util.CreateLoftPatch(arc1, arc2)
    patch = patch_ptr.GetReference()
    patch.Id = 1

    # transpose the patch
    bsplines_patch_util.Transpose(patch_ptr)

    print(patch)

    return patch_ptr

def main():
    patch_ptr = CreatePatch()
    patch = patch_ptr.GetReference()
    patch.Id = 1
    patch.Prefix = "surface"
    mpatch_export3.Export(patch, "trans_surface.m")

if __name__ == "__main__":
    main()

