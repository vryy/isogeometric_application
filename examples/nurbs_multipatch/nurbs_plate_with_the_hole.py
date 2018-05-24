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
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
multipatch_refine_util = MultiPatchRefinementUtility()
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()

import geometry_factory

def CreateMultiPatch():
    ####### create arc 1
    r1 = 1.0
    arc1_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r1, 0.0, 45.0)
    arc1 = arc1_ptr.GetReference()

    # create line 1
    b1 = 4.0
    line1_ptr = geometry_factory.CreateLine([b1, 0.0, 0.0], [b1, b1, 0.0], arc1.Order(0))
    line1 = line1_ptr.GetReference()

    # create patch 1
    patch1_ptr = bsplines_patch_util.CreateConnectedPatch(line1, arc1)
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1

    ####### create arc 2
    arc2_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r1, 45.0, 90.0)
    arc2 = arc2_ptr.GetReference()

    # create line 2
    line2_ptr = geometry_factory.CreateLine([b1, b1, 0.0], [0.0, b1, 0.0], arc2.Order(0))
    line2 = line2_ptr.GetReference()

    # create patch 2
    patch2_ptr = bsplines_patch_util.CreateConnectedPatch(line2, arc2)
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2

    ######create multipatch
    mpatch = MultiPatch2D()
    mpatch.AddPatch(patch1_ptr)
    mpatch.AddPatch(patch2_ptr)
    bsplines_patch_util.MakeInterface(patch1, BoundarySide.Right, patch2, BoundarySide.Left, BoundaryDirection.Forward)
    multipatch_refine_util.DegreeElevate(patch1_ptr, [1, 1])
    multipatch_refine_util.InsertKnots(patch1_ptr, [[0.5], [0.5]])

    return mpatch

def main():
    mpatch = CreateMultiPatch()
    mpatch.Enumerate()
    print(mpatch)

    mpatch_export.Export(mpatch, "plate_with_hole.m")

if __name__ == "__main__":
    main()


