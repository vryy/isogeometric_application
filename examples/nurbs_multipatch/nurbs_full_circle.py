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

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
multipatch_refine_util = MultiPatchRefinementUtility()
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()
mpatch_export2 = MultiNURBSPatchGeoExporter()

import geometry_factory

def CreateMultiPatch(r):
    ####### create arc 1
    arc1_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r, -45.0, 45.0)
    arc1 = arc1_ptr.GetReference()

    # create line 1
    b = r/5
    line1_ptr = geometry_factory.CreateLine([b, -b, 0.0], [b, b, 0.0], arc1.Order(0))
    line1 = line1_ptr.GetReference()

    # create patch 1
    patch1_ptr = bsplines_patch_util.CreateLoftPatch(arc1, line1)
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1

    ####### create arc 2
    arc2_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r, 45.0, 135.0)
    arc2 = arc2_ptr.GetReference()

    # create line 2
    line2_ptr = geometry_factory.CreateLine([b, b, 0.0], [-b, b, 0.0], arc2.Order(0))
    line2 = line2_ptr.GetReference()

    # create patch 2
    patch2_ptr = bsplines_patch_util.CreateLoftPatch(arc2, line2)
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2

    ####### create arc 3
    arc3_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r, 135.0, 225.0)
    arc3 = arc3_ptr.GetReference()

    # create line 3
    line3_ptr = geometry_factory.CreateLine([-b, b, 0.0], [-b, -b, 0.0], arc3.Order(0))
    line3 = line3_ptr.GetReference()

    # create patch 3
    patch3_ptr = bsplines_patch_util.CreateLoftPatch(arc3, line3)
    patch3 = patch3_ptr.GetReference()
    patch3.Id = 3

    ####### create arc 4
    arc4_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r, 225.0, 315.0)
    arc4 = arc4_ptr.GetReference()

    # create line 4
    line4_ptr = geometry_factory.CreateLine([-b, -b, 0.0], [b, -b, 0.0], arc4.Order(0))
    line4 = line4_ptr.GetReference()

    # create patch 4
    patch4_ptr = bsplines_patch_util.CreateLoftPatch(arc4, line4)
    patch4 = patch4_ptr.GetReference()
    patch4.Id = 4

    ####### create line 5
    line5_ptr = geometry_factory.CreateLine([-b, -b, 0.0], [-b, b, 0.0], arc3.Order(0))
    line5 = line5_ptr.GetReference()

    # create patch 5
    patch5_ptr = bsplines_patch_util.CreateLoftPatch(line1, line5)
    multipatch_refine_util.DegreeElevate(patch5_ptr, [0, 1])
    patch5 = patch5_ptr.GetReference()
    patch5.Id = 5

    # # print(patch2)
    # print(patch3)
    # print("line2:")
    # print(line2)

    ######create multipatch
    mpatch = MultiPatch2D()
    mpatch.AddPatch(patch1_ptr)
    mpatch.AddPatch(patch2_ptr)
    mpatch.AddPatch(patch3_ptr)
    mpatch.AddPatch(patch4_ptr)
    mpatch.AddPatch(patch5_ptr)
    bsplines_patch_util.MakeInterface(patch1, BoundarySide.Right, patch2, BoundarySide.Left, BoundaryDirection.Forward)
    bsplines_patch_util.MakeInterface(patch2, BoundarySide.Right, patch3, BoundarySide.Left, BoundaryDirection.Forward)
    bsplines_patch_util.MakeInterface(patch3, BoundarySide.Right, patch4, BoundarySide.Left, BoundaryDirection.Forward)
    bsplines_patch_util.MakeInterface(patch4, BoundarySide.Right, patch1, BoundarySide.Left, BoundaryDirection.Forward)
    bsplines_patch_util.MakeInterface(patch1, BoundarySide.Top, patch5, BoundarySide.Bottom, BoundaryDirection.Forward)
    bsplines_patch_util.MakeInterface(patch2, BoundarySide.Top, patch5, BoundarySide.Right, BoundaryDirection.Forward)
    bsplines_patch_util.MakeInterface(patch3, BoundarySide.Top, patch5, BoundarySide.Top, BoundaryDirection.Reversed)
    bsplines_patch_util.MakeInterface(patch4, BoundarySide.Top, patch5, BoundarySide.Left, BoundaryDirection.Reversed)
    ## multipatch_refine_util.DegreeElevate(patch1_ptr, [1, 1])
    ## multipatch_refine_util.InsertKnots(patch1_ptr, [[0.5], [0.5]])

    return mpatch

def main():
    mpatch = CreateMultiPatch(0.025)
    mpatch.Enumerate()
    print(mpatch)

    mpatch_export.Export(mpatch, "full_circle.m")
    # mpatch_export2.Export(mpatch, "full_circle.geo")

if __name__ == "__main__":
    main()
