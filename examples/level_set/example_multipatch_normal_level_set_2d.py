##################################################################
##################################################################
import sys
import os
import math
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.IsogeometricBRepApplication import *
kernel = Kernel()   #defining kernel

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export1 = MultiNURBSPatchGeoExporter()
mpatch_export2 = MultiNURBSPatchGLVisExporter()
mpatch_export3 = MultiNURBSPatchMatlabExporter()

import geometry_factory

def CreateMultiPatch(R = 10.5, L = 50.0, scale_value = 1.0):
    #############################################
    ############### patch 1 #####################
    ####### create arc 1
    arc1_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'y', R, 180.0, 90.0)
    arc1 = arc1_ptr.GetReference()

    ####### create arc 2
    arc2_ptr = geometry_factory.CreateSmallArc([0.0, L/2.0, 0.0], 'y', R*scale_value, 180.0, 90.0)
    arc2 = arc2_ptr.GetReference()

    # create patch 1
    patch1_ptr = bsplines_patch_util.CreateLoftPatch(arc1, arc2)
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1

    ############################################
    ############### patch 2 ####################
    ####### create arc 2_1
    arc3_ptr = geometry_factory.CreateSmallArc([0.0, L/2.0, 0.0], 'y', R*scale_value, 180.0, 90.0)
    arc3 = arc3_ptr.GetReference()

    ####### create arc 2_2
    arc4_ptr = geometry_factory.CreateSmallArc([0.0, L , 0.0], 'y', R, 180.0, 90.0)
    arc4 = arc4_ptr.GetReference()

    # create patch 2
    patch2_ptr = bsplines_patch_util.CreateLoftPatch(arc3, arc4)
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2

    ############################################
    ############### patch 3 ####################
    ####### create arc 3_1
    arc5_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'y', R, 90.0, 0.0)
    arc5 = arc5_ptr.GetReference()

    ####### create arc 3_2
    arc6_ptr = geometry_factory.CreateSmallArc([0.0, L/2.0, 0.0], 'y', R*scale_value, 90.0, 0.0)
    arc6 = arc6_ptr.GetReference()

    # create patch 2
    patch3_ptr = bsplines_patch_util.CreateLoftPatch(arc5, arc6)
    patch3 = patch3_ptr.GetReference()
    patch3.Id = 3

    ############################################
    ############### patch 4 ####################
    ####### create arc 4_1
    arc7_ptr = geometry_factory.CreateSmallArc([0.0, L/2.0, 0.0], 'y', R*scale_value, 90.0, 0.0)
    arc7 = arc7_ptr.GetReference()

    ####### create arc 4_2
    arc8_ptr = geometry_factory.CreateSmallArc([0.0, L, 0.0], 'y', R, 90.0, 0.0)
    arc8 = arc8_ptr.GetReference()

    # create patch 2
    patch4_ptr = bsplines_patch_util.CreateLoftPatch(arc7, arc8)
    patch4 = patch4_ptr.GetReference()
    patch4.Id = 4

    #############################################
    ############### patch 5 #####################
    ####### create arc 1

    arc9_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'y', R, 0.0, -90.0)
    arc9 = arc9_ptr.GetReference()

    ####### create arc 2
    arc10_ptr = geometry_factory.CreateSmallArc([0.0, L/2.0 , 0.0], 'y', R*scale_value, 0.0, -90.0)
    arc10 = arc10_ptr.GetReference()

    # create patch 1
    patch5_ptr = bsplines_patch_util.CreateLoftPatch(arc9, arc10)
    patch5 = patch5_ptr.GetReference()
    patch5.Id = 5

    ############################################
    ############### patch 6 ####################
    ####### create arc 2_1
    arc11_ptr = geometry_factory.CreateSmallArc([0.0, L/2.0, 0.0], 'y', R*scale_value, 0.0, -90.0)
    arc11 = arc11_ptr.GetReference()

    ####### create arc 2_2
    arc12_ptr = geometry_factory.CreateSmallArc([0.0, L, 0.0], 'y', R, 0.0, -90.0)
    arc12 = arc12_ptr.GetReference()

    # create patch 2
    patch6_ptr = bsplines_patch_util.CreateLoftPatch(arc11, arc12)
    patch6 = patch6_ptr.GetReference()
    patch6.Id = 6

    ############################################
    ############### patch 7 ####################
    ####### create arc 3_1
    arc13_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'y', R, -90.0, -180.0)
    arc13 = arc13_ptr.GetReference()

    ####### create arc 3_2
    arc14_ptr = geometry_factory.CreateSmallArc([0.0, L/2.0 , 0.0], 'y', R*scale_value, -90.0, -180.0)
    arc14 = arc14_ptr.GetReference()

    # create patch 2
    patch7_ptr = bsplines_patch_util.CreateLoftPatch(arc13, arc14)
    patch7 = patch7_ptr.GetReference()
    patch7.Id = 7

    ############################################
    ############### patch 8 ####################
    ####### create arc 4_1
    arc15_ptr = geometry_factory.CreateSmallArc([0.0, L/2.0, 0.0], 'y', R*scale_value, -90.0, -180.0)
    arc15 = arc15_ptr.GetReference()

    ####### create arc 4_2
    arc16_ptr = geometry_factory.CreateSmallArc([0.0, L, 0.0], 'y', R, -90.0, -180.0)
    arc16 = arc16_ptr.GetReference()

    # create patch 2
    patch8_ptr = bsplines_patch_util.CreateLoftPatch(arc15, arc16)
    patch8 = patch8_ptr.GetReference()
    patch8.Id = 8

    ######create multipatch
    mpatch = MultiPatch2D()
    mpatch.AddPatch(patch1_ptr)
    mpatch.AddPatch(patch2_ptr)
    mpatch.AddPatch(patch3_ptr)
    mpatch.AddPatch(patch4_ptr)
    mpatch.AddPatch(patch5_ptr)
    mpatch.AddPatch(patch6_ptr)
    mpatch.AddPatch(patch7_ptr)
    mpatch.AddPatch(patch8_ptr)

    bsplines_patch_util.MakeInterface(patch1, BoundarySide.Top, patch2, BoundarySide.Bottom, BoundaryDirection.Forward)
    bsplines_patch_util.MakeInterface(patch3, BoundarySide.Top, patch4, BoundarySide.Bottom, BoundaryDirection.Forward)
    bsplines_patch_util.MakeInterface(patch5, BoundarySide.Top, patch6, BoundarySide.Bottom, BoundaryDirection.Forward)
    bsplines_patch_util.MakeInterface(patch7, BoundarySide.Top, patch8, BoundarySide.Bottom, BoundaryDirection.Forward)

    bsplines_patch_util.MakeInterface(patch7, BoundarySide.Right, patch1, BoundarySide.Left, BoundaryDirection.Forward)
    bsplines_patch_util.MakeInterface(patch8, BoundarySide.Right, patch2, BoundarySide.Left, BoundaryDirection.Forward)
    bsplines_patch_util.MakeInterface(patch1, BoundarySide.Right, patch3, BoundarySide.Left, BoundaryDirection.Forward)
    bsplines_patch_util.MakeInterface(patch2, BoundarySide.Right, patch4, BoundarySide.Left, BoundaryDirection.Forward)

    bsplines_patch_util.MakeInterface(patch3, BoundarySide.Right, patch5, BoundarySide.Left, BoundaryDirection.Forward)
    bsplines_patch_util.MakeInterface(patch4, BoundarySide.Right, patch6, BoundarySide.Left, BoundaryDirection.Forward)
    bsplines_patch_util.MakeInterface(patch5, BoundarySide.Right, patch7, BoundarySide.Left, BoundaryDirection.Forward)
    bsplines_patch_util.MakeInterface(patch6, BoundarySide.Right, patch8, BoundarySide.Left, BoundaryDirection.Forward)

    return mpatch

######################################################################

def main():
    R = 1.5
    mpatch = CreateMultiPatch(R = R, L = 1.0, scale_value = 1.0)
    print("Multipatch is created")

    # mpatch_export2.Export(mpatch, "mpatch.mesh")
    mpatch_export3.Export(mpatch, "mpatch.m")

    system_size = mpatch.Enumerate()
    print("system_size:", system_size)
    # print(mpatch)

    multipatch_util.CheckInterfaces(mpatch)

    ls = MultiPatchNormalLevelSet2D(mpatch)
    # ls.SetEchoLevel(5)
    # # print(ls)

    P = Point3D(0.25, 0.25, 2.1)
    l = math.sqrt(0.25**2 + 2.1**2)
    lsv = ls.GetValue(P)
    lsv_ref = math.sqrt((0.25 - 0.25/l*R)**2 + (2.1 - 2.1/l*R)**2)
    print("level set value: %f, ref = %f" % (lsv, lsv_ref))
    assert(abs(lsv - lsv_ref < 1e-13))

    P = Point3D(0.0, 0.25, 2.1)
    l = 2.1
    lsv = ls.GetValue(P)
    lsv_ref = 2.1 - R
    print("level set value: %f, ref = %f" % (lsv, lsv_ref))
    assert(abs(lsv - lsv_ref < 1e-13))

    print("All tests passed")

if __name__ == "__main__":
    main()
