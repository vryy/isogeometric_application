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

import geometry_factory

def CreateMultiPatch(r):
    ####### create arc 1
    arc1_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r, -45.0, 45.0)
    arc1 = arc1_ptr.GetReference()
    arc1.Id = 1

    ####### create arc 2
    arc2_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r, 45.0, 135.0)
    arc2 = arc2_ptr.GetReference()
    arc2.Id = 2

    ####### create arc 3
    arc3_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r, 135.0, 225.0)
    arc3 = arc3_ptr.GetReference()
    arc3.Id = 3

    ####### create arc 4
    arc4_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r, 225.0, 315.0)
    arc4 = arc4_ptr.GetReference()
    arc4.Id = 4

    ###### create multipatch
    mpatch = MultiPatch1D()
    mpatch.AddPatch(arc1_ptr)
    mpatch.AddPatch(arc2_ptr)
    mpatch.AddPatch(arc3_ptr)
    mpatch.AddPatch(arc4_ptr)
    bsplines_patch_util.MakeInterface(arc1, BoundarySide.Right, arc2, BoundarySide.Left)
    bsplines_patch_util.MakeInterface(arc2, BoundarySide.Right, arc3, BoundarySide.Left)
    bsplines_patch_util.MakeInterface(arc3, BoundarySide.Right, arc4, BoundarySide.Left)
    bsplines_patch_util.MakeInterface(arc4, BoundarySide.Right, arc1, BoundarySide.Left)

    return mpatch

######################################################################

def main():
    mpatch = CreateMultiPatch(1.0)
    print("Multipatch is created")

    # mpatch_export2.Export(mpatch, "mpatch.mesh")

    system_size = mpatch.Enumerate()
    print("system_size:", system_size)
    # print(mpatch)

    multipatch_util.CheckInterfaces(mpatch)

    ls = MultiPatchNormalLevelSet1D(mpatch)
    ls.SetEchoLevel(0)
    # print(ls)

    P = Point3D(1.6, 0.0, 0.0)
    # print("level set value: " + str(ls.GetValue(P)))
    assert(abs(ls.GetValue(P) - 0.6) < 1e-13)

    P = Point3D(0.6, 0.6, 0.0)
    lsv = ls.GetValue(P)
    lsv_ref = math.sqrt(2.0)*(0.6-1.0/math.sqrt(2.0))
    # print("level set value: " + str(lsv))
    # print("lsv_ref: " + str(lsv_ref))
    assert(abs(lsv - lsv_ref) < 1e-13)

if __name__ == "__main__":
    main()
