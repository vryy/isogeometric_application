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
    order_u = 2
    order_v = 2
    l = 500.0

    fes1 = nurbs_fespace_library.CreateRectangularFESpace(order_u, order_v)
    ctrl_grid_1 = grid_lib.CreateRectangularControlPointGrid(0.0, 0.0, fes1.Number(0), fes1.Number(1), l, l)
#    print(ctrl_grid_1)
    id1 = 1
    patch1_ptr = multipatch_util.CreatePatchPointer(id1, fes1)
    patch1 = patch1_ptr.GetReference()
    patch1.CreateControlPointGridFunction(ctrl_grid_1)
#    print(patch1)
#    mpatch_export1.Export(patch1, "patch1.txt")

    fes2 = nurbs_fespace_library.CreateRectangularFESpace(order_u, order_v)
    ctrl_grid_2 = grid_lib.CreateRectangularControlPointGrid(0.0, l, fes2.Number(0), fes2.Number(1), l, 2*l)
    id2 = 2
    patch2_ptr = multipatch_util.CreatePatchPointer(id2, fes2)
    patch2 = patch2_ptr.GetReference()
    patch2.CreateControlPointGridFunction(ctrl_grid_2)
#    print(patch2)
#    mpatch_export1.Export(patch2, "patch2.txt")

    fes3 = nurbs_fespace_library.CreateRectangularFESpace(order_u, order_v)
    ctrl_grid_3 = grid_lib.CreateRectangularControlPointGrid(l, l, fes3.Number(0), fes3.Number(1), 2*l, 2*l)
    id3 = 3
    patch3_ptr = multipatch_util.CreatePatchPointer(id3, fes3)
    patch3 = patch3_ptr.GetReference()
    patch3.CreateControlPointGridFunction(ctrl_grid_3)
#    print(patch3)
#    mpatch_export.Export(patch3, "patch3.txt")

    mpatch = MultiPatch2D()
    mpatch.AddPatch(patch1_ptr)
    mpatch.AddPatch(patch2_ptr)
    mpatch.AddPatch(patch3_ptr)
    patch1.Id = 1
    patch2.Id = 2
    patch3.Id = 3
    bsplines_patch_util.MakeInterface(patch1, BoundarySide.Top, patch2, BoundarySide.Bottom, BoundaryDirection.Forward)
    bsplines_patch_util.MakeInterface(patch2, BoundarySide.Right, patch3, BoundarySide.Left, BoundaryDirection.Forward)

    return mpatch

def Refine(mpatch):
    print("############REFINEMENT###############")
    multipatch_refine_util = MultiPatchRefinementUtility()
    patch1_ptr = mpatch[1]
    multipatch_refine_util.InsertKnots(patch1_ptr, [[0.25, 0.5, 0.75], [0.5]])

    patch2_ptr = mpatch[2]
    multipatch_refine_util.InsertKnots(patch2_ptr, [[], [0.5]])

    return mpatch

def main():
    mpatch = CreateMultiPatch()
    mpatch = Refine(mpatch)
    mpatch.Enumerate()
    print(mpatch)

    mpatch_export.Export(mpatch, "l_shape.m")

if __name__ == "__main__":
    main()


