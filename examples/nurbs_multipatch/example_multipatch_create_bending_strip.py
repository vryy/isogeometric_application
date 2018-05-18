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
mpatch_export1 = MultiNURBSPatchGeoExporter()
mpatch_export2 = MultiNURBSPatchGLVisExporter()
bending_strip_utility = BendingStripUtility()

def CreateMultiPatch():
    mpatch = MultiPatch2D()

    fes1 = nurbs_fespace_library.CreateRectangularFESpace(3, 3)
    ctrl_grid_1 = grid_lib.CreateRectangularControlPointGrid(0.0, 0.0, fes1.Number(0), fes1.Number(1), 1.0, 1.0)
#    print(ctrl_grid_1)
    patch1_ptr = multipatch_util.CreatePatchPointer(1, fes1)
    patch1 = patch1_ptr.GetReference()
    patch1.CreateControlPointGridFunction(ctrl_grid_1)
#    print(patch1)
    mpatch_export1.Export(patch1, "patch1.txt")

    fes2 = nurbs_fespace_library.CreateRectangularFESpace(3, 3)
    ctrl_grid_2 = grid_lib.CreateRectangularControlPointGrid(1.0, 0.0, fes1.Number(0), fes1.Number(1), 2.0, 1.0)
    patch2_ptr = multipatch_util.CreatePatchPointer(2, fes2)
    patch2 = patch2_ptr.GetReference()
    patch2.CreateControlPointGridFunction(ctrl_grid_2)
#    print(patch2)
    mpatch_export1.Export(patch2, "patch2.txt")

    mpatch.AddPatch(patch1_ptr)
    mpatch.AddPatch(patch2_ptr)
    patch1.Id = 1
    patch2.Id = 2
    mpatch.MakeNeighbor(patch1, BoundarySide.Right, patch2, BoundarySide.Left)
    #print(mpatch)

    print("############REFINEMENT###############")
    multipatch_refine_util = MultiPatchRefinementUtility()
    multipatch_refine_util.InsertKnots(patch1_ptr, [[0.5], [0.5]])
    patch1 = mpatch[1].GetReference()
    patch2 = mpatch[2].GetReference()
#    print(patch1)
#    print(patch2)
#    print(mpatch)
    system_size = mpatch.Enumerate()
    print("system_size:", system_size)
    #print("############RESULTS###############")

    print("############CREATE BENDING STRIP PATCH###############")
    strip_patch_id = 3
    strip_order = 2 # quadratic
#    print(patch1)
#    print(patch2)
    bspatch_ptr = bending_strip_utility.CreateBendingStripNURBSPatch(strip_patch_id, patch1, BoundarySide.Right, patch2, BoundarySide.Left, strip_order)
    mpatch.AddPatch(bspatch_ptr)

    return mpatch

######################################################################

def main():
    mpatch = CreateMultiPatch()
#    print(mpatch)

    bsp_ptr = mpatch[3]
    bsp = bsp_ptr.GetReference()
    print(bsp)

#    mpatch_export2.Export(mpatch, "mpatch.mesh")

#    print(mpatch)

if __name__ == "__main__":
    main()

#################RESULTS#####################(Validated with Matlab)
##>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<
##-------------Begin MultiPatchInfo-------------
##MultiPatch overview: Number of patches = 2
##MultiPatch details:
##-------------Begin PatchInfo-------------
##Patch2D, Id = 1, Add = 0x24416d0
##-------------Begin FESpaceInfo-------------
##NURBSFESpace2D, Add = 0x2440f90, n = ( 5 5), p = ( 3 3)
## knot vector 0: 0 0 0 0 0.5 1 1 1 1
## knot vector 1: 0 0 0 0 0.5 1 1 1 1
##-------------End FESpaceInfo-------------
##Grid CONTROL_POINT: [5, 5]
## Data:
## (
##  ( (0, 0, 0, 1) (0, 0.166667, 0, 1) (0, 0.5, 0, 1) (0, 0.833333, 0, 1) (0, 1, 0, 1))
##  ( (0.166667, 0, 0, 1) (0.166667, 0.166667, 0, 1) (0.166667, 0.5, 0, 1) (0.166667, 0.833333, 0, 1) (0.166667, 1, 0, 1))
##  ( (0.5, 0, 0, 1) (0.5, 0.166667, 0, 1) (0.5, 0.5, 0, 1) (0.5, 0.833333, 0, 1) (0.5, 1, 0, 1))
##  ( (0.833333, 0, 0, 1) (0.833333, 0.166667, 0, 1) (0.833333, 0.5, 0, 1) (0.833333, 0.833333, 0, 1) (0.833333, 1, 0, 1))
##  ( (1, 0, 0, 1) (1, 0.166667, 0, 1) (1, 0.5, 0, 1) (1, 0.833333, 0, 1) (1, 1, 0, 1))
## )
##Neighbors =  right:2
##-------------End PatchInfo-------------
##-------------Begin PatchInfo-------------
##Patch2D, Id = 2, Add = 0x2440e20
##-------------Begin FESpaceInfo-------------
##NURBSFESpace2D, Add = 0x24ca0a0, n = ( 4 5), p = ( 3 3)
## knot vector 0: 0 0 0 0 1 1 1 1
## knot vector 1: 0 0 0 0 0.5 1 1 1 1
##-------------End FESpaceInfo-------------
##Grid CONTROL_POINT: [4, 5]
## Data:
## (
##  ( (1, 0, 0, 1) (1, 0.166667, 0, 1) (1, 0.5, 0, 1) (1, 0.833333, 0, 1) (1, 1, 0, 1))
##  ( (1.33333, 0, 0, 1) (1.33333, 0.166667, 0, 1) (1.33333, 0.5, 0, 1) (1.33333, 0.833333, 0, 1) (1.33333, 1, 0, 1))
##  ( (1.66667, 0, 0, 1) (1.66667, 0.166667, 0, 1) (1.66667, 0.5, 0, 1) (1.66667, 0.833333, 0, 1) (1.66667, 1, 0, 1))
##  ( (2, 0, 0, 1) (2, 0.166667, 0, 1) (2, 0.5, 0, 1) (2, 0.833333, 0, 1) (2, 1, 0, 1))
## )
##Neighbors =  left:1
##-------------End PatchInfo-------------
##-------------End MultiPatchInfo-------------
##>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<


