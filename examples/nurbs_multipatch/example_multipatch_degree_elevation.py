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
import pdb

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
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export1 = MultiNURBSPatchGLVisExporter()
mpatch_export2 = MultiNURBSPatchMatlabExporter()
mpatch = MultiPatch2D()

def CreateMultiPatch():

    fes1 = nurbs_fespace_library.CreateRectangularFESpace(3, 3)
    ctrl_grid_1 = grid_lib.CreateRectangularControlPointGrid(0.0, 0.0, fes1.Number(0), fes1.Number(1), 1.0, 1.0)
    patch1_ptr = multipatch_util.CreatePatchPointer(1, fes1)
    patch1 = patch1_ptr.GetReference()
    patch1.CreateControlPointGridFunction(ctrl_grid_1)
    #print(patch1)

    fes2 = nurbs_fespace_library.CreateRectangularFESpace(3, 3)
    ctrl_grid_2 = grid_lib.CreateRectangularControlPointGrid(1.0, 0.0, fes1.Number(0), fes1.Number(1), 2.0, 1.0)
    patch2_ptr = multipatch_util.CreatePatchPointer(2, fes2)
    patch2 = patch2_ptr.GetReference()
    patch2.CreateControlPointGridFunction(ctrl_grid_2)
    #print(patch2)

    mpatch.AddPatch(patch1_ptr)
    mpatch.AddPatch(patch2_ptr)
    bsplines_patch_util.MakeInterface(patch1, BoundarySide.Right, patch2, BoundarySide.Left, BoundaryDirection.Forward)

    print("############REFINEMENT###############")
    multipatch_refine_util = MultiPatchRefinementUtility()

    multipatch_refine_util.DegreeElevate(patch1_ptr, [1, 1])

    #patch1 = patch1_ptr.GetReference()
    #patch2 = patch2_ptr.GetReference()
    print("new patch 1 is at address " + str(mpatch[1]))
    print("new patch 2 is at address " + str(mpatch[2]))
    print("############REFINEMENT COMPLETED###############")

    return mpatch

#print("############RESULTS###############")
def main():
#    pdb.set_trace()
    mpatch = CreateMultiPatch()
    mpatch.Enumerate()
    print(mpatch)
    mpatch_export1.Export(mpatch, "mpatch.mesh")
    mpatch_export2.Export(mpatch, "mpatch.m")

if __name__ == "__main__":
    main()

#################RESULTS#####################(Validated with Matlab)
##>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<
##-------------Begin MultiPatchInfo-------------
##MultiPatch overview: Number of patches = 2
##MultiPatch details:
##-------------Begin PatchInfo-------------
##Patch2D, Id = 1, Add = 0x12d56d0
##-------------Begin FESpaceInfo-------------
##NURBSFESpace2D, Add = 0x12d4f90, n = ( 5 5), p = ( 4 4)
## knot vector 0: 0 0 0 0 0 1 1 1 1 1
## knot vector 1: 0 0 0 0 0 1 1 1 1 1
##-------------End FESpaceInfo-------------
##Grid CONTROL_POINT: [5, 5]
## Data:
## (
##  ( (0, 0, 0, 1) (0, 0.25, 0, 1) (0, 0.5, 0, 1) (0, 0.75, 0, 1) (0, 1, 0, 1))
##  ( (0.25, 0, 0, 1) (0.25, 0.25, 0, 1) (0.25, 0.5, 0, 1) (0.25, 0.75, 0, 1) (0.25, 1, 0, 1))
##  ( (0.5, 0, 0, 1) (0.5, 0.25, 0, 1) (0.5, 0.5, 0, 1) (0.5, 0.75, 0, 1) (0.5, 1, 0, 1))
##  ( (0.75, 0, 0, 1) (0.75, 0.25, 0, 1) (0.75, 0.5, 0, 1) (0.75, 0.75, 0, 1) (0.75, 1, 0, 1))
##  ( (1, 0, 0, 1) (1, 0.25, 0, 1) (1, 0.5, 0, 1) (1, 0.75, 0, 1) (1, 1, 0, 1))
## )
##Neighbors =  right:2
##-------------End PatchInfo-------------
##-------------Begin PatchInfo-------------
##Patch2D, Id = 2, Add = 0x12d4e20
##-------------Begin FESpaceInfo-------------
##NURBSFESpace2D, Add = 0x135df50, n = ( 4 5), p = ( 3 4)
## knot vector 0: 0 0 0 0 1 1 1 1
## knot vector 1: 0 0 0 0 0 1 1 1 1 1
##-------------End FESpaceInfo-------------
##Grid CONTROL_POINT: [4, 5]
## Data:
## (
##  ( (1, 0, 0, 1) (1, 0.25, 0, 1) (1, 0.5, 0, 1) (1, 0.75, 0, 1) (1, 1, 0, 1))
##  ( (1.33333, 0, 0, 1) (1.33333, 0.25, 0, 1) (1.33333, 0.5, 0, 1) (1.33333, 0.75, 0, 1) (1.33333, 1, 0, 1))
##  ( (1.66667, 0, 0, 1) (1.66667, 0.25, 0, 1) (1.66667, 0.5, 0, 1) (1.66667, 0.75, 0, 1) (1.66667, 1, 0, 1))
##  ( (2, 0, 0, 1) (2, 0.25, 0, 1) (2, 0.5, 0, 1) (2, 0.75, 0, 1) (2, 1, 0, 1))
## )
##Neighbors =  left:1
##-------------End PatchInfo-------------
##-------------End MultiPatchInfo-------------
##>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<


