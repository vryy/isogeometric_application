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
    ctrl_grid_2 = grid_lib.CreateRectangularControlPointGrid(0.0, 1.0, fes1.Number(0), fes1.Number(1), 1.0, 2.0)
    patch2_ptr = multipatch_util.CreatePatchPointer(2, fes2)
    patch2 = patch2_ptr.GetReference()
    patch2.CreateControlPointGridFunction(ctrl_grid_2)
    #print(patch2)

    mpatch.AddPatch(patch1_ptr)
    mpatch.AddPatch(patch2_ptr)
    bsplines_patch_util.MakeInterface(patch1, BoundarySide.Top, patch2, BoundarySide.Bottom, BoundaryDirection.Forward)

    print("############REFINEMENT###############")
    multipatch_refine_util = MultiPatchRefinementUtility()

    multipatch_refine_util.DegreeElevate(patch1_ptr, [1, 1])

    #patch1 = patch1_ptr.GetReference()
    #patch2 = patch2_ptr.GetReference()

    return mpatch

def main():
    #print("############RESULTS###############")
    mpatch = CreateMultiPatch()
    mpatch.Enumerate()
    print(mpatch)
    mpatch_export1.Export(mpatch, "mpatch.mesh")
    mpatch_export2.Export(mpatch, "mpatch.m")

if __name__ == "__main__":
    main()

