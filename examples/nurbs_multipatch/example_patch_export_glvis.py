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
mpatch = MultiPatch2D()
mpatch_export = MultiNURBSPatchGLVisExporter()

fes1 = nurbs_fespace_library.CreateRectangularFESpace(1, 1)
ctrl_grid_1 = grid_lib.CreateRectangularControlPointGrid(0.0, 0.0, fes1.Number(0), fes1.Number(1), 1.0, 1.0)
patch1_ptr = multipatch_util.CreatePatchPointer(1, fes1)
patch1 = patch1_ptr.GetReference()
patch1.CreateControlPointGridFunction(ctrl_grid_1)
#print(patch1)

mpatch.AddPatch(patch1_ptr)

print("############REFINEMENT###############")
multipatch_refine_util = MultiPatchRefinementUtility()
#multipatch_refine_util.InsertKnots(patch1_ptr, [[0.5], [0.5]])
multipatch_refine_util.DegreeElevate(patch1_ptr, [0, 1])
mpatch.Enumerate()
print(mpatch)
#patch1 = patch1_ptr.GetReference()
#patch2 = patch2_ptr.GetReference()
#print("############RESULTS###############")
mpatch_export.Export(mpatch, "mpatch.mesh")

