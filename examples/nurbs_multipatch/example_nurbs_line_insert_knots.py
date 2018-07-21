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
mpatch_export2 = MultiNURBSPatchMatlabExporter()

fes1 = nurbs_fespace_library.CreateLinearFESpace(3)
ctrl_grid_1 = grid_lib.CreateLinearControlPointGrid(0.0, 0.0, 0.0, fes1.Number(0), 1.0, 0.0, 0.0)
#print(fes1)
#print(ctrl_grid_1)

patch1_ptr = multipatch_util.CreatePatchPointer(1, fes1)
patch1 = patch1_ptr.GetReference()
patch1.CreateControlPointGridFunction(ctrl_grid_1)
#print(patch1)

print("######BEFORE REFINEMENT#######")
#print(patch1)
mpatch_export1.Export(patch1, "line.txt")

print("######AFTER REFINEMENT#######")
mpatch = MultiPatch1D()
mpatch.AddPatch(patch1_ptr)
mpatch.Enumerate()

multipatch_refine_util = MultiPatchRefinementUtility()
multipatch_refine_util.InsertKnots(patch1_ptr, [[0.4, 0.5]])
mpatch.Enumerate()
print(mpatch) # result is the same as matlab code below
#mpatch_export2.Export(mpatch, "line.mesh")

## REFERENCE matlab code
##lin = nrbdegelev(nrbline,2)
##lin1 = nrbkntins(lin,[0.4 0.5])
##lin1.coefs


