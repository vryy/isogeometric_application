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
mpatch_export3 = MultiNURBSPatchMatlabExporter()

print("######CREATE PATCH#######")
fes1 = nurbs_fespace_library.CreateLinearFESpace(3)
ctrl_grid_1 = StructuredControlPointGrid1D(4)
#print(fes1)
ctrl_grid_1.SetValue(0, ControlPoint(0.0, 0.0, 0.0, 1.0))
ctrl_grid_1.SetValue(1, ControlPoint(1.0, 0.0, 0.0, 1.0))
ctrl_grid_1.SetValue(2, ControlPoint(1.0, 1.0, 0.0, 1.0))
ctrl_grid_1.SetValue(3, ControlPoint(2.0, 1.0, 0.0, 1.0))
#print(ctrl_grid_1)

patch1_ptr = multipatch_util.CreatePatchPointer(1, fes1)
patch1 = patch1_ptr.GetReference()
patch1.CreateControlPointGridFunction(ctrl_grid_1)
patch1.Id = 1

print("######CREATE MULTIPATCH#######")
mpatch = MultiPatch1D()
mpatch.AddPatch(patch1_ptr)
mpatch.Enumerate()

print("######REFINEMENT#######")
multipatch_refine_util = MultiPatchRefinementUtility()
multipatch_refine_util.InsertKnots(mpatch[1], [[0.4]])
multipatch_refine_util.InsertKnots(mpatch[1], [[0.4]])
multipatch_refine_util.InsertKnots(mpatch[1], [[0.4]])
multipatch_refine_util.InsertKnots(mpatch[1], [[0.4]])
mpatch.Enumerate()
print(mpatch) # result is the same as matlab code below

mpatch_export3.Export(mpatch, "curve.m")

