##################################################################
import sys
import os
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel
##################################################################

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
ctrl_grid_1.SetValue(2, ControlPoint(2.0, 0.0, 0.0, 1.0))
ctrl_grid_1.SetValue(3, ControlPoint(3.0, 0.0, 0.0, 1.0))
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
multipatch_refine_util.InsertKnots(mpatch[1], [[0.1]])
multipatch_refine_util.InsertKnots(mpatch[1], [[0.2]])
multipatch_refine_util.InsertKnots(mpatch[1], [[0.3]])
multipatch_refine_util.InsertKnots(mpatch[1], [[0.4]])
multipatch_refine_util.InsertKnots(mpatch[1], [[0.5]])
mpatch.Enumerate()
print(mpatch) # result is the same as matlab code below

patch = mpatch[1].GetReference()
gf = patch.GridFunction(CONTROL_POINT_COORDINATES)
print(gf.GetValue([0.0]))
print(gf.GetValue([0.1]))
print(gf.GetValue([0.2]))
print(gf.GetValue([0.3]))
print(gf.GetValue([0.4]))
print(gf.GetValue([1.0]))

mpatch_export3.Export(mpatch, "curve.m")

