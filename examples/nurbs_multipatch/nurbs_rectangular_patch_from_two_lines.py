##################################################################
##################################################################
import sys
import os
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
mpatch_export = MultiNURBSPatchMatlabExporter()

## create line 1
fes1 = nurbs_fespace_library.CreateLinearFESpace(3)
ctrl_grid_1 = grid_lib.CreateLinearControlPointGrid(0.0, 0.0, 0.0, fes1.Number(0), 1.0, 0.0, 0.0)

patch1_ptr = multipatch_util.CreatePatchPointer(1, fes1)
patch1 = patch1_ptr.GetReference()
patch1.CreateControlPointGridFunction(ctrl_grid_1)
print(patch1)

## create line 2
fes2 = nurbs_fespace_library.CreateLinearFESpace(3)
ctrl_grid_2 = grid_lib.CreateLinearControlPointGrid(0.0, 1.0, 0.0, fes2.Number(0), 1.0, 1.0, 0.0)

patch2_ptr = multipatch_util.CreatePatchPointer(2, fes2)
patch2 = patch2_ptr.GetReference()
patch2.CreateControlPointGridFunction(ctrl_grid_2)
print(patch2)

## create rectangular patch by connect the two lines
rec_patch_ptr = bsplines_patch_util.CreateLoftPatch(patch1, patch2)
rec_patch = rec_patch_ptr.GetReference()
rec_patch.Id = 3
print(rec_patch)
mpatch_export.Export(rec_patch, "rec.m")


