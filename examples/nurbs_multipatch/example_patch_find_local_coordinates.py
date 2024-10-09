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

patch_util = BSplinesPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()

patches_ptr = patch_util.CreatePatchFromGeo("../../tests/geo_ring.txt")
patch = patches_ptr[0].GetReference()
# print(patch)

grid_func = patch.GridFunction(CONTROL_POINT_COORDINATES)
[stat, xi] = patch.LocalCoordinates([1.0, 1.0, 0.0], [0.5, 0.5, 0.0])
print("stat:", stat)
print("xi:", str(xi))
P1 = grid_func.GetValue([xi[0], xi[1]])
print("P1:", str(P1))

print("([1.0, 1.0, 0.0] is inside:", patch.IsInside([1.0, 1.0, 0.0], [0.5, 0.5, 0.0]))