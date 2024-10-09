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

patches_ptr = patch_util.CreatePatchFromGeo("../../tests/infinite_plate.txt")
patch = patches_ptr[0].GetReference()
print(patch)

mpatch_export.Export(patch, "infinite_plate.m")

## To validate the results of the test, firstly the geometry shall be read by geo_load
# geom = geo_load('../../tests/infinite_plate.txt')
# load the nurbs geometry by typing infinite_plate
# check the differences: nurbs.coefs - geom.nurbs.coefs
# Results shall be:
#ans(:,:,1) =

#     0     0     0
#     0     0     0
#     0     0     0
#     0     0     0


#ans(:,:,2) =

#     0     0     0
#     0     0     0
#     0     0     0
#     0     0     0


#ans(:,:,3) =

#     0     0     0
#     0     0     0
#     0     0     0
#     0     0     0


#ans(:,:,4) =

#     0     0     0
#     0     0     0
#     0     0     0
#     0     0     0
