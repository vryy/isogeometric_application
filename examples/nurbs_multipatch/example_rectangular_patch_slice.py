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

fes1 = nurbs_fespace_library.CreateRectangularFESpace(2, 1)
ctrl_grid_1 = grid_lib.CreateRectangularControlPointGrid(0.0, 0.0, fes1.Number(0), fes1.Number(1), 1.0, 1.0)
patch1_ptr = multipatch_util.CreatePatchPointer(1, fes1)
patch1 = patch1_ptr.GetReference()
patch1.CreateControlPointGridFunction(ctrl_grid_1)
print(patch1)

spatch = patch1.ConstructSlicedPatch(1, 0.3)
print(spatch.GetReference())
