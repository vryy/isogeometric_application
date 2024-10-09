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

fes1 = nurbs_fespace_library.CreateRectangularFESpace(3, 3)
ctrl_grid_1 = grid_lib.CreateRectangularControlPointGrid(0.0, 0.0, fes1.Number(0), fes1.Number(1), 1.0, 1.0)
print(fes1)
print(ctrl_grid_1)

patch1 = Patch2D(1, fes1)
patch1.CreateControlPointGridFunction(ctrl_grid_1)
print(patch1)

fes2 = nurbs_fespace_library.CreateCubicFESpace(3, 3, 2)
ctrl_grid_2 = grid_lib.CreateCubicControlPointGrid(0.0, 0.0, 0.0, fes2.Number(0), fes2.Number(1), fes2.Number(2), 1.0, 1.0, 1.0)
print(fes2)
print(ctrl_grid_2)

patch2 = Patch3D(2, fes2)
patch2.CreateControlPointGridFunction(ctrl_grid_2)
print(patch2)

