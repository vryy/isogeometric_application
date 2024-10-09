##################################################################
##################################################################
import sys
import os
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
kernel = Kernel()   #defining kernel

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()

import example_multipatch_insert_knots

mpatch = example_multipatch_insert_knots.CreateMultiPatch()

patch1_ptr = mpatch[1]
patch1 = patch1_ptr.GetReference()
patch2_ptr = mpatch[2]
patch2 = patch2_ptr.GetReference()

element_name = "KinematicLinearBezier2D"

mpatch_mp = MultiPatchModelPart2D(mpatch)
mpatch_mp.BeginModelPart()
mpatch_mp.CreateNodes()
mpatch_mp.AddElements(patch1, element_name, 1, 1)
mpatch_mp.AddElements(patch2, element_name, 100, 1)
mpatch_mp.EndModelPart()
print(mpatch_mp)

