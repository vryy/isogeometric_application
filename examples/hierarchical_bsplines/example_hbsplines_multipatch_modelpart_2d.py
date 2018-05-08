##################################################################
# test the domain manager in 2D
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
kernel = Kernel()   #defining kernel

bsplines_patch_util = BSplinesPatchUtility()
hbsplines_patch_util = HBSplinesPatchUtility()
hbsplines_refinement_util = HBSplinesRefinementUtility()

from example_hbsplines_refine import *

patch1_ptr = CreatePatch()
patch1 = patch1_ptr.GetReference()
#print(patch1)

mpatch = MultiPatch2D()
mpatch.AddPatch(patch1_ptr)
print(mpatch)

element_name = "KinematicLinearBezier2D"

mpatch_mp = MultiPatchModelPart2D(mpatch)
mpatch_mp.BeginModelPart()
mpatch_mp.CreateNodes()
mpatch_mp.AddElements(patch1, element_name, 1, 1)
mpatch_mp.EndModelPart()
print(mpatch_mp)

