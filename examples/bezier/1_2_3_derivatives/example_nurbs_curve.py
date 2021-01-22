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
# multipatch_refine_util.InsertKnots(mpatch[1], [[0.4]])
# multipatch_refine_util.InsertKnots(mpatch[1], [[0.4]])
# multipatch_refine_util.InsertKnots(mpatch[1], [[0.4]])
# multipatch_refine_util.InsertKnots(mpatch[1], [[0.4]])
mpatch.Enumerate()
# print(mpatch)

mpatch_export3.Export(mpatch, "curve.m")

mpatch_mp = MultiPatchModelPart1D(mpatch)

mpatch_mp.BeginModelPart()
model_part = mpatch_mp.GetModelPart()

mpatch_mp.CreateNodes()

patch_ids = [1]
element_name = "DummyElementBezier"
for sid in patch_ids:
    patch_ptr = mpatch[sid]
    patch = patch_ptr.GetReference()

    ## add line elements
    last_elem_id = 0
    prop = model_part.Properties[1]
    elems = mpatch_mp.AddElements(patch, element_name, last_elem_id+1, prop)

mpatch_mp.EndModelPart()

element = model_part.Elements[1]
test_utils = IsogeometricTestUtils()
p = 0.4
test_utils.ProbeShapeFunctionValues(element, p)
print("--------------")
test_utils.ProbeShapeFunctionDerivatives(element, p)
print("--------------")
test_utils.ProbeShapeFunctionSecondDerivatives(element, p)
print("--------------")
test_utils.ProbeShapeFunctionThirdDerivatives(element, p)

### matlab code to verify the results:
# curve
# s=findspan(Patch_1.number-1, Patch_1.order-1, 0.4, Patch_1.knots)
# d=basisfunder(s, Patch_1.order-1, 0.4, Patch_1.knots, 3)
