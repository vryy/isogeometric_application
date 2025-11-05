from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
multipatch_refine_util = MultiPatchRefinementUtility()
bsplines_patch_util = BSplinesPatchUtility()

def CreateLine(start_point, end_point, order = 1):
    Id = 0
    fes = nurbs_fespace_library.CreatePrimitiveFESpace(order)
    ctrl_grid = grid_lib.CreateLinearControlPointGrid(start_point[0], start_point[1], start_point[2], fes.Number(0), end_point[0], end_point[1], end_point[2])
    # patch_ptr = Patch1DSelector.RealPatch.Create(Id, fes)
    patch_ptr = Patch1DSelector.ComplexPatch.Create(Id, fes)
    patch = patch_ptr.GetReference()
    patch.CreateControlPointGridFunction(ctrl_grid)
    return patch_ptr

start_point = [0.0, 0.0, 0.0]
end_point = [1.0, 1.0, 1.0]
line_ptr = CreateLine(start_point, end_point)
line = line_ptr.GetReference()
print(line)

# Patch1DSelector.RealPatch.Create(Id, fes)
# Patch1DSelector.ComplexPatch
