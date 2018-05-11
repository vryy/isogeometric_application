##################################################################
# test the domain manager in 2D
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

bsplines_patch_util = BSplinesPatchUtility()
hbsplines_patch_util = HBSplinesPatchUtility()
hbsplines_refinement_util = HBSplinesRefinementUtility()
ctrl_grid_util = PointBasedControlGridUtility()

import example_hbsplines_refine

def main():
    patch_ptr = example_hbsplines_refine.CreatePatch()

    patch = patch_ptr.GetReference()
    temp_ctrl_grid = ctrl_grid_util.CreatePointBasedControlGrid(TEMPERATURE, patch.FESpace())
    for i in range(0, temp_ctrl_grid.Size()):
        temp_ctrl_grid[i] = 1.0
#    print(temp_ctrl_grid)
    temperature_grid_func = patch.CreateGridFunction(temp_ctrl_grid)
    print("temperature_grid_func:")
    print(temperature_grid_func)
    print("temperature_grid_func.GetValue(0.5, 0.5):", temperature_grid_func.GetValue([0.5, 0.5]))
##    print(patch)
    echo_level = 1
    refine_bf_id = 9
    hbsplines_refinement_util.Refine(patch, refine_bf_id, echo_level)
    patch.FESpace().Enumerate()
    print("temperature_grid_func:")
    print(temperature_grid_func)
    print("temperature_grid_func.GetValue(0.5, 0.5):", temperature_grid_func.GetValue([0.5, 0.5]))
#    print(patch)

if __name__ == "__main__":
    main()

