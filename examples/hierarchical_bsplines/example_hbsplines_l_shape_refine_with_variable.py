##################################################################
# test the domain manager in 2D
##################################################################
#importing Kratos modules
import sys
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

bsplines_patch_util = BSplinesPatchUtility()
hbsplines_patch_util = HBSplinesPatchUtility()
hbsplines_refinement_util = HBSplinesRefinementUtility()
ctrl_grid_util = PointBasedControlGridUtility()

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()

hmpatch_export = MultiHBSplinesPatchMatlabExporter()

sys.path.append("../nurbs_multipatch")
import nurbs_l_shape

def main():
    mpatch = nurbs_l_shape.CreateMultiPatch()
    mpatch.Enumerate()

    ###############

    patch1_ptr = mpatch[1]
    patch1 = patch1_ptr.GetReference()
    hpatch1_ptr = hbsplines_patch_util.CreatePatchFromBSplines(patch1)
    hpatch1 = hpatch1_ptr.GetReference()
#    print(hpatch1)

    patch2_ptr = mpatch[2]
    patch2 = patch2_ptr.GetReference()
    hpatch2_ptr = hbsplines_patch_util.CreatePatchFromBSplines(patch2)
    hpatch2 = hpatch2_ptr.GetReference()
#    print(hpatch2)

    patch3_ptr = mpatch[3]
    patch3 = patch3_ptr.GetReference()
    hpatch3_ptr = hbsplines_patch_util.CreatePatchFromBSplines(patch3)
    hpatch3 = hpatch3_ptr.GetReference()
#    print(hpatch2)

    hmpatch = MultiPatch2D()
    hmpatch.AddPatch(hpatch1_ptr)
    hmpatch.AddPatch(hpatch2_ptr)
    hmpatch.AddPatch(hpatch3_ptr)
    multipatch_util.MakeInterface(hpatch1, BoundarySide.Top, hpatch2, BoundarySide.Bottom)
    multipatch_util.MakeInterface(hpatch2, BoundarySide.Right, hpatch3, BoundarySide.Left)
    hmpatch.Enumerate()
    print("Create hierarchical B-Splines multipatch completed")

#    print(hmpatch)

    hpatch1_ptr = hmpatch[1]
    hpatch1 = hpatch1_ptr.GetReference()

    hpatch2_ptr = hmpatch[2]
    hpatch2 = hpatch2_ptr.GetReference()

    hpatch3_ptr = hmpatch[3]
    hpatch3 = hpatch3_ptr.GetReference()

    temp_ctrl_grid = ctrl_grid_util.CreatePointBasedControlGrid(TEMPERATURE, hpatch1.FESpace())
    for i in range(0, temp_ctrl_grid.Size()):
        temp_ctrl_grid[i] = 1.0
#    print(temp_ctrl_grid)
    temperature_grid_func = hpatch1.CreateGridFunction(temp_ctrl_grid)
    print(temperature_grid_func)
    print("temperature_grid_func.GetValue(0.5, 0.5):", temperature_grid_func.GetValue([0.5, 0.5]))

    ###############

    echo_level = IsogeometricEchoFlags.ECHO_REFINEMENT

    for step in range(0, 3):
        corner_bf_list = hpatch1.FESpace().GetBoundaryBfs(BoundaryFlag.Right + BoundaryFlag.Top)
        refine_bf_id = corner_bf_list[0].Id
        print("refine_bf_id:", refine_bf_id)
        hbsplines_refinement_util.Refine(hpatch1, refine_bf_id, echo_level)
        hmpatch.Enumerate()

    ###

    ################

    boundary_basis_1 = hpatch1.FESpace().GetBoundaryBfs(BoundaryFlag.Top)
    boundary_basis_2 = hpatch2.FESpace().GetBoundaryBfs(BoundaryFlag.Bottom)

    for i in range(0, len(boundary_basis_1)):
        bf1 = boundary_basis_1[i]
        bf2 = boundary_basis_2[i]
        print(bf1.Id, bf1.EquationId)
        print(bf2.Id, bf2.EquationId)
#        print(bf1)
#        print(bf2)
        print("-----------")

    boundary_basis_3 = hpatch2.FESpace().GetBoundaryBfs(BoundaryFlag.Right)
    boundary_basis_4 = hpatch3.FESpace().GetBoundaryBfs(BoundaryFlag.Left)

    for i in range(0, len(boundary_basis_3)):
        bf1 = boundary_basis_3[i]
        bf2 = boundary_basis_4[i]
        print(bf1.Id, bf1.EquationId)
        print(bf2.Id, bf2.EquationId)
#        print(bf1)
#        print(bf2)
        print("------------------")

    ################

    print(temperature_grid_func)
    print("temperature_grid_func.GetValue(0.5, 0.5):", temperature_grid_func.GetValue([0.5, 0.5]))

    ################

    bfespace = hpatch1.FESpace().ConstructBoundaryFESpace(BoundarySide.Top)
    print("patch 1 boundary FESpace on top:")
    print(bfespace)

    ################

    hmpatch_export.Export(hmpatch, "hierarchical_b_splines_l_shape.m")

if __name__ == "__main__":
    main()

