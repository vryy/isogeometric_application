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

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()

hmpatch_export = MultiHBSplinesPatchMatlabExporter()

def CreateMultiPatch():

    ### create B-Splines multipatch

    fes1 = nurbs_fespace_library.CreateRectangularFESpace(3, 3)
    ctrl_grid_1 = grid_lib.CreateRectangularControlPointGrid(0.0, 0.0, fes1.Number(0), fes1.Number(1), 1.0, 1.0)
    patch1_ptr = multipatch_util.CreatePatchPointer(1, fes1)
    patch1 = patch1_ptr.GetReference()
    patch1.CreateControlPointGridFunction(ctrl_grid_1)
    #print(patch1)

    fes2 = nurbs_fespace_library.CreateRectangularFESpace(3, 3)
    ctrl_grid_2 = grid_lib.CreateRectangularControlPointGrid(1.2, 0.0, fes1.Number(0), fes1.Number(1), 2.2, 1.0)
    patch2_ptr = multipatch_util.CreatePatchPointer(2, fes2)
    patch2 = patch2_ptr.GetReference()
    patch2.CreateControlPointGridFunction(ctrl_grid_2)
    #print(patch2)

    mpatch = MultiPatch2D()
    mpatch.AddPatch(patch1_ptr)
    mpatch.AddPatch(patch2_ptr)
    multipatch_util.MakeInterface(patch1, BoundarySide.Right, patch2, BoundarySide.Left)
    mpatch.Enumerate()

#    ### create hierarchical B-Splines multipatch

    hpatch1_ptr = hbsplines_patch_util.CreatePatchFromBSplines(patch1)
    hpatch1 = hpatch1_ptr.GetReference()
#    print(hpatch1)

    hpatch2_ptr = hbsplines_patch_util.CreatePatchFromBSplines(patch2)
    hpatch2 = hpatch2_ptr.GetReference()
#    print(hpatch2)

    hmpatch = MultiPatch2D()
    hmpatch.AddPatch(hpatch1_ptr)
    hmpatch.AddPatch(hpatch2_ptr)
    multipatch_util.MakeInterface(hpatch1, BoundarySide.Right, hpatch2, BoundarySide.Left)
    hmpatch.Enumerate()
    print("Create multipatch completed")

    return hmpatch

def main():
    hmpatch = CreateMultiPatch()
#    print(hmpatch)

    hpatch1_ptr = hmpatch[1]
    hpatch1 = hpatch1_ptr.GetReference()

    hpatch2_ptr = hmpatch[2]
    hpatch2 = hpatch2_ptr.GetReference()

    ###############

    echo_level = IsogeometricEchoFlags.ECHO_REFINEMENT
    refine_bf_id = 4
    hbsplines_refinement_util.Refine(hpatch1, refine_bf_id, echo_level)
    hmpatch.Enumerate()
#    print(hmpatch)

    ###

    refine_bf_id = 20
    hbsplines_refinement_util.Refine(hpatch1, refine_bf_id, echo_level)
    hmpatch.Enumerate()

    ###

    refine_bf_id = 26
    hbsplines_refinement_util.Refine(hpatch1, refine_bf_id, echo_level)
    hmpatch.Enumerate()

    ################

    boundary_basis_1 = hpatch1.FESpace().GetBoundaryBfs(BoundaryFlag.Right)
    boundary_basis_2 = hpatch2.FESpace().GetBoundaryBfs(BoundaryFlag.Left)

    for i in range(0, len(boundary_basis_1)):
        bf1 = boundary_basis_1[i]
        bf2 = boundary_basis_2[i]
        print(bf1.Id, bf1.EquationId)
        print(bf2.Id, bf2.EquationId)
#        print(bf1)
#        print(bf2)
        print("-----------")

    hmpatch_export.Export(hmpatch, "hierarchical_b_splines_two_rectangles.m")

if __name__ == "__main__":
    main()

