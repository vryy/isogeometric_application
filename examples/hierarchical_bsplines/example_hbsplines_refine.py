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
multipatch_refine_util = MultiPatchRefinementUtility()

def CreatePatch():
    fn = "square.txt"
    patches_ptr = bsplines_patch_util.CreatePatchFromGeo(fn)

    multipatch_refine_util.DegreeElevate(patches_ptr[0], [1, 1])

    patch = patches_ptr[0].GetReference()
    patch.FESpace().Enumerate()
    #print(patch)

    hpatch_ptr = hbsplines_patch_util.CreatePatchFromBSplines(patch)
    hpatch = hpatch_ptr.GetReference()
    #print(hpatch)

    return hpatch_ptr

def main():
    patch_ptr = CreatePatch()
    patch = patch_ptr.GetReference()
    echo_level = 1
    hbsplines_refinement_util.Refine(patch, 1, echo_level)
    patch.FESpace().Enumerate()
    print(patch)

if __name__ == "__main__":
    main()

