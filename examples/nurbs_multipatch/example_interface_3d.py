from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.IsogeometricApplication import *

multipatch_util = MultiPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()
bsplines_patch_util = BSplinesPatchUtility()

import geometry_factory

def is_same(l1, l2):
    for i in range(0, len(l1)):
        if l1[i] != l2[i]:
            return False
    return True

def test():
    ref_interface_indices = {}
    ref_interface_indices[1] = {}
    ref_interface_indices[1][1] = [1, 3, 5, 7]
    ref_interface_indices[1][2] = [4, 6, 0, 2]
    ref_interface_indices[1][3] = [2, 0, 6, 4]
    ref_interface_indices[1][4] = [7, 5, 3, 1]
    ref_interface_indices[2] = {}
    ref_interface_indices[2][1] = [0, 4, 2, 6]
    ref_interface_indices[2][2] = [5, 1, 7, 3]
    ref_interface_indices[2][3] = [3, 7, 1, 5]
    ref_interface_indices[2][4] = [6, 2, 4, 0]
    ref_interface_indices[3] = {}
    ref_interface_indices[3][1] = [0, 1, 4, 5]
    ref_interface_indices[3][2] = [6, 7, 2, 3]
    ref_interface_indices[3][3] = [3, 2, 7, 6]
    ref_interface_indices[3][4] = [5, 4, 1, 0]
    ref_interface_indices[4] = {}
    ref_interface_indices[4][1] = [2, 6, 3, 7]
    ref_interface_indices[4][2] = [4, 0, 5, 1]
    ref_interface_indices[4][3] = [1, 5, 0, 4]
    ref_interface_indices[4][4] = [7, 3, 6, 2]
    ref_interface_indices[5] = {}
    ref_interface_indices[5][1] = [4, 5, 6, 7]
    ref_interface_indices[5][2] = [2, 3, 0, 1]
    ref_interface_indices[5][3] = [1, 0, 3, 2]
    ref_interface_indices[5][4] = [7, 6, 5, 4]
    ref_interface_indices[6] = {}
    ref_interface_indices[6][1] = [0, 2, 1, 3]
    ref_interface_indices[6][2] = [6, 4, 7, 5]
    ref_interface_indices[6][3] = [5, 7, 4, 6]
    ref_interface_indices[6][4] = [3, 1, 2, 0]

    settings = [1, 2, 3, 4, 5, 6]
    configurations = [1, 2, 3, 4]
    for setting in settings:
        for configuration in configurations:
            mpatch = geometry_factory.CreateTwoSlabs(setting=setting, configuration = configuration, dist = 0.0, L=1.0, w=1.0, h=1.0)
            mpatch.Enumerate()
            is_valid = mpatch.Validate()
            if not is_valid:
                raise Exception("MultiPatch is invalid")
            # mpatch_export.Export(mpatch, "two_cubes.m")
            multipatch_util.CheckInterfaces(mpatch, True) # check the indices on patch interface

            # bpatch_ptr = mpatch[2].GetReference().ConstructBoundaryPatch(BoundarySide3D.U0)
            # func_indices = bpatch_ptr.GetReference().FESpace().FunctionIndices()
            # # print(func_indices)
            # assert(is_same(func_indices, ref_interface_indices[setting][configuration]))

    print("Test passed")

if __name__ == '__main__':
    test()
