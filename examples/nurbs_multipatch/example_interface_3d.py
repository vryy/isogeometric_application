from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.IsogeometricApplication import *

multipatch_util = MultiPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()
bsplines_patch_util = BSplinesPatchUtility()

import geometry_factory

def CreateMultiPatch(setting=1, configuration = 1, dist = 0.1):
    if setting == 1:
        return CreateMultiPatch1(configuration = configuration, dist = dist)
    elif setting == 2:
        return CreateMultiPatch2(configuration = configuration, dist = dist)
    elif setting == 3:
        return CreateMultiPatch3(configuration = configuration, dist = dist)
    elif setting == 4:
        return CreateMultiPatch4(configuration = configuration, dist = dist)
    elif setting == 5:
        return CreateMultiPatch5(configuration = configuration, dist = dist)
    elif setting == 6:
        return CreateMultiPatch6(configuration = configuration, dist = dist)
    else:
        raise Exception("Unknown setting %d", setting)

## in this setting, the two surface patch matches on v->v and w->w, that means uv_or_vu == True
def CreateMultiPatch1(configuration = 1, dist = 0.1):

    # create patch 1
    p1 = [0.0, 0.0, 0.0]
    p2 = [1.0, 1.0, 1.0]
    patch1_ptr = geometry_factory.CreateSlab(p1, p2)
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1
    patch1.LayerIndex = 1

    # create patch 2
    p3 = [1.0 + dist, 0.0, 0.0]
    p4 = [2.0 + dist, 1.0, 1.0]
    patch2_ptr = geometry_factory.CreateSlab(p3, p4)
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2
    patch2.LayerIndex = 1

    if configuration == 1:

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U1, patch2, BoundarySide3D.U0, True, BoundaryDirection.Forward, BoundaryDirection.Forward)

    elif configuration == 2:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U0, patch2, BoundarySide3D.U0, True, BoundaryDirection.Forward, BoundaryDirection.Reversed)

    elif configuration == 3:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 1)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U0, patch2, BoundarySide3D.U0, True, BoundaryDirection.Reversed, BoundaryDirection.Forward)

    elif configuration == 4:

        bsplines_patch_util.Reverse(patch1, 1)
        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U1, patch2, BoundarySide3D.U0, True, BoundaryDirection.Reversed, BoundaryDirection.Reversed)

    else:
        raise Exception("Unknown configuration %d" % (configuration))

    return mpatch

## in this setting, the two surface patch matches on v->w and w->v, that means uv_or_vu == False
def CreateMultiPatch2(configuration = 1, dist = 0.1):

    # create patch 1
    p1 = [0.0, 0.0, 0.0]
    d1 = [1.0, 0.0, 0.0]
    d2 = [0.0, 0.0, 1.0]
    d3 = [0.0, 1.0, 0.0]
    patch1_ptr = geometry_factory.CreateParallelepiped(p1, d1, d2, d3)
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1
    patch1.LayerIndex = 1

    # create patch 2
    p3 = [1.0 + dist, 0.0, 0.0]
    p4 = [2.0 + dist, 1.0, 1.0]
    patch2_ptr = geometry_factory.CreateSlab(p3, p4)
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2
    patch2.LayerIndex = 1

    if configuration == 1:

        bsplines_patch_util.Reverse(patch1, 0)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U0, patch2, BoundarySide3D.U0, False, BoundaryDirection.Forward, BoundaryDirection.Forward)

    elif configuration == 2:

        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U1, patch2, BoundarySide3D.U0, False, BoundaryDirection.Forward, BoundaryDirection.Reversed)

    elif configuration == 3:

        bsplines_patch_util.Reverse(patch1, 1)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U1, patch2, BoundarySide3D.U0, False, BoundaryDirection.Reversed, BoundaryDirection.Forward)

    elif configuration == 4:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 1)
        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.U0, patch2, BoundarySide3D.U0, False, BoundaryDirection.Reversed, BoundaryDirection.Reversed)

    else:
        raise Exception("Unknown configuration %d" % (configuration))

    return mpatch

## in this setting, the two surface patch matches on u->v and w->w, that means uv_or_vu == True
def CreateMultiPatch3(configuration = 1, dist = 0.1):

    # create patch 1
    p1 = [0.0, 0.0, 0.0]
    d1 = [0.0, 1.0, 0.0]
    d2 = [1.0, 0.0, 0.0]
    d3 = [0.0, 0.0, 1.0]
    patch1_ptr = geometry_factory.CreateParallelepiped(p1, d1, d2, d3)
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1
    patch1.LayerIndex = 1

    # create patch 2
    p3 = [1.0 + dist, 0.0, 0.0]
    p4 = [2.0 + dist, 1.0, 1.0]
    patch2_ptr = geometry_factory.CreateSlab(p3, p4)
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2
    patch2.LayerIndex = 1

    if configuration == 1:

        bsplines_patch_util.Reverse(patch1, 1)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V0, patch2, BoundarySide3D.U0, True, BoundaryDirection.Forward, BoundaryDirection.Forward)

    elif configuration == 2:

        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V1, patch2, BoundarySide3D.U0, True, BoundaryDirection.Forward, BoundaryDirection.Reversed)

    elif configuration == 3:

        bsplines_patch_util.Reverse(patch1, 0)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V1, patch2, BoundarySide3D.U0, True, BoundaryDirection.Reversed, BoundaryDirection.Forward)

    elif configuration == 4:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 1)
        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V0, patch2, BoundarySide3D.U0, True, BoundaryDirection.Reversed, BoundaryDirection.Reversed)

    else:
        raise Exception("Unknown configuration %d" % (configuration))

    return mpatch

## in this setting, the two surface patch matches on w->v and u->w, that means uv_or_vu == False
def CreateMultiPatch4(configuration = 1, dist = 0.1):

    # create patch 1
    p1 = [0.0, 0.0, 0.0]
    d1 = [0.0, 0.0, 1.0]
    d2 = [1.0, 0.0, 0.0]
    d3 = [0.0, 1.0, 0.0]
    patch1_ptr = geometry_factory.CreateParallelepiped(p1, d1, d2, d3)
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1
    patch1.LayerIndex = 1

    # create patch 2
    p3 = [1.0 + dist, 0.0, 0.0]
    p4 = [2.0 + dist, 1.0, 1.0]
    patch2_ptr = geometry_factory.CreateSlab(p3, p4)
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2
    patch2.LayerIndex = 1

    if configuration == 1:

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V1, patch2, BoundarySide3D.U0, False, BoundaryDirection.Forward, BoundaryDirection.Forward)

    elif configuration == 2:

        bsplines_patch_util.Reverse(patch1, 2)
        bsplines_patch_util.Reverse(patch1, 1)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V0, patch2, BoundarySide3D.U0, False, BoundaryDirection.Forward, BoundaryDirection.Reversed)

    elif configuration == 3:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 1)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V0, patch2, BoundarySide3D.U0, False, BoundaryDirection.Reversed, BoundaryDirection.Forward)

    elif configuration == 4:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.V1, patch2, BoundarySide3D.U0, False, BoundaryDirection.Reversed, BoundaryDirection.Reversed)

    else:
        raise Exception("Unknown configuration %d" % (configuration))

    return mpatch

## in this setting, the two surface patch matches on u->v and v->w, that means uv_or_vu == True
def CreateMultiPatch5(configuration = 1, dist = 0.1):

    # create patch 1
    p1 = [0.0, 0.0, 0.0]
    d1 = [0.0, 1.0, 0.0]
    d2 = [0.0, 0.0, 1.0]
    d3 = [1.0, 0.0, 0.0]
    patch1_ptr = geometry_factory.CreateParallelepiped(p1, d1, d2, d3)
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1
    patch1.LayerIndex = 1

    # create patch 2
    p3 = [1.0 + dist, 0.0, 0.0]
    p4 = [2.0 + dist, 1.0, 1.0]
    patch2_ptr = geometry_factory.CreateSlab(p3, p4)
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2
    patch2.LayerIndex = 1

    if configuration == 1:

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W1, patch2, BoundarySide3D.U0, True, BoundaryDirection.Forward, BoundaryDirection.Forward)

    elif configuration == 2:

        bsplines_patch_util.Reverse(patch1, 1)
        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W0, patch2, BoundarySide3D.U0, True, BoundaryDirection.Forward, BoundaryDirection.Reversed)

    elif configuration == 3:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W0, patch2, BoundarySide3D.U0, True, BoundaryDirection.Reversed, BoundaryDirection.Forward)

    elif configuration == 4:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 1)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W1, patch2, BoundarySide3D.U0, True, BoundaryDirection.Reversed, BoundaryDirection.Reversed)

    else:
        raise Exception("Unknown configuration %d" % (configuration))

    return mpatch

## in this setting, the two surface patch matches on v->v and u->w, that means uv_or_vu == False
def CreateMultiPatch6(configuration = 1, dist = 0.1):

    # create patch 1
    p1 = [0.0, 0.0, 0.0]
    d1 = [0.0, 0.0, 1.0]
    d2 = [0.0, 1.0, 0.0]
    d3 = [1.0, 0.0, 0.0]
    patch1_ptr = geometry_factory.CreateParallelepiped(p1, d1, d2, d3)
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1
    patch1.LayerIndex = 1

    # create patch 2
    p3 = [1.0 + dist, 0.0, 0.0]
    p4 = [2.0 + dist, 1.0, 1.0]
    patch2_ptr = geometry_factory.CreateSlab(p3, p4)
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2
    patch2.LayerIndex = 1

    if configuration == 1:

        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W0, patch2, BoundarySide3D.U0, False, BoundaryDirection.Forward, BoundaryDirection.Forward)

    elif configuration == 2:

        bsplines_patch_util.Reverse(patch1, 1)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W1, patch2, BoundarySide3D.U0, False, BoundaryDirection.Forward, BoundaryDirection.Reversed)

    elif configuration == 3:

        bsplines_patch_util.Reverse(patch1, 0)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W1, patch2, BoundarySide3D.U0, False, BoundaryDirection.Reversed, BoundaryDirection.Forward)

    elif configuration == 4:

        bsplines_patch_util.Reverse(patch1, 0)
        bsplines_patch_util.Reverse(patch1, 1)
        bsplines_patch_util.Reverse(patch1, 2)

        ######create multipatch
        mpatch = MultiPatch3D([patch1_ptr, patch2_ptr])
        bsplines_patch_util.MakeInterface(patch1, BoundarySide3D.W0, patch2, BoundarySide3D.U0, False, BoundaryDirection.Reversed, BoundaryDirection.Reversed)

    else:
        raise Exception("Unknown configuration %d" % (configuration))

    return mpatch

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
            mpatch = CreateMultiPatch(setting=setting, configuration = configuration, dist = 0.0)
            mpatch.Enumerate()
            # mpatch_export.Export(mpatch, "two_cubes.m")
            multipatch_util.CheckInterfaces(mpatch, True) # check the indices on patch interface

            # bpatch_ptr = mpatch[2].GetReference().ConstructBoundaryPatch(BoundarySide3D.U0)
            # func_indices = bpatch_ptr.GetReference().FESpace().FunctionIndices()
            # # print(func_indices)
            # assert(is_same(func_indices, ref_interface_indices[setting][configuration]))

    print("Test passed")

if __name__ == '__main__':
    test()
