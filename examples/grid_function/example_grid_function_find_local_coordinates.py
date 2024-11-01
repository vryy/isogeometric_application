##################################################################
##################################################################
import sys
import os
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

def main():
    patch_util = BSplinesPatchUtility()
    mpatch_export = MultiNURBSPatchMatlabExporter()

    patches_ptr = patch_util.CreatePatchFromGeo("../../tests/geo_ring.txt")
    patch = patches_ptr[0].GetReference()
    # print(patch)

    grid_func = patch.GridFunction(CONTROL_POINT_COORDINATES)
    print(grid_func)
    P = Vector(3)
    P[0] = 1.0
    P[1] = 1.0
    P[2] = 0.0
    xi0 = Vector(3)
    xi0[0] = 0.5
    xi0[1] = 0.5
    xi0[2] = 0.0
    [stat, xi] = grid_func.LocalCoordinates(P, xi0)
    print("stat:", stat)
    print("xi:", str(xi))
    P1 = grid_func.GetValue([xi[0], xi[1]])
    print("P1:", str(P1))

    assert(abs(P1[0] - 1.0) < 1e-13)
    assert(abs(P1[1] - 1.0) < 1e-13)

if __name__ == "__main__":
    main()
