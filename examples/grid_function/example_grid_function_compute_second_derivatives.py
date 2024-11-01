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

    d2P = grid_func.GetSecondDerivative([0.1, 0.1])
    print(d2P[0])
    print(d2P[1])

if __name__ == "__main__":
    main()
