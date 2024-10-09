##################################################################
##################################################################
import sys
import os
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

import geometry_factory

mpatch_import = MultiNURBSPatchGeoImporter3D()
mpatch_export3 = MultiNURBSPatchMatlabExporter()

def main():
    fname = os.path.basename(__file__)
    name = os.path.splitext(fname)[0]

    patch_ptr = mpatch_import.ImportSingle("../../tests/geo_ring3d.txt")
    patch = patch_ptr.GetReference()

    intersect_util = IsogeometricIntersectionUtility()
    A = 1.0
    B = 0.0
    C = 0.0
    D = -0.5

    # firstly, check the intersection
    stat = intersect_util.CheckIntersection(patch_ptr, A, B, C, D)
    print("intersection info: " + str(stat))

    # secondly, compute the intersection
    [stats, points] = intersect_util.ComputeIntersection([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], patch_ptr, A, B, C, D, 30, 1.0e-6)
    print("stat:", stats) # 0 is OK
    print("point:", points)

    gf = patch.GridFunction(CONTROL_POINT)
    for i in range(0, len(stats)):
        if stats[i] == 0:
            print("point:", points[i])
            P = gf.GetValue(points[i])
            print("P:", str(P))
            print("P.X:", P.X())


if __name__ == "__main__":
    main()

