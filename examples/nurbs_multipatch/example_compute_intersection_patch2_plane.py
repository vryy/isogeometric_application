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

mpatch_import = MultiNURBSPatchGeoImporter2D()
mpatch_export3 = MultiNURBSPatchMatlabExporter()

def main():
    fname = os.path.basename(__file__)
    name = os.path.splitext(fname)[0]

    patch_ptr = mpatch_import.ImportSingle("../../tests/geo_ring.txt")
    patch = patch_ptr.GetReference()

    intersect_util = IsogeometricIntersectionUtility()
    A = 1.0
    B = 0.0
    C = 0.0
    D = -0.5
    [stats, points] = intersect_util.ComputeIntersection([0.5, 0.5, 0.5, 0.5], patch_ptr, A, B, C, D, 30, 1.0e-6)
    print("stat:", stats) # 0 is OK
    print("point:", points)

    gf = patch.GridFunction(CONTROL_POINT)
    P1 = gf.GetValue(points[0])
    print("P1:", str(P1))
    print("P1.X:", str(P1.X()))
    P2 = gf.GetValue(points[1])
    print("P2:", str(P2))
    print("P2.X:", str(P2.X()))

if __name__ == "__main__":
    main()

