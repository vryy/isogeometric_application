##################################################################
######################## include.py   ############################
##################################################################
##### ekate - Enhanced KRATOS for Advanced Tunnel Enineering #####
##### copyright by CIMNE, Barcelona, Spain                   #####
#####          and Institute for Structural Mechanics, RUB   #####
##### all rights reserved                                    #####
##################################################################
##################################################################
##################################################################
##################################################################
import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

import geometry_factory

mpatch_export3 = MultiNURBSPatchMatlabExporter()

def main():
    fname = os.path.basename(__file__)
    name = os.path.splitext(fname)[0]

    points1 = [[0.0, 0.0, 0.0], [1.0, 0.1, 0.0], [2.0, 0.1, 0.0], [3.0, 0.1, 0.0]]
    curve1_ptr = geometry_factory.CreateCurve(points1, 2)
    curve1 = curve1_ptr.GetReference()
    curve1.Id = 1
    mpatch_export3.Export(curve1_ptr, name + "_curve1.m")

    points2 = [[0.0, 3.0, 0.0], [0.5, 1.0, 0.0], [1.0, 0.1, 0.0], [2.0, -1.0, 0.0], [3.0, -3.0, 0.0]]
    curve2_ptr = geometry_factory.CreateCurve(points2, 2)
    curve2 = curve2_ptr.GetReference()
    curve2.Id = 2
    mpatch_export3.Export(curve2_ptr, name + "_curve2.m")

    intersect_util = IsogeometricIntersectionUtility()
    [stat, points] = intersect_util.ComputeIntersection(0.5, 0.5, curve1_ptr, curve2_ptr, 30, 1.0e-6, 0)

    print("stat:", stat) # 0 is OK
    print("points:", points)

    gf1 = curve1.GridFunction(CONTROL_POINT)
    P1 = gf1.GetValue([points[0]])
    print("P1:", str(P1))

    gf2 = curve2.GridFunction(CONTROL_POINT)
    P2 = gf2.GetValue([points[1]])
    print("P2:", str(P2))

if __name__ == "__main__":
    main()

