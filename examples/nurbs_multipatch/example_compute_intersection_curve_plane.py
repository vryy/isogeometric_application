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

    points = [[0.0, 3.0, 0.0], [0.5, 1.0, 0.0], [1.0, 0.1, 0.0], [2.0, -1.0, 0.0], [3.0, -3.0, 0.0]]
    curve_ptr = geometry_factory.CreateCurve(points, 2)
    curve = curve_ptr.GetReference()
    curve.Id = 1
    mpatch_export3.Export(curve_ptr, name + "_curve2.m")

    intersect_util = IsogeometricIntersectionUtility()
    A = 0.0
    B = 0.0
    C = 1.0
    D = -1.0
    [stat, point] = intersect_util.ComputeIntersection(0.5, curve_ptr, A, B, C, D, 30, 1.0e-6)

    print("stat:", stat)
    print("point:", point)

if __name__ == "__main__":
    main()

