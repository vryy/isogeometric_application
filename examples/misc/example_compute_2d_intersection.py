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

def main():

    p1 = Vector(3)
    p1[0] = 0.0
    p1[1] = 0.0
    p1[2] = 0.0

    p2 = Vector(3)
    p2[0] = 1.0
    p2[1] = 1.0
    p2[2] = 0.0

    p3 = Vector(3)
    p3[0] = 1.0
    p3[1] = 0.0
    p3[2] = 0.0

    p4 = Vector(3)
    p4[0] = 0.0
    p4[1] = 1.0
    p4[2] = 0.0

    math_utils = IsogeometricMathUtils()
    output = math_utils.ComputeIntersection2D(p1, p2, p3, p4, 1e-13)
    assert(output[0] == 0)
    assert(abs(output[1][0] - 0.5) < 1e-13)
    assert(abs(output[1][1] - 0.5) < 1e-13)

    p3[0] = 2.0
    p3[1] = 0.0
    p3[2] = 0.0

    p4[0] = 2.0
    p4[1] = 3.0
    p4[2] = 0.0

    output = math_utils.ComputeIntersection2D(p1, p2, p3, p4, 1e-13)
    assert(output[0] == 2)
    assert(abs(output[1][0] - 2.0) < 1e-13)
    assert(abs(output[1][1] - 2.0) < 1e-13)

if __name__ == "__main__":
    main()

