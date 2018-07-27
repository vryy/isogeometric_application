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

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()

import geometry_factory

def CreateMultiPatch():
    ## create arc 1
    r1 = 1.0
    arc1_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r1, 0.0, 90.0) # put 'y' will give interesting ring
    arc1 = arc1_ptr.GetReference()
    arc1.Id = 1

    ## create arc 2
    r2 = 2.0
    arc2_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], 'z', r2, 0.0, 90.0)
    arc2 = arc2_ptr.GetReference()
    arc2.Id = 2

    ## create ring patch by connect the two arcs
    ring_patch_ptr = bsplines_patch_util.CreateConnectedPatch(arc2, arc1)
    ring_patch = ring_patch_ptr.GetReference()
    ring_patch.Id = 1

    ######create multipatch
    mpatch = MultiPatch2D()
    mpatch.AddPatch(ring_patch_ptr)

    return mpatch

def main():
    mpatch = CreateMultiPatch()
    mpatch.Enumerate()
    print(mpatch)
    mpatch_export.Export(mpatch, "ring.m")

if __name__ == "__main__":
    main()

