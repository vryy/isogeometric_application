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
multipatch_refine_util = MultiPatchRefinementUtility()
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export1 = MultiNURBSPatchGeoExporter()
mpatch_export2 = MultiNURBSPatchGLVisExporter()
mpatch_export3 = MultiNURBSPatchMatlabExporter()

import example_multipatch_insert_knots

def main():
    mpatch = example_multipatch_insert_knots.CreateMultiPatch()
    multipatch_refine_util.InsertKnots(mpatch[1], [[], [0.25]])
    print("Multipatch is created")
    system_size = mpatch.Enumerate()
    print("system_size:", system_size)
    print(mpatch)

    mpatch_export3.Export(mpatch, "mpatch.m")

    patch1_ptr = mpatch[1]
    patch1 = patch1_ptr.GetReference()
    bsplines_patch_util.Reverse(patch1, 1)

    print("Multipatch after reverse:")
    print(mpatch)

    mpatch_export3.Export(mpatch, "mpatch_reversed.m")

if __name__ == "__main__":
    main()

#################RESULTS#####################(Validated with Matlab)
##>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<


