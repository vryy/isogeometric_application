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
from KratosMultiphysics.BRepApplication import *
from KratosMultiphysics.IsogeometricApplication import *
kernel = Kernel()   #defining kernel

nurbs_fespace_library = BSplinesFESpaceLibrary()
grid_lib = ControlGridLibrary()
multipatch_util = MultiPatchUtility()
multipatch_refine_util = MultiPatchRefinementUtility()
bsplines_patch_util = BSplinesPatchUtility()
mpatch_export = MultiNURBSPatchMatlabExporter()

import geometry_factory

axis = 'z'

####### create arc 1
r1 = 1.0
arc1_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], axis, r1, 0.0, 90.0)
arc1 = arc1_ptr.GetReference()

####### create arc 2
r2 = 2.0
arc2_ptr = geometry_factory.CreateSmallArc([0.0, 0.0, 0.0], axis, r2, 0.0, 90.0)
arc2 = arc2_ptr.GetReference()

# create patch 2
patch_ptr = bsplines_patch_util.CreateLoftPatch(arc1, arc2)
patch = patch_ptr.GetReference()
patch.Id = 1

######create multipatch
mpatch = MultiPatch2D()
mpatch.AddPatch(patch)
#multipatch_refine_util.DegreeElevate(patch_ptr, [1, 1])
print(mpatch)

mpatch_export.Export(mpatch, "internal_pressurized_cylinder_2d.m")


