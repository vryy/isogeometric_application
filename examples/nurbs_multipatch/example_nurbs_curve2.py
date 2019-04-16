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

p = 2
nurbs_fespace_library = BSplinesFESpaceLibrary()
fes = nurbs_fespace_library.CreateLinearFESpace(p)

ctrl_grid = StructuredControlPointGrid1D(3)
ctrl_grid.SetValue(0, ControlPoint(0.0, 0.0, 0.0, 1.0))
ctrl_grid.SetValue(1, ControlPoint(0.5, 0.0, 0.0, 1.0))
ctrl_grid.SetValue(2, ControlPoint(1.0, 1.0, 0.0, 1.0))

multipatch_util = MultiPatchUtility()
patch_ptr = multipatch_util.CreatePatchPointer(1, fes)
patch = patch_ptr.GetReference()
patch.CreateControlPointGridFunction(ctrl_grid)

print("###H-REFINEMENT###")

mpatch = MultiPatch1D()
mpatch.AddPatch(patch_ptr)

multipatch_refine_util = MultiPatchRefinementUtility()

#multipatch_refine_util.DegreeElevate(mpatch[1],[1])

multipatch_refine_util.InsertKnots(mpatch[1], [[0.5]]) # -> C1

print("###EXPORT###")
mpatch_export = MultiNURBSPatchMatlabExporter()
mpatch_export.Export(mpatch, "curve2.m")


