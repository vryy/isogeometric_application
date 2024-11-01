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
from KratosMultiphysics.LayerApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.IsogeometricApplication import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.IsogeometricStructuralApplication import *
kernel = Kernel()   #defining kernel

import geometry_factory
import model_iga_include

multipatch_util = MultiPatchUtility()

def main():

    patch1_ptr = geometry_factory.CreateSlab([0.0, 0.0, 0.0], [1.0, 1.0,1.0])
    patch1 = patch1_ptr.GetReference()
    patch1.Id = 1
    patch1.LayerIndex = 1
    patch2_ptr = geometry_factory.CreateSlab([1.0, 0.0, 0.0], [2.0, 1.0,1.0])
    patch2 = patch2_ptr.GetReference()
    patch2.Id = 2
    patch2.LayerIndex = 2

    mpatch = MultiPatch3D()
    mpatch.AddPatch(patch1)
    mpatch.AddPatch(patch2)

    mpatch.Enumerate()
    print(mpatch)

    post_params = {}
    post_params['name'] = "slab"
    post_params['variables list'] = []
    # model_iga_include.PostMultiPatch(mpatch, 3, 0.0, post_params)

    cmmp = ConformingMultipatchLagrangeModelPart3D(mpatch)
    cmmp.SetEchoLevel(3)
    cmmp.SetUniformSampling(1, 0, 2)
    cmmp.SetUniformSampling(1, 1, 2)
    cmmp.SetUniformSampling(1, 2, 2)
    cmmp.SetUniformSampling(2, 0, 2)
    cmmp.SetUniformSampling(2, 1, 2)
    cmmp.SetUniformSampling(2, 2, 2)
    cmmp.BeginModelPart()
    mp = cmmp.GetModelPart()
    cmmp.CreateNodes()
    last_elem_id = multipatch_util.GetLastElementId(mp)
    cmmp.AddElements(patch1, "KinematicLinear3D8N", last_elem_id+1, mp.Properties[1])
    last_elem_id = multipatch_util.GetLastElementId(mp)
    cmmp.AddElements(patch2, "KinematicLinear3D8N", last_elem_id+1, mp.Properties[2])
    last_cond_id = multipatch_util.GetLastConditionId(mp)
    cmmp.AddConditions(patch1, BoundarySide3D.U0, "FaceForce3D4N", last_cond_id+1, mp.Properties[1])
    last_cond_id = multipatch_util.GetLastConditionId(mp)
    cmmp.AddConditions(patch1, 0, 0.5, "FaceForce3D4N", last_cond_id+1, mp.Properties[1])
    last_cond_id = multipatch_util.GetLastConditionId(mp)
    cmmp.AddConditions(patch2, BoundarySide3D.U1, "FaceForce3D4N", last_cond_id+1, mp.Properties[2])
    cmmp.EndModelPart()

    print(mp)

    model_iga_include.WriteGiD(mp, 0.0, post_params)

if __name__ == "__main__":
    main()

