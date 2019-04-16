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
from KratosMultiphysics.StructuralApplication import *
kernel = Kernel()   #defining kernel

def main():
    model_part = ModelPart("test")

    p1 = Vector(3)
    p1[0] = 0.0
    p1[1] = 0.0
    p1[2] = 0.0

    p2 = Vector(3)
    p2[0] = 1.0
    p2[1] = 0.0
    p2[2] = 0.0

    p3 = Vector(3)
    p3[0] = 1.0
    p3[1] = 1.0
    p3[2] = 0.0

    p4 = Vector(3)
    p4[0] = 0.0
    p4[1] = 1.0
    p4[2] = 0.0

    points = [p1, p2, p3, p4]

    center = Vector(3)
    center[0] = 0.0
    center[1] = 0.0
    center[2] = 0.0

    normal = Vector(3)
    normal[0] = 0.0
    normal[1] = 0.0
    normal[2] = 1.0

    t1 = Vector(3)
    t1[0] = 1.0
    t1[1] = 0.0
    t1[2] = 0.0

    t2 = Vector(3)
    t2[0] = 0.0
    t2[1] = 1.0
    t2[2] = 0.0

    post_util = IsogeometricPostUtility()
    prop = model_part.Properties[1]
    nrefine = 3
    last_node_id = 0
    last_cond_id = 0
    [nodes, conds] = post_util.CreateConditions(points, center, normal, t1, t2, model_part, "Face3D3N", last_node_id, last_cond_id, prop, nrefine)
    for node in nodes:
        print(node)
    for cond in conds:
        model_part.AddCondition(cond)

    print(model_part)

    #######WRITE TO GID
    time = 0.0
    write_deformed_flag = WriteDeformedMeshFlag.WriteUndeformed
    write_elements = WriteConditionsFlag.WriteConditions
    #write_elements = WriteConditionsFlag.WriteElementsOnly
    post_mode = GiDPostMode.GiD_PostBinary
    multi_file_flag = MultiFileFlag.MultipleFiles
    gid_io = StructuralGidIO("triangle", post_mode, multi_file_flag, write_deformed_flag, write_elements)
    gid_io.InitializeMesh( time )
    post_mesh = model_part.GetMesh()
    gid_io.WriteMesh( post_mesh )
    print("mesh written...")
    gid_io.FinalizeMesh()

if __name__ == "__main__":
    main()

